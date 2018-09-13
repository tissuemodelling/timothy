/* **************************************************************************
 * Timothy - Tissue Modelling Framework
 * Copyright (C) 2014-2018 Maciej Cytowski
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * *************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "global.h"
#include "octree.h"
#include "utils.h"

/*! \file octree.c
 *  \brief contains functions that build octree and functions
 *         that provide a convenient way of operating on octree.
 */

/*!
 * This function creates an empty octree node in a given position.
 */
void octemptynode(int64_t father,int level,unsigned int xbit,unsigned int ybit,unsigned int zbit,systeminfo_t systeminfo, cellsinfo_t *cellsinfo)
{
        int i;
        unsigned int bit;
        unsigned int len;
        octnode_t *octree=cellsinfo->octree;
        int64_t octsize=cellsinfo->octsize;

        if(father==-1) {
                octree[octsize].xcode=0;
                octree[octsize].ycode=0;
                octree[octsize].zcode=0;
                octree[octsize].xlimit=1<<level;
                octree[octsize].ylimit=1<<level;
                octree[octsize].zlimit=1<<level;
        } else {
                bit=1<<level;
                octree[octsize].xcode=octree[father].xcode|xbit;
                octree[octsize].ycode=octree[father].ycode|ybit;
                octree[octsize].zcode=octree[father].zcode|zbit;
                len=(octree[father].xlimit-octree[father].xcode)/2;
                octree[octsize].xlimit=octree[octsize].xcode+len;
                octree[octsize].ylimit=octree[octsize].ycode+len;
                octree[octsize].zlimit=octree[octsize].zcode+len;
        }
        octree[octsize].father=father;
        for(i=0; i<8; i++)
                octree[octsize].child[i]=-1;
        octree[octsize].data=-1;
        octree[octsize].level=level;
        if(cellsinfo->octsize+1==cellsinfo->octmaxsize) {
                cellsinfo->octmaxsize+=64;
                if(!(cellsinfo->octree=(octnode_t*)realloc(cellsinfo->octree,sizeof(octnode_t)*cellsinfo->octmaxsize))) {
                        terminate(systeminfo,"cannot reallocate cellsinfo->octree", __FILE__, __LINE__);
                }
        }
        cellsinfo->octsize++;

        return;
}

/*!
 * This function inserts a given cell c into a proper position in octree.
 */
void octinsertcell(int64_t c,uint3dv_t *loccode,systeminfo_t systeminfo, cellsinfo_t *cellsinfo)
{
        int64_t octIdx=0;
        int64_t father,child;
        unsigned int level;
        unsigned int childBranchBit,childIndex;
        unsigned int xbit,ybit,zbit;
        int inserted=0;
        octnode_t *octree=cellsinfo->octree;

        level=ROOT_LEVEL-1;
        while(!inserted) {
                if(octree[octIdx].data==-1) { /* insert cell into node */
                        octree[octIdx].data=c;
                        inserted=1;
                } else {
                        if(octree[octIdx].data==-2) { /* this is not a leaf */
                                childBranchBit = 1 << (level);
                                xbit=(loccode[c].x) & childBranchBit;
                                ybit=(loccode[c].y) & childBranchBit;
                                zbit=(loccode[c].z) & childBranchBit;
                                childIndex = (   ((xbit) >> (level))
                                                 + ((ybit) >> (level-1))
                                                 + ((zbit) >> (level-2)) );
                                level--;
                                father=octIdx;
                                child = octree[octIdx].child[childIndex];
                                if(child==-1) {
                                        octree[octIdx].child[childIndex]=cellsinfo->octsize;
                                        child=cellsinfo->octsize;
                                        octemptynode(father,level,xbit,ybit,zbit,systeminfo,cellsinfo);
                                }
                                octIdx=child;
                        } else { /* occupied leaf */
                                int64_t oc;
                                oc=octree[octIdx].data;
                                childBranchBit = 1 << (level);
                                xbit=(loccode[oc].x) & childBranchBit;
                                ybit=(loccode[oc].y) & childBranchBit;
                                zbit=(loccode[oc].z) & childBranchBit;
                                childIndex = (   ((xbit) >> (level))
                                                 + ((ybit) >> (level-1))
                                                 + ((zbit) >> (level-2)) );
                                level--;
                                father=octIdx;
                                octree[octIdx].child[childIndex]=cellsinfo->octsize;
                                child=cellsinfo->octsize;
                                octemptynode(father,level,xbit,ybit,zbit,systeminfo,cellsinfo);
                                octree[child].data=octree[father].data;
                                octree[father].data=-2;
                                octIdx=father;
                                level++;
                        }
                }
        }
        return;
}

/*!
 * This is a driving function for creating an octree.
 */
void octbuild(systeminfo_t systeminfo, cellsinfo_t *cellsinfo,celltype_t* celltype)
{
        int i,c;
        double3dv_t bMin,bMax;
        double epsilon=0.01;
        uint3dv_t *loccode;

        if(cellsinfo->localcount.n==0)
                return;

        bMin.x= DBL_MAX;
        bMin.y= DBL_MAX;
        bMin.z= DBL_MAX;
        bMax.x=-DBL_MAX;
        bMax.y=-DBL_MAX;
        bMax.z=-DBL_MAX;
        for(i=0; i<cellsinfo->localcount.n; i++) {
                double x,y,z;
                double e;
                x=cellsinfo->cells[i].x;
                y=cellsinfo->cells[i].y;
                z=cellsinfo->cells[i].z;
                e=celltype[cellsinfo->cells[i].ctype].h+epsilon;
                bMin.x=(x-e<bMin.x ? x-e : bMin.x);
                bMax.x=(x+e>bMax.x ? x+e : bMax.x);
                bMin.y=(y-e<bMin.y ? y-e : bMin.y);
                bMax.y=(y+e>bMax.y ? y+e : bMax.y);
                bMin.z=(z-e<bMin.z ? z-e : bMin.z);
                bMax.z=(z+e>bMax.z ? z+e : bMax.z);
        }
        affScale=bMax.x-bMin.x;
        affScale=(affScale>bMax.y-bMin.y ? affScale : bMax.y-bMin.y);
        affScale=(affScale>bMax.z-bMin.z ? affScale : bMax.z-bMin.z);
        affShift.x=bMin.x;
        affShift.y=bMin.y;
        affShift.z=bMin.z;

        /* each cell coordinate will be shifted by affShift and scaled by affScale */

        if(!(loccode=(uint3dv_t*) malloc(sizeof(uint3dv_t)*cellsinfo->localcount.n)))
                terminate(systeminfo,"cannot allocate loccode", __FILE__, __LINE__);
        //#pragma omp parallel for
        for(c=0; c<cellsinfo->localcount.n; c++) {
                loccode[c].x=(unsigned int)( ((cellsinfo->cells[c].x-affShift.x)/affScale) * MAXVALUE );
                loccode[c].y=(unsigned int)( ((cellsinfo->cells[c].y-affShift.y)/affScale) * MAXVALUE );
                loccode[c].z=(unsigned int)( ((cellsinfo->cells[c].z-affShift.z)/affScale) * MAXVALUE );
        }
        /* memory space required to store octree (to be reviewed again!) */
        cellsinfo->octmaxsize=cellsinfo->localcount.n*64;
        if(!(cellsinfo->octree=(octnode_t*) malloc(sizeof(octnode_t)*cellsinfo->octmaxsize)))
                terminate(systeminfo,"cannot allocate cellsinfo->octree", __FILE__, __LINE__);

        cellsinfo->octsize=0;
        octemptynode(-1,ROOT_LEVEL,0,0,0,systeminfo,cellsinfo);

        for(c=0; c<cellsinfo->localcount.n; c++) {
                octinsertcell(c,loccode,systeminfo,cellsinfo);
        }

        /* a security check for space available in the octree buffer should be implemented */
        free(loccode);
        return;
}

/*!
 * This function computes bounding box of a remote cell to be used for neighbour searching.
 * Difference between remote and local versions is the box croping which is applied
 * to some remote cells.
 */
void octcomputeboxr(int64_t c,uint3dv_t *minloccode,uint3dv_t *maxloccode,cellcommdata_t cellcommdata,celltype_t* celltype)
{
        double3dv_t mincor,maxcor;
        double h=celltype[cellcommdata.recvcellindata[c].ctype].h;
        /* compute corners */
        mincor.x=(cellcommdata.recvcellindata[c].x-h-affShift.x)/affScale;
        mincor.y=(cellcommdata.recvcellindata[c].y-h-affShift.y)/affScale;
        mincor.z=(cellcommdata.recvcellindata[c].z-h-affShift.z)/affScale;
        maxcor.x=(cellcommdata.recvcellindata[c].x+h-affShift.x)/affScale;
        maxcor.y=(cellcommdata.recvcellindata[c].y+h-affShift.y)/affScale;
        maxcor.z=(cellcommdata.recvcellindata[c].z+h-affShift.z)/affScale;
        /* for remote cells - box crop */
        mincor.x=(mincor.x<0.0 ? 0.0 : mincor.x);
        mincor.y=(mincor.y<0.0 ? 0.0 : mincor.y);
        mincor.z=(mincor.z<0.0 ? 0.0 : mincor.z);
        maxcor.x=(maxcor.x>1.0 ? 1.0 : maxcor.x);
        maxcor.y=(maxcor.y>1.0 ? 1.0 : maxcor.y);
        maxcor.z=(maxcor.z>1.0 ? 1.0 : maxcor.z);
        /* compute location codes of corners */
        minloccode[0].x=(unsigned int)(mincor.x*MAXVALUE);
        minloccode[0].y=(unsigned int)(mincor.y*MAXVALUE);
        minloccode[0].z=(unsigned int)(mincor.z*MAXVALUE);
        maxloccode[0].x=(unsigned int)(maxcor.x*MAXVALUE);
        maxloccode[0].y=(unsigned int)(maxcor.y*MAXVALUE);
        maxloccode[0].z=(unsigned int)(maxcor.z*MAXVALUE);

        return;
}

/*!
 * This function computes bounding box of a local cell to be used for neighbour searching.
 */
void octcomputebox(int64_t c,uint3dv_t *minloccode,uint3dv_t *maxloccode,cellsinfo_t cellsinfo,celltype_t* celltype)
{
        double3dv_t mincor,maxcor;
        double h=celltype[cellsinfo.cells[c].ctype].h;
        /* compute corners */
        mincor.x=(cellsinfo.cells[c].x-h-affShift.x)/affScale;
        mincor.y=(cellsinfo.cells[c].y-h-affShift.y)/affScale;
        mincor.z=(cellsinfo.cells[c].z-h-affShift.z)/affScale;
        maxcor.x=(cellsinfo.cells[c].x+h-affShift.x)/affScale;
        maxcor.y=(cellsinfo.cells[c].y+h-affShift.y)/affScale;
        maxcor.z=(cellsinfo.cells[c].z+h-affShift.z)/affScale;
        /* compute location codes of corners */
        minloccode[0].x=(unsigned int)(mincor.x*MAXVALUE);
        minloccode[0].y=(unsigned int)(mincor.y*MAXVALUE);
        minloccode[0].z=(unsigned int)(mincor.z*MAXVALUE);
        maxloccode[0].x=(unsigned int)(maxcor.x*MAXVALUE);
        maxloccode[0].y=(unsigned int)(maxcor.y*MAXVALUE);
        maxloccode[0].z=(unsigned int)(maxcor.z*MAXVALUE);

        return;
}

/*!
 * This function traverses the tree to a given level and returns node index.
 */
int octtraversetolevel(unsigned int level,unsigned int xloc,unsigned int yloc,unsigned int zloc,unsigned int lmin,cellsinfo_t cellsinfo)
{
        int cellidx=0;
        int parent=0;
        unsigned int n;
        unsigned int childbranchbit,childindex;
        n=(level)-(lmin)+1;
        while(n--) {
                childbranchbit = 1 << (level);
                childindex = (   (((xloc) & childbranchbit) >> (level))
                                 + (((yloc) & childbranchbit) >> (level-1))
                                 + (((zloc) & childbranchbit) >> ((level-2))) );
                level--;
                parent=cellidx;
                cellidx=cellsinfo.octree[cellidx].child[childindex];
                if(cellsinfo.octree[cellidx].data != -2) break; /* a leaf */
        }
        if(cellidx==-1) cellidx=parent;
        return cellidx;
}

/*!
 * This function determines the level of the smallest cell containing a given region.
 */
int octlocateregion(uint3dv_t minloccode,uint3dv_t maxloccode,cellsinfo_t cellsinfo)
{
        unsigned int l1,l2,lmin;
        unsigned int level;
        int cellidx;
        uint3dv_t locdiff;
        /* XOR of location codes */
        locdiff.x=minloccode.x ^ maxloccode.x;
        locdiff.y=minloccode.y ^ maxloccode.y;
        locdiff.z=minloccode.z ^ maxloccode.z;
        /* Determining the level of the smallest cell containing the region */
        l1=ROOT_LEVEL;
        l2=ROOT_LEVEL;
        lmin=ROOT_LEVEL;
        while(!(locdiff.x & (1<<l1)) && l1) l1--;
        while(!(locdiff.y & (1<<l2)) && (l2>l1)) l2--;
        while(!(locdiff.z & (1<<lmin)) && (lmin>l2)) lmin--;
        lmin++;
        level=ROOT_LEVEL-1;
        cellidx=octtraversetolevel(level,minloccode.x,minloccode.y,minloccode.z,lmin,cellsinfo);
        return cellidx;
}

/*!
 * This function deallocates memory used by octree.
 */
void octfree(cellsinfo_t *cellsinfo)
{
        if(cellsinfo->localcount.n==0) return;
        free(cellsinfo->octree);
        return;
        //  free(locCode);
}

/*!
 * This function initializes heap for tree traversal.
 */
void octheapinit(systeminfo_t systeminfo,octheap_t *ttheap)
{
        ttheap->size=64;
        ttheap->count=0;
        if(!(ttheap->data=malloc(sizeof(int)*ttheap->size)))
                terminate(systeminfo,"cannot allocate ttheap->data", __FILE__, __LINE__);
        return;
}

/*!
 * This function adds an element to the heap.
 */
void octheappush(systeminfo_t systeminfo, octheap_t *ttheap,int idx)
{
        if(ttheap->count==ttheap->size) {
                ttheap->size+=64;
                if(!(ttheap->data=realloc(ttheap->data,sizeof(int)*ttheap->size)))
                        terminate(systeminfo,"cannot reallocate ttheap->data", __FILE__, __LINE__);
        }
        ttheap->data[ttheap->count]=idx;
        ttheap->count+=1;
        return;
}

/*!
 * This function removes an element from the heap.
 */
int octheappop(octheap_t *ttheap)
{
        ttheap->count-=1;
        return ttheap->data[ttheap->count];
}

/*!
 * This functione deallocates memory used by the heap.
 */
void octheapfree(octheap_t *ttheap)
{
        free(ttheap->data);
        return;
}
