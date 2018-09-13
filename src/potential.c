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
#include "inline.h"

/*! \file potential.c
 *  \brief contains functions that compute the potential
 */

/*!
 * This function computes potential for two neighbour cells.
 */
double potential(int dimension,celldata_t* c1,celldata_t* c2,celltype_t celltype)
{
        double dist=0.0;
        double h, h2, h3;
        double xc, D;
        double poisson = 0.33;
        double young;
        double csize;

        dist = sqrt( (c1->x-c2->x)*(c1->x-c2->x) + (c1->y-c2->y)*(c1->y-c2->y) + (c1->z-c2->z)*(c1->z-c2->z));
        c1->mindist=(dist<c1->mindist ? dist : c1->mindist);

        h = celltype.h;
        h2 = celltype.h2;
        h3 = celltype.h3;
        csize = celltype.size;

        if (dist <= h) {
                double r01, r02;
                double area;
                double sc=1.0;
                sc=(dimension==2 ? h2 : h3);

                c1->density+=sc * (c2->size / csize) * sphkernel(dimension,dist,h,h2,h3);
                xc = c1->size + c2->size - dist;
                if (xc <= 0.0)
                        return 0.0;
                D = 0.75 * ((1.0 - poisson * poisson) / c1->young + (1.0 - poisson * poisson / c2->young));

                /* adhesion */
                r01 = (c1->size * c1->size - c2->size * c2->size + dist * dist) / (2 * dist);
                r02 = dist - r01;
                area = M_PI * ((c1->size * c1->size * (c1->size - r01) - (c1->size * c1->size * c1->size - r01 * r01 * r01) / 3) +
                               (c2->size * c2->size * (c2->size - r02) - (c2->size * c2->size * c2->size - r02 * r02 * r02) / 3));

                /* compute potential */
                return (2.0 * pow(xc, 5 / 2) / (5.0 * D)) * sqrt((c1->size * c2->size) / (c1->size + c2->size)) + area * 0.1;
        } else
                return 0.0;
}




/*!
 * This function implements tree walk algorithm for each local cell.
 * Function ccPot(..) is called for each pair of neighbours.
 */
void computepotential(systeminfo_t systeminfo, cellsinfo_t *cellsinfo,celltype_t* celltype,cellcommdata_t cellcommdata)
{
        if(cellsinfo->localcount.n<=1) return;
        #pragma omp parallel
        {
                int p;
                #pragma omp for schedule(dynamic,64)
                for (p = 0; p < cellsinfo->localcount.n; p++) {
                        int64_t cellIdx;
                        int newIdx;
                        int s;
                        int octIndex;
                        uint3dv_t minLocCode,maxLocCode;
                        octheap_t octh;

                        octheapinit(systeminfo,&octh);

                        cellsinfo->cells[p].density = 0.0;
                        cellsinfo->cells[p].v = 0.0;
                        cellsinfo->cells[p].mindist = DBL_MAX;

                        octcomputebox(p,&minLocCode,&maxLocCode,*cellsinfo,celltype);
                        octIndex=octlocateregion(minLocCode,maxLocCode,*cellsinfo);

                        for(s=0; s<8; s++) {
                                int idx;
                                idx=cellsinfo->octree[octIndex].child[s];
                                if(idx!=-1 && octnodeintersection(idx,minLocCode,maxLocCode,*cellsinfo))
                                        octheappush(systeminfo,&octh,idx);
                        }

                        while(octh.count>0) {
                                int idx;
                                idx=octheappop(&octh);
                                cellIdx=cellsinfo->octree[idx].data;
                                if(cellIdx>=0) {
                                        if(cellIdx!=p) {
                                                cellsinfo->cells[p].v = potential(cellsinfo->dimension,&(cellsinfo->cells[p]),&(cellsinfo->cells[cellIdx]),celltype[cellsinfo->cells[p].ctype]);
                                        }
                                } else {
                                        for(s=0; s<8; s++) {
                                                newIdx=cellsinfo->octree[idx].child[s];
                                                if(newIdx!=-1 && octnodeintersection(newIdx,minLocCode,maxLocCode,*cellsinfo))
                                                        octheappush(systeminfo,&octh,newIdx);
                                        }
                                }
                        }
                        octheapfree(&octh);
                }
        }
        return;
}

/*!
 * This function implements tree walk algorithm for each remote cell.
 * Function ccPot(..) is called for each pair of neighbours.
 */
void computeremotepotential(systeminfo_t systeminfo, cellsinfo_t *cellsinfo,celltype_t* celltype,cellcommdata_t cellcommdata)
{
        if(cellcommdata.numimp<=0) return;
        #pragma omp parallel
        {
                int rp;
                #pragma omp for schedule(dynamic,64)
                for (rp = 0; rp < cellcommdata.numimp; rp++) {
                        int64_t cellIdx;
                        int newIdx;
                        int s;
                        double v;
                        uint3dv_t minLocCode,maxLocCode;
                        octheap_t octh;
                        celldata_t cell;
                        cell.x=cellcommdata.recvcellindata[rp].x;
                        cell.y=cellcommdata.recvcellindata[rp].y;
                        cell.z=cellcommdata.recvcellindata[rp].z;
                        cell.size=cellcommdata.recvcellindata[rp].size;
                        cell.young=cellcommdata.recvcellindata[rp].young;
                        cell.ctype=cellcommdata.recvcellindata[rp].ctype;
                        octheapinit(systeminfo,&octh);
                        octcomputeboxr(rp,&minLocCode,&maxLocCode,cellcommdata,celltype);
                        octheappush(systeminfo,&octh,0);
                        while(octh.count>0) {
                                int idx;
                                idx=octheappop(&octh);
                                cellIdx=cellsinfo->octree[idx].data;
                                if(cellIdx>=0) {
                                        v=potential(cellsinfo->dimension,&(cellsinfo->cells[cellIdx]),&(cell),celltype[cellsinfo->cells[cellIdx].ctype]);
                                        #pragma omp atomic
                                        cellsinfo->cells[cellIdx].v += v;
                                } else {
                                        for(s=0; s<8; s++) {
                                                newIdx=cellsinfo->octree[idx].child[s];
                                                if(newIdx!=-1 && octnodeintersection(newIdx,minLocCode,maxLocCode,*cellsinfo))
                                                        octheappush(systeminfo,&octh,newIdx);
                                        }
                                }
                        }
                        octheapfree(&octh);
                }
        }
        return;
}

/*!
 *  This function computes potential gradient for two neighbour cells.
 */
void pgradient(int dimension,celldata_t c1,celldata_t c2,double3dv_t* f,celltype_t celltype)
{
        double h,h3,h4;
        h = celltype.h;
        h3 = celltype.h3;
        h4 = celltype.h4;
        /* compute the gradient of SPH kernel function */
        sphgradient(dimension,c1,c2,f,h,h3,h4);
        return;
}


/*!
 * This function implements tree walk algorithm for each local cell to compute potential gradient.
 * Function ccPotGrad(..) is called for each pair of neighbours.
 */
void computegradient(systeminfo_t systeminfo, cellsinfo_t *cellsinfo,celltype_t* celltype,cellcommdata_t cellcommdata)
{
        int p;
        if(cellsinfo->localcount.n<=1) return;
  #pragma omp parallel for schedule(dynamic,64)
        for(p=0; p<cellsinfo->localcount.n; p++) {
                int64_t cellIdx;
                int newIdx;
                int s;
                int octIndex;
                uint3dv_t minLocCode,maxLocCode;
                octheap_t octh;

                octheapinit(systeminfo,&octh);
                cellsinfo->forces[p].x=0.0;
                cellsinfo->forces[p].y=0.0;
                cellsinfo->forces[p].z=0.0;

                octcomputebox(p,&minLocCode,&maxLocCode,*cellsinfo,celltype);
                octIndex=octlocateregion(minLocCode,maxLocCode,*cellsinfo);

                for(s=0; s<8; s++) {
                        int idx;
                        idx=cellsinfo->octree[octIndex].child[s];
                        if(idx!=-1 && octnodeintersection(idx,minLocCode,maxLocCode,*cellsinfo))
                                octheappush(systeminfo,&octh,idx);
                }
                while(octh.count>0) {
                        int idx;
                        idx=octheappop(&octh);
                        cellIdx=cellsinfo->octree[idx].data;
                        if(cellIdx>=0) {
                                if(cellIdx!=p) {
                                        double3dv_t f;
                                        pgradient(cellsinfo->dimension,cellsinfo->cells[p],cellsinfo->cells[cellIdx],&f,celltype[cellsinfo->cells[p].ctype]);
                                        cellsinfo->forces[p].x+=f.x;
                                        cellsinfo->forces[p].y+=f.y;
                                        cellsinfo->forces[p].z+=f.z;
                                }
                        } else {
                                for(s=0; s<8; s++) {
                                        newIdx=cellsinfo->octree[idx].child[s];
                                        if(newIdx!=-1 && octnodeintersection(newIdx,minLocCode,maxLocCode,*cellsinfo))
                                                octheappush(systeminfo,&octh,newIdx);
                                }
                        }
                }
                octheapfree(&octh);
        }
        return;
}

/*!
 * This function implements tree walk algorithm for each remote cell to compute potential gradient.
 * Function ccPotGrad(..) is called for each pair of neighbours.
 */
void computeremotegradient(systeminfo_t systeminfo, cellsinfo_t *cellsinfo,celltype_t* celltype,cellcommdata_t cellcommdata)
{
        int rp;
        if(cellcommdata.numimp<=0) return;
        #pragma omp parallel for schedule(dynamic,64)
        for (rp = 0; rp < cellcommdata.numimp; rp++) {
                int64_t cellIdx;
                int newIdx;
                int s;
                uint3dv_t minLocCode,maxLocCode;
                octheap_t octh;
                celldata_t cell;
                cell.x=cellcommdata.recvcellindata[rp].x;
                cell.y=cellcommdata.recvcellindata[rp].y;
                cell.z=cellcommdata.recvcellindata[rp].z;
                cell.size=cellcommdata.recvcellindata[rp].size;
                cell.young=cellcommdata.recvcellindata[rp].young;
                cell.ctype=cellcommdata.recvcellindata[rp].ctype;
                octheapinit(systeminfo,&octh);
                octcomputeboxr(rp,&minLocCode,&maxLocCode,cellcommdata,celltype);
                octheappush(systeminfo,&octh,0);
                while(octh.count>0) {
                        int idx;
                        idx=octheappop(&octh);
                        cellIdx=cellsinfo->octree[idx].data;
                        if(cellIdx>=0) {
                                double3dv_t f;
                                pgradient(cellsinfo->dimension,cellsinfo->cells[cellIdx],cell,&f,celltype[cellsinfo->cells[cellIdx].ctype]);
                                #pragma omp atomic
                                cellsinfo->forces[cellIdx].x+=f.x;
                                #pragma omp atomic
                                cellsinfo->forces[cellIdx].y+=f.y;
                                #pragma omp atomic
                                cellsinfo->forces[cellIdx].z+=f.z;
                        } else {
                                for(s=0; s<8; s++) {
                                        newIdx=cellsinfo->octree[idx].child[s];
                                        if(newIdx!=-1 && octnodeintersection(newIdx,minLocCode,maxLocCode,*cellsinfo))
                                                octheappush(systeminfo,&octh,newIdx);
                                }
                        }
                }
                octheapfree(&octh);
        }
        return;
}
