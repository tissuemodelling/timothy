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
#include <mpi.h>
#include <math.h>
#include <inttypes.h>
#include <float.h>

#include "global.h"
#include "utils.h"
#include "inline.h"

/*! \file cells.c
 *  \brief contains functions which control current states and evolution of cells
 */

/*!
 * This function checks whether the cell p is outside the computational box.
 */
static inline int outsidethebox(cellsinfo_t *cellsinfo,int c)
{
        double x, y, z, r;

        x = cellsinfo->cells[c].x;
        y = cellsinfo->cells[c].y;
        z = cellsinfo->cells[c].z;
        r = cellsinfo->cells[c].size;

        if (x - r <= -BOXSIZEX/2.0 || x + r >= BOXSIZEX/2.0 )
                return 1;
        if (y - r <= -BOXSIZEY/2.0 || y + r >= BOXSIZEY/2.0 )
                return 1;
        if ( cellsinfo->dimension == 3 && (z - r <= -BOXSIZEZ/2.0 || z + r >= BOXSIZEZ/2.0 ) )
                return 1;

        return 0;
}

/* normal distribution (Box-Muller transform) */
void cellsrandominit(int nrandom,int ctype,systeminfo_t systeminfo,settings_t settings,celltype_t *celltype,cellsinfo_t *cellsinfo) {
        int i;
        double x1, x2, x3;
        double z1, z2, z3;
        double r1, r2;
        double l;
        double D=1.0;
        uint64_t idx;

        idx=cellsinfo->localcount.n-1;

        if (settings.dimension == 2)
                D = celltype[ctype].size*pow(8.0 * nrandom, 1.0 / 2.0);
        if (settings.dimension == 3)
                D = celltype[ctype].size*pow(8.0 * nrandom, 1.0 / 3.0);

        for (i = 0; i < nrandom; i++) {

                r2 = 1.1;
                while (r2 >= 1.0) {
                        r1 = 1.1;
                        while (r1 == 0 || r1 >= 1.0) {
                                x1 = ((float)rand()/RAND_MAX)*2.0 - 1.0;
                                x2 = ((float)rand()/RAND_MAX)*2.0 - 1.0;
                                x3 = ((float)rand()/RAND_MAX)*2.0 - 1.0;
                                r1 = x1 * x1 + x2 * x2 + x3 * x3;
                        }
                        l = sqrt(-2 * log(r1) / r1);
                        z1 = x1 * l;
                        z2 = x2 * l;
                        z3 = x3 * l;

                        r2 = z1 * z1 + z2 * z2 + z3 * z3;
                }

                cellsinfo->cells[cellsinfo->localcount.n].x=z1 * D + cellsinfo->cells[idx].x;
                cellsinfo->cells[cellsinfo->localcount.n].y=z2 * D + cellsinfo->cells[idx].y;
                if(settings.dimension>2)
                        cellsinfo->cells[cellsinfo->localcount.n].z=z3 * D + cellsinfo->cells[idx].z;
                else
                        cellsinfo->cells[cellsinfo->localcount.n].z=0.0;

                cellsinfo->cells[cellsinfo->localcount.n].ctype=ctype;
                cellsinfo->cells[cellsinfo->localcount.n].density=0.0;
                cellsinfo->cells[cellsinfo->localcount.n].size=celltype[ctype].size;
                cellsinfo->localcount.n+=1;
                cellsinfo->localcount.g0phase+=1;
                cellsinfo->localtypecount[ctype].n+=1;
                cellsinfo->localtypecount[ctype].g0phase+=1;
                cellsinfo->cells[cellsinfo->localcount.n].young=2100.0 + (float)rand() * 100.0;


                if(cellsinfo->localcount.n==settings.maxlocalcells) {
                        printf("warning: too many local cells, skipping rest of file %s\n",celltype[i].inputfile);
                        break;
                }

        }
        return;
}

/*!
 * This function implements mitosis of cells.
 */
void mitosis(systeminfo_t systeminfo, settings_t settings, celltype_t *celltype, cellsinfo_t *cellsinfo, int c, unsigned char *removecell,int64_t *removecount)
{

        double sc;
        double shift[3];

        if (cellsinfo->localcount.n + 1 > settings.maxlocalcells)
                terminate(systeminfo,"cannot allocate systeminfo->coords", __FILE__, __LINE__);

        sc = sqrt(cellsinfo->velocity[c].x * cellsinfo->velocity[c].x + cellsinfo->velocity[c].y * cellsinfo->velocity[c].y +
                  cellsinfo->velocity[c].z * cellsinfo->velocity[c].z);

        /* daughter cells are shifted away from the center of parent cell */
        /* direction of shift chosen randomly */
        int accept = 0;
        while (accept == 0) {
                shift[0] = ((double)rand()/RAND_MAX)*2.0 - 1.0;
                shift[1] = ((double)rand()/RAND_MAX)*2.0 - 1.0;
                if (settings.dimension == 3)
                        shift[2] = ((double)rand()/RAND_MAX)*2.0 - 1.0;
                else
                        shift[2] = 0.0;
                sc = sqrt(pow(shift[0], 2) + pow(shift[1], 2) + pow(shift[2], 2));
                if (sc == 0)
                        continue;
                sc = cellsinfo->cells[c].size / (2 * sc);
                shift[0] = sc * shift[0];
                shift[1] = sc * shift[1];
                shift[2] = sc * shift[2];
                accept = 1;
        }

        /* 1st daughter cell position, size, type and age */
        cellsinfo->cells[cellsinfo->localcount.n].x = cellsinfo->cells[c].x + shift[0];
        cellsinfo->cells[cellsinfo->localcount.n].y = cellsinfo->cells[c].y + shift[1];
        cellsinfo->cells[cellsinfo->localcount.n].z = cellsinfo->cells[c].z + shift[2];
        cellsinfo->cells[cellsinfo->localcount.n].size = pow(2.0, -(1.0 / 3.0)) * cellsinfo->cells[c].size;;
        cellsinfo->cells[cellsinfo->localcount.n].ctype = cellsinfo->cells[c].ctype;
        cellsinfo->cells[cellsinfo->localcount.n].age = cellsinfo->cells[c].age + 1;

        /* 2nd daughter cell position, size, type and age */
        cellsinfo->cells[c].x -= shift[0];
        cellsinfo->cells[c].y -= shift[1];
        cellsinfo->cells[c].z -= shift[2];
        cellsinfo->cells[c].size = cellsinfo->cells[cellsinfo->localcount.n].size;;
        cellsinfo->cells[c].age += 1;

        /* 2nd daughter cell cycle phases lenghts */
        cellsinfo->cells[c].g1 = celltype[cellsinfo->cells[c].ctype].g1 * (1 + (((float)rand()/RAND_MAX)*2.0 - 1.0) * celltype[cellsinfo->cells[c].ctype].v);
        cellsinfo->cells[c].g2 = celltype[cellsinfo->cells[c].ctype].g2 * (1 + (((float)rand()/RAND_MAX)*2.0 - 1.0) * celltype[cellsinfo->cells[c].ctype].v);
        cellsinfo->cells[c].s = celltype[cellsinfo->cells[c].ctype].s * (1 + (((float)rand()/RAND_MAX)*2.0 - 1.0) * celltype[cellsinfo->cells[c].ctype].v);
        cellsinfo->cells[c].m = celltype[cellsinfo->cells[c].ctype].m * (1 + (((float)rand()/RAND_MAX)*2.0 - 1.0) * celltype[cellsinfo->cells[c].ctype].v);

        /* 1st daughter cell global ID */
        cellsinfo->cells[cellsinfo->localcount.n].gid =
                (unsigned long long int) systeminfo.rank *(unsigned long long int)
                settings.maxlocalcells + (unsigned long long int) (cellsinfo->localcount.n);

        /* 1st daughter cell parameters */
        cellsinfo->cells[cellsinfo->localcount.n].v = 0.0;
        cellsinfo->cells[cellsinfo->localcount.n].density = cellsinfo->cells[c].density;
        cellsinfo->cells[cellsinfo->localcount.n].ctype = cellsinfo->cells[c].ctype;
        cellsinfo->cells[cellsinfo->localcount.n].young = 2100.0 + (((float)rand()/RAND_MAX)*2.0 - 1.0) * 100.0;
        cellsinfo->cells[cellsinfo->localcount.n].phase = 1;
        cellsinfo->cells[cellsinfo->localcount.n].death = 0;
        cellsinfo->cells[cellsinfo->localcount.n].phasetime = 0.0;
        cellsinfo->cells[cellsinfo->localcount.n].mindist = DBL_MAX;

        /* 1st daughter cell cycle phases lenghts */
        cellsinfo->cells[cellsinfo->localcount.n].g1 = celltype[cellsinfo->cells[cellsinfo->localcount.n].ctype].g1 * (1 + (((float)rand()/RAND_MAX)*2.0 - 1.0) * cellsinfo->cells[cellsinfo->localcount.n].v);
        cellsinfo->cells[cellsinfo->localcount.n].g2 = celltype[cellsinfo->cells[cellsinfo->localcount.n].ctype].g2 * (1 + (((float)rand()/RAND_MAX)*2.0 - 1.0) * cellsinfo->cells[cellsinfo->localcount.n].v);
        cellsinfo->cells[cellsinfo->localcount.n].s = celltype[cellsinfo->cells[cellsinfo->localcount.n].ctype].s * (1 + (((float)rand()/RAND_MAX)*2.0 - 1.0) * cellsinfo->cells[cellsinfo->localcount.n].v);
        cellsinfo->cells[cellsinfo->localcount.n].m = celltype[cellsinfo->cells[cellsinfo->localcount.n].ctype].m * (1 + (((float)rand()/RAND_MAX)*2.0 - 1.0) * cellsinfo->cells[cellsinfo->localcount.n].v);

        /* update local cell counters */

        cellsinfo->localtypecount[cellsinfo->cells[cellsinfo->localcount.n].ctype].n += 1;
        cellsinfo->localtypecount[cellsinfo->cells[cellsinfo->localcount.n].ctype].g1phase += 1;
        cellsinfo->localcount.n += 1;
        cellsinfo->localcount.g1phase += 1;

        return;

}

/*!
 * This function finds locates cell closest to the center of mass of the systeminfo
 * and marks this cell as a cancer cell.
 */
void markMiddleCancerCell(systeminfo_t systeminfo,cellsinfo_t *cellsinfo)
{
        int c;
        int middle = 0;
        double dist;
        struct {
                double val;
                int rank;
        } lmdist, gmdist;
        double center[3];
        double gcenter[3];

        /* each process computes its local center of mass */
        center[0] = 0.0;
        center[1] = 0.0;
        center[2] = 0.0;

        for (c = 0; c < cellsinfo->localcount.n; c++) {
                center[0] += cellsinfo->cells[c].x / cellsinfo->globalcount.n;
                center[1] += cellsinfo->cells[c].y / cellsinfo->globalcount.n;
                center[2] += cellsinfo->cells[c].z / cellsinfo->globalcount.n;
        }

        /* MPI Reduce operation computes global center of mass */
        MPI_Allreduce(center, gcenter, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        /* intialization */
        lmdist.rank = systeminfo.rank;
        lmdist.val = INT_MAX;

        /* each process finds local cell closest to the global center of mass */
        for (c = 0; c < cellsinfo->localcount.n; c++) {
                dist =  sqrt((cellsinfo->cells[c].x - gcenter[0]) * (cellsinfo->cells[c].x - gcenter[0]) +
                             (cellsinfo->cells[c].y - gcenter[1]) * (cellsinfo->cells[c].y - gcenter[1]) +
                             (cellsinfo->cells[c].z - gcenter[2]) * (cellsinfo->cells[c].z - gcenter[2]));
                if (dist < lmdist.val) {
                        lmdist.val = dist;
                        middle = c;
                }
        }

        /* MPI_Allreduce locates the cell closest to the global center of mass */
        MPI_Allreduce(&lmdist, &gmdist, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
        /* mark the found cell as cancer one */
        if (systeminfo.rank == gmdist.rank) {
                printf("Most middle cells id %d located on rank %d\n",middle,systeminfo.rank);
                /* you can implement action for middle cell here */
        }

        return;
}

/*!
 * This function dealocates all tables allocated during initialization of cell data
 */
/*void cellsdestroy()
   {
        int f;
        free(tlnc);
 #ifdef __MIC__
        _mm_free(cells);
 #else
        free(cells);
 #endif
        for (f = 0; f < NFIELDS; f++)
                free(cellFields[f]);
        free(cellFields);
 #ifdef __MIC__
        _mm_free(velocity);
 #else
        free(velocity);
 #endif
   } */

/*!
 * This function removes a dead cell from the simulation.
 */
void cellsdeath(int oldlnc,cellsinfo_t *cellsinfo,unsigned char *removecell,int64_t removecount)
{
        int c, pos;

        if(cellsinfo->globalcount.n==1) return;

        pos = 0;
        for (c = 0; c < cellsinfo->localcount.n; c++) {
                /* shift cells after dead cell removal */
                if (c >= oldlnc) {
                        cellsinfo->cells[pos] = cellsinfo->cells[c];
                        pos++;
                        continue;
                }
                if (c != pos && removecell[c] == 0)
                        cellsinfo->cells[pos] = cellsinfo->cells[c];
                if (removecell[c] == 0)
                        pos++;
                /* update cell counters */
                if (removecell[c] == 1) {
                        switch (cellsinfo->cells[c].phase) {
                        case 0:
                                cellsinfo->localcount.g0phase -= 1;
                                cellsinfo->localtypecount[cellsinfo->cells[c].ctype].g0phase -= 1;
                                break;
                        case 1:
                                cellsinfo->localcount.g1phase -= 1;
                                cellsinfo->localtypecount[cellsinfo->cells[c].ctype].g1phase -= 1;
                                break;
                        case 2:
                                cellsinfo->localcount.sphase -= 1;
                                cellsinfo->localtypecount[cellsinfo->cells[c].ctype].sphase -= 1;
                                break;
                        case 3:
                                cellsinfo->localcount.g2phase -= 1;
                                cellsinfo->localtypecount[cellsinfo->cells[c].ctype].g2phase -= 1;
                                break;
                        case 4:
                                cellsinfo->localcount.mphase -= 1;
                                cellsinfo->localtypecount[cellsinfo->cells[c].ctype].mphase -= 1;
                                break;
                        }
                        cellsinfo->localtypecount[cellsinfo->cells[c].ctype].n -= 1;
                }
        }
        cellsinfo->localcount.n -= removecount;
}


/*!
 * This function updates cells' positions.
 */
void updatepositions(settings_t settings,celltype_t *celltype,cellsinfo_t *cellsinfo,unsigned char *removecell, int64_t *removecount)
{
        int c;
        double alpha=0.00001;

        /* move cells */
        for (c = 0; c < cellsinfo->localcount.n; c++) {
                float randommove=celltype[cellsinfo->cells[c].ctype].randommove;
                cellsinfo->cells[c].x += cellsinfo->forces[c].x;
                cellsinfo->cells[c].y += cellsinfo->forces[c].y;
                if(settings.dimension==3) cellsinfo->cells[c].z += cellsinfo->forces[c].z;
                /* random movement */
                cellsinfo->cells[c].x += randommove*(((float)rand()/RAND_MAX)*2.0 - 1.0);
                cellsinfo->cells[c].y += randommove*(((float)rand()/RAND_MAX)*2.0 - 1.0);
                if(settings.dimension==3) cellsinfo->cells[c].z += randommove*(((float)rand()/RAND_MAX)*2.0 - 1.0);

                /* mark cells outside the box for removal */
                if( outsidethebox(cellsinfo,c) && removecell[c]==0) {
                        removecell[c]=1;
                        (*removecount)=(*removecount)+1;
                }
        }
        return;
}

/*!
 * This function updates cells' cycle phases.
 */
int updatecellcycles(systeminfo_t systeminfo,settings_t settings,celltype_t *celltype, cellsinfo_t *cellsinfo,cellenvdata_t **cellenvdata,unsigned char *removecell,int64_t *removecount)
{

        int c,f;
        double eps, csize, rd;
        int lnc;
        int lowenvlevel1,lowenvlevel2;

        lnc = cellsinfo->localcount.n;

        for (c = 0; c < lnc; c++) {

                eps = celltype[cellsinfo->cells[c].ctype].criticaldensity;
                csize = celltype[cellsinfo->cells[c].ctype].size;
                rd = celltype[cellsinfo->cells[c].ctype].rd;

                if (removecell[c])
                        continue;

                lowenvlevel1=0;
                lowenvlevel2=0;
                for(f=0; f<settings.numberoffields; f++) {
                        if(cellenvdata[f][c].value<celltype[cellsinfo->cells[c].ctype].criticallevel1[f]) { lowenvlevel1=1; break; }
                        if(cellenvdata[f][c].value<celltype[cellsinfo->cells[c].ctype].criticallevel2[f]) { lowenvlevel2=1; break; }
                }

                if (settings.simulationstart) {

                        if (cellsinfo->cells[c].phase != 0 && cellsinfo->cells[c].density <= eps)
                                cellsinfo->cells[c].phasetime += settings.secondsperstep / 3600.0;

                        switch (cellsinfo->cells[c].phase) {

                        case 0: /* G0 phase */

                                if (lowenvlevel2)                    /* low environemnt level */
                                {
                                        cellsinfo->cells[c].phase = 5;
                                        cellsinfo->cells[c].phasetime = 0;
                                        cellsinfo->localcount.g0phase -=1;
                                        cellsinfo->localcount.necroticphase +=1;
                                        break;
                                }
                                /* transition to G1 phase */
                                if(cellsinfo->cells[c].density <= eps /* enough space for growth */
                                   || cellsinfo->globalcount.n==1     /* single cell in the simulation */
                                   || (!lowenvlevel1) )                  /* sufficient level of nutrients */
                                {
                                        cellsinfo->cells[c].phase = 1;
                                        cellsinfo->localcount.g0phase -=1;
                                        cellsinfo->localcount.g1phase +=1;
                                        break;
                                }
                                break;

                        case 1: /* G1 phase */
                                /* transition to G0 or Necrotic phase */
                                if(cellsinfo->cells[c].density > eps  /* not enough space for growth */
                                   || lowenvlevel1 )                  /* environment level below level1 */
                                {
                                        if (lowenvlevel2 )            /* environment level below level2  -> necrotic phase */
                                        {
                                                cellsinfo->cells[c].phase = 5;
                                                cellsinfo->cells[c].phasetime = 0;
                                                cellsinfo->localcount.g1phase -=1;
                                                cellsinfo->localcount.necroticphase +=1;
                                        } else {                      /* environment level above level2 -> G0 phase */
                                                cellsinfo->cells[c].phase = 0;
                                                cellsinfo->localcount.g1phase -=1;
                                                cellsinfo->localcount.g0phase +=1;
                                        }
                                        break;
                                }

                                /* cells grow in phase G1 */
                                if (cellsinfo->cells[c].size < csize) {
                                        cellsinfo->cells[c].size += (csize - pow(2.0,-(1.0 / 3.0)) * csize) * (settings.secondsperstep) / (3600.0 * cellsinfo->cells[c].g1);
                                }
                                if (cellsinfo->cells[c].size > csize)
                                        cellsinfo->cells[c].size = csize;

                                /* change to S phase or die */
                                if (cellsinfo->cells[c].phasetime >= cellsinfo->cells[c].g1) {
                                        int death;
                                        cellsinfo->cells[c].phase = 2;
                                        cellsinfo->cells[c].phasetime = 0;
                                        cellsinfo->localcount.g1phase -=1;
                                        cellsinfo->localcount.sphase +=1;
                                        death = ((float)rand()/RAND_MAX < rd ? 1 : 0);
                                        if (death && removecell[c]==0) {
                                                removecell[c] = 1;
                                                (*removecount)=(*removecount)+1;
                                        }
                                }
                                break;
                        case 2: /* S phase */
                                if (lowenvlevel2) {
                                        cellsinfo->cells[c].phase = 5;
                                        cellsinfo->cells[c].phasetime = 0;
                                        cellsinfo->localcount.sphase -=1;
                                        cellsinfo->localcount.necroticphase +=1;
                                        break;
                                }
                                if (cellsinfo->cells[c].phasetime >= cellsinfo->cells[c].s) {
                                        cellsinfo->cells[c].phase = 3;
                                        cellsinfo->cells[c].phasetime = 0;
                                        cellsinfo->localcount.sphase -=1;
                                        cellsinfo->localcount.g2phase +=1;
                                        break;
                                }
                                break;
                        case 3: /* G2 phase */
                                if (lowenvlevel2) {
                                        cellsinfo->cells[c].phase = 5;
                                        cellsinfo->cells[c].phasetime = 0;
                                        cellsinfo->localcount.g2phase -=1;
                                        cellsinfo->localcount.necroticphase +=1;
                                        break;
                                }
                                if (cellsinfo->cells[c].phasetime >= cellsinfo->cells[c].g2) {
                                        int death;
                                        cellsinfo->cells[c].phase = 4;
                                        cellsinfo->cells[c].phasetime = 0;
                                        cellsinfo->localcount.g2phase -=1;
                                        cellsinfo->localcount.mphase +=1;
                                        death = ((float)rand()/RAND_MAX < rd ? 1 : 0);
                                        if (death && removecell[c]==0) {
                                                removecell[c] = 1;
                                                (*removecount)=(*removecount)+1;
                                        }
                                        break;
                                }
                                break;
                        case 4: /* M phase */
                                if (lowenvlevel2) {
                                        cellsinfo->cells[c].phase = 5;
                                        cellsinfo->cells[c].phasetime = 0;
                                        cellsinfo->localcount.mphase -=1;
                                        cellsinfo->localcount.necroticphase +=1;

                                } else if (cellsinfo->cells[c].phasetime >= cellsinfo->cells[c].m) {
                                        mitosis(systeminfo,settings,celltype,cellsinfo,c,removecell,removecount);
                                        cellsinfo->cells[c].phase = 1;
                                        cellsinfo->cells[c].phasetime = 0;
                                        cellsinfo->localcount.mphase -=1;
                                        cellsinfo->localcount.g1phase +=1;
                                }
                                break;
                        } // switch
                } // if
        } // for loop

        return 0;
}

/*!
 * This function drives the whole cell cycle update.
 */
void cellsupdate(systeminfo_t systeminfo, settings_t settings,celltype_t *celltype,cellenvdata_t **cellenvdata,cellsinfo_t *cellsinfo)
{
        int oldlnc;
        unsigned char *removecell;
        int64_t removecount;
        /* number of local cells might change during the update */
        oldlnc = cellsinfo->localcount.n;
        removecell = (unsigned char *) calloc(oldlnc, sizeof(unsigned char));
        removecount = 0;

        updatepositions(settings,celltype,cellsinfo,removecell,&removecount);
        updatecellcycles(systeminfo,settings,celltype,cellsinfo,cellenvdata,removecell,&removecount);
        cellsdeath(oldlnc,cellsinfo,removecell,removecount);
        updateglobalcounts(cellsinfo);

        free(removecell);
}
