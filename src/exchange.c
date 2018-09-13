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

#include <stdlib.h>

#include "global.h"
#include "utils.h"


/*!
 * This function is a comparison function used
 * for sorting export list table.
 */
int explistcompare(const void *a, const void *b)
{
        return ((explist_t *) a)->proc - ((explist_t *) b)->proc;
}

/*!
 * This function uses Zoltan's library function Zoltan_LB_Box_Assign
 * to find possible intersections of cells' neighbourhoods
 * and other processes' geometries.
 */
void createexportlist(systeminfo_t systeminfo,settings_t settings,cellsinfo_t cellsinfo,grid_t grid,struct Zoltan_Struct *ztn,celltype_t* celltype,cellcommdata_t *cellcommdata,fieldcommdata_t *fieldcommdata)
{

        int i, p;
        int procs[systeminfo.size];
        int numprocs;

        if (cellsinfo.globalcount.n < systeminfo.size*MIN_CELLS_PER_PROC || systeminfo.size == 1)
                return;

        /* 1. create export list for cellular dynamics computations */
        /* reset counters */
        cellcommdata->numexp = 0;
        cellcommdata->numimp = 0;
        cellcommdata->explistmaxsize=settings.maxlocalcells;
        /* allocate tables */
        if(!(cellcommdata->explist = (explist_t*) malloc(sizeof(explist_t) * cellcommdata->explistmaxsize)))
                terminate(systeminfo,"cannot allocate cellcommdata->explist", __FILE__, __LINE__);
        if(!(cellcommdata->recvcount = (int *) calloc(systeminfo.size, sizeof(int))))
                terminate(systeminfo,"cannot allocate cellcommdata->recvcount", __FILE__, __LINE__);
        if(!(cellcommdata->sendcount = (int *) calloc(systeminfo.size, sizeof(int))))
                terminate(systeminfo,"cannot allocate cellcommdata->sendcount", __FILE__, __LINE__);
        if(!(cellcommdata->sendoffset = (int64_t *) calloc(systeminfo.size, sizeof(int64_t))))
                terminate(systeminfo,"cannot allocate cellcommdata->sendoffset", __FILE__, __LINE__);
        if(!(cellcommdata->recvoffset = (int64_t *) calloc(systeminfo.size, sizeof(int64_t))))
                terminate(systeminfo,"cannot allocate cellcommdata->recvoffset", __FILE__, __LINE__);
        /* loop over local cells */
        /*#pragma omp parallel for private(procs) */
        for (p = 0; p < cellsinfo.localcount.n; p++) {
                double xmin, xmax, ymin, ymax, zmin, zmax;
                double r;
                if (cellsinfo.globalcount.n < systeminfo.size*MIN_CELLS_PER_PROC)
                        continue;
                r = celltype[cellsinfo.cells[p].ctype].h * 1.5;
                /* compute neighbourhood box */
                xmin = cellsinfo.cells[p].x - r;
                xmax = cellsinfo.cells[p].x + r;
                ymin = cellsinfo.cells[p].y - r;
                ymax = cellsinfo.cells[p].y + r;
                if (cellsinfo.dimension == 3) {
                        zmin = cellsinfo.cells[p].z - r;
                        zmax = cellsinfo.cells[p].z + r;
                } else {
                        zmin = 0.0;
                        zmax = 0.0;
                }
                /* look for possible neighbours */
                Zoltan_LB_Box_Assign(ztn, xmin, ymin, zmin, xmax, ymax, zmax, procs, &numprocs);
                /* loop over receivers */
                for (i = 0; i < numprocs; i++) {
                        if (procs[i] == systeminfo.rank || cellsinfo.cellsperproc[procs[i]] == 0)
                                continue;
                        cellcommdata->explist[cellcommdata->numexp].cell = p;
                        cellcommdata->explist[cellcommdata->numexp].proc = procs[i];
                        cellcommdata->sendcount[procs[i]]++;
                        cellcommdata->numexp++;
                        /* too many refugees  - reallocate */
                        if (cellcommdata->numexp >= cellcommdata->explistmaxsize) {
                                cellcommdata->explistmaxsize+=64;
                                if(!(cellcommdata->explist=(explist_t*)realloc(cellcommdata->explist,sizeof(explist_t)*cellcommdata->explistmaxsize))) {
                                        terminate(systeminfo,"cannot reallocate cellcommdata->explist", __FILE__, __LINE__);
                                }
                        }
                }
        }
        /* sort export list with respect to process number */
        qsort(cellcommdata->explist, cellcommdata->numexp, sizeof(explist_t), explistcompare);
        /* distribute the information on transfer sizes between each process */
        MPI_Alltoall(cellcommdata->sendcount, 1, MPI_INT, cellcommdata->recvcount, 1, MPI_INT, MPI_COMM_WORLD);
        /* compute send offsets */
        cellcommdata->sendoffset[0] = 0;
        for (i = 1; i < systeminfo.size; i++)
                cellcommdata->sendoffset[i] = cellcommdata->sendoffset[i - 1] + cellcommdata->sendcount[i - 1];
        /* compute receive offsets */
        cellcommdata->recvoffset[0] = 0;
        for (i = 1; i < systeminfo.size; i++)
                cellcommdata->recvoffset[i] = cellcommdata->recvoffset[i - 1] + cellcommdata->recvcount[i - 1];
        /* count cells to be imported */
        for (i = 0; i < systeminfo.size; i++)
                cellcommdata->numimp += cellcommdata->recvcount[i];

        /* 2. create export list for environment computations */
        /* reset counters */
        fieldcommdata->numexp = 0;
        fieldcommdata->numimp = 0;
        fieldcommdata->explistmaxsize=settings.maxlocalcells;
        /* allocate tables */
        if(!(fieldcommdata->explist = (explist_t*) malloc(sizeof(explist_t) * fieldcommdata->explistmaxsize)))
                terminate(systeminfo,"cannot allocate cellcommdata->explist", __FILE__, __LINE__);
        if(!(fieldcommdata->recvcount = (int *) calloc(systeminfo.size, sizeof(int))))
                terminate(systeminfo,"cannot allocate cellcommdata->recvcount", __FILE__, __LINE__);
        if(!(fieldcommdata->sendcount = (int *) calloc(systeminfo.size, sizeof(int))))
                terminate(systeminfo,"cannot allocate cellcommdata->sendcount", __FILE__, __LINE__);
        if(!(fieldcommdata->sendoffset = (int64_t *) calloc(systeminfo.size, sizeof(int64_t))))
                terminate(systeminfo,"cannot allocate cellcommdata->sendoffset", __FILE__, __LINE__);
        if(!(fieldcommdata->recvoffset = (int64_t *) calloc(systeminfo.size, sizeof(int64_t))))
                terminate(systeminfo,"cannot allocate cellcommdata->recvoffset", __FILE__, __LINE__);
        /* loop over local cells */
        /*#pragma omp parallel for private(procs) */
        for (p = 0; p < cellsinfo.localcount.n; p++) {
                double x,y,z;
                x=cellsinfo.cells[p].x;
                y=cellsinfo.cells[p].y;
                z=cellsinfo.cells[p].z;
                for(i=0; i<systeminfo.size; i++) {
                        if(i==systeminfo.rank) continue;
                        if(x>=grid.lowleftnear[i].x && x<grid.uprightfar[i].x &&
                           y>=grid.lowleftnear[i].y && y<grid.uprightfar[i].y &&
                           z>=grid.lowleftnear[i].z && z<grid.uprightfar[i].z ) {
                                fieldcommdata->explist[fieldcommdata->numexp].cell = p;
                                fieldcommdata->explist[fieldcommdata->numexp].proc = i;
                                fieldcommdata->sendcount[i]++;
                                fieldcommdata->numexp++;
                                break;
                        }
                }
        }
        /* sort export list with respect to process number */
        qsort(fieldcommdata->explist, fieldcommdata->numexp, sizeof(explist_t), explistcompare);
        /* distribute the information on transfer sizes between each process */
        MPI_Alltoall(fieldcommdata->sendcount, 1, MPI_INT, fieldcommdata->recvcount, 1, MPI_INT, MPI_COMM_WORLD);
        /* compute send offsets */
        fieldcommdata->sendoffset[0] = 0;
        for (i = 1; i < systeminfo.size; i++)
                fieldcommdata->sendoffset[i] = fieldcommdata->sendoffset[i - 1] + fieldcommdata->sendcount[i - 1];
        /* compute receive offsets */
        fieldcommdata->recvoffset[0] = 0;
        for (i = 1; i < systeminfo.size; i++)
                fieldcommdata->recvoffset[i] = fieldcommdata->recvoffset[i - 1] + fieldcommdata->recvcount[i - 1];
        /* count cells to be imported */
        for (i = 0; i < systeminfo.size; i++)
                fieldcommdata->numimp += fieldcommdata->recvcount[i];

        return;
}

/*!
 * This function deallocates all communication buffers and
 * auxiliary tables.
 */
void exchangecleanup(systeminfo_t systeminfo,cellsinfo_t cellsinfo,cellcommdata_t *cellcommdata,fieldcommdata_t *fieldcommdata)
{
        if (cellsinfo.globalcount.n < systeminfo.size*MIN_CELLS_PER_PROC || systeminfo.size == 1)
                return;
        /* deallocate cell dynamics exchange arrays */
        free(cellcommdata->recvcellindata);
        free(cellcommdata->recvcelloutdata);
        free(cellcommdata->explist);
        free(cellcommdata->recvcount);
        free(cellcommdata->sendcount);
        free(cellcommdata->sendoffset);
        free(cellcommdata->recvoffset);
        /* deallocate environment exchange arrays */
        //free(fieldcommdata->recvfieldindata);
        //free(fieldcommdata->recvfieldoutdata);
        free(fieldcommdata->explist);
        free(fieldcommdata->recvcount);
        free(fieldcommdata->sendcount);
        free(fieldcommdata->sendoffset);
        free(fieldcommdata->recvoffset);
        return;
}

/*!
 * This function initiate sending and receiving cells' data between processes.
 */
void cellssendrecv(systeminfo_t systeminfo, cellsinfo_t cellsinfo, cellcommdata_t *cellcommdata)
{
        int i;

        if (cellsinfo.globalcount.n < systeminfo.size*MIN_CELLS_PER_PROC || systeminfo.size == 1)
                return;

        cellindata_t *sendcellindata;
        cellindata_t *recvcellindata;
        celloutdata_t *sendcelloutdata;
        celloutdata_t *recvcelloutdata;
        /* allocate communication buffers */
        cellcommdata->sendcellindata = (cellindata_t *) malloc(cellcommdata->numexp * sizeof(cellindata_t));
        cellcommdata->recvcellindata = (cellindata_t *) malloc(cellcommdata->numimp * sizeof(cellindata_t));
        cellcommdata->reqsend = (MPI_Request *) malloc(sizeof(MPI_Request) * systeminfo.size);
        cellcommdata->reqrecv = (MPI_Request *) malloc(sizeof(MPI_Request) * systeminfo.size);

        /* create reduced particle data buffer for exporting */
        for (i = 0; i < cellcommdata->numexp; i++) {
                int64_t cellidx=cellcommdata->explist[i].cell;
                cellcommdata->sendcellindata[i].x = cellsinfo.cells[cellidx].x;
                cellcommdata->sendcellindata[i].y = cellsinfo.cells[cellidx].y;
                cellcommdata->sendcellindata[i].z = cellsinfo.cells[cellidx].z;
                cellcommdata->sendcellindata[i].size = cellsinfo.cells[cellidx].size;
                cellcommdata->sendcellindata[i].young = (double) cellsinfo.cells[cellidx].young;
                cellcommdata->sendcellindata[i].ctype = cellsinfo.cells[cellidx].ctype;
        }

        /* send cells - asynchronous MPI call */
        for (i = 0; i < systeminfo.size; i++) {
                if (cellcommdata->sendcount[i] == 0 || cellsinfo.cellsperproc[i] == 0 || cellsinfo.localcount.n == 0)
                        continue;
                MPI_Isend(&(cellcommdata->sendcellindata[cellcommdata->sendoffset[i]]),
                          cellcommdata->sendcount[i] * sizeof(cellindata_t), MPI_BYTE, i, systeminfo.rank,
                          MPI_COMM_WORLD, &(cellcommdata->reqsend[i]));
        }

        /* receive cells - asynchronous MPI call */
        for (i = 0; i < systeminfo.size; i++) {
                if (cellcommdata->recvcount[i] == 0 || cellsinfo.cellsperproc[i] == 0 || cellsinfo.localcount.n == 0)
                        continue;
                MPI_Irecv(&(cellcommdata->recvcellindata[cellcommdata->recvoffset[i]]),
                          cellcommdata->recvcount[i] * sizeof(cellindata_t), MPI_BYTE, i, i,
                          MPI_COMM_WORLD, &(cellcommdata->reqrecv[i]));
        }

        return;
}

/*!
 * This function waits for cells' data exchange completion.
 */
void cellswait(systeminfo_t systeminfo, cellsinfo_t cellsinfo, cellcommdata_t *cellcommdata)
{
        int i;
        MPI_Status status;

        if (cellsinfo.globalcount.n < systeminfo.size*MIN_CELLS_PER_PROC || systeminfo.size == 1)
                return;

        /* wait for send completion */
        for (i = 0; i < systeminfo.size; i++) {
                if (cellcommdata->sendcount[i] == 0 || cellsinfo.cellsperproc[i] == 0 || cellsinfo.localcount.n == 0)
                        continue;
                if (MPI_Wait(&(cellcommdata->reqsend[i]), &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication error", __FILE__, __LINE__);
        }

        /* wait for receive completion */
        for (i = 0; i < systeminfo.size; i++) {
                if (cellcommdata->recvcount[i] == 0 || cellsinfo.cellsperproc[i] == 0 || cellsinfo.localcount.n == 0)
                        continue;
                if (MPI_Wait(&(cellcommdata->reqrecv[i]), &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication error", __FILE__, __LINE__);
        }

        /* some of the buffers can be deallocated here */
        free(cellcommdata->sendcellindata);
        free(cellcommdata->reqsend);
        free(cellcommdata->reqrecv);
}

/*!
 * This function initiate sending and receiving density
 * and potential values between processes.
 */
void datasendrecv(systeminfo_t systeminfo, cellsinfo_t cellsinfo, cellcommdata_t *cellcommdata)
{
        int i;

        if (cellsinfo.globalcount.n < systeminfo.size*MIN_CELLS_PER_PROC || systeminfo.size == 1)
                return;

        /* allocate communication buffers */
        cellcommdata->sendcelloutdata =
                (celloutdata_t *) malloc(cellcommdata->numexp * sizeof(celloutdata_t));
        cellcommdata->recvcelloutdata =
                (celloutdata_t *) malloc(cellcommdata->numimp * sizeof(celloutdata_t));
        cellcommdata->reqsend = (MPI_Request *) malloc(sizeof(MPI_Request) * systeminfo.size);
        cellcommdata->reqrecv = (MPI_Request *) malloc(sizeof(MPI_Request) * systeminfo.size);

        /* create density and potential buffer for exporting */
        for (i = 0; i < cellcommdata->numexp; i++) {
                int64_t cellidx=cellcommdata->explist[i].cell;
                cellcommdata->sendcelloutdata[i].v = cellsinfo.cells[cellidx].v;
                cellcommdata->sendcelloutdata[i].density = cellsinfo.cells[cellidx].density;
        }

        /* send data - asynchronous MPI call */
        for (i = 0; i < systeminfo.size; i++) {
                if (cellcommdata->sendcount[i] == 0 || cellsinfo.cellsperproc[i] == 0 || cellsinfo.localcount.n == 0)
                        continue;
                MPI_Isend(&(cellcommdata->sendcelloutdata[cellcommdata->sendoffset[i]]),
                          cellcommdata->sendcount[i] * sizeof(celloutdata_t), MPI_BYTE, i,
                          systeminfo.rank, MPI_COMM_WORLD, &(cellcommdata->reqsend[i]));
        }

        /* receive data - asynchronous MPI call */
        for (i = 0; i < systeminfo.size; i++) {
                if (cellcommdata->recvcount[i] == 0 || cellsinfo.cellsperproc[i] == 0 || cellsinfo.localcount.n == 0)
                        continue;
                MPI_Irecv(&(cellcommdata->recvcelloutdata[cellcommdata->recvoffset[i]]),
                          cellcommdata->recvcount[i] * sizeof(celloutdata_t), MPI_BYTE, i, i,
                          MPI_COMM_WORLD, &(cellcommdata->reqrecv[i]));
        }

        return;
}

/*!
 * This function waits for density and potential data exchange completion.
 */
void datawait(systeminfo_t systeminfo, cellsinfo_t cellsinfo, cellcommdata_t *cellcommdata)
{
        int i;
        MPI_Status status;

        if (cellsinfo.globalcount.n < systeminfo.size*MIN_CELLS_PER_PROC || systeminfo.size == 1)
                return;

        // Wait for send completion
        for (i = 0; i < systeminfo.size; i++) {
                if (cellcommdata->sendcount[i] == 0 || cellsinfo.cellsperproc[i] == 0 || cellsinfo.localcount.n == 0)
                        continue;
                if (MPI_Wait(&(cellcommdata->reqsend[i]), &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication error", __FILE__, __LINE__);
        }

        // Wait for receive completion
        for (i = 0; i < systeminfo.size; i++) {
                if (cellcommdata->recvcount[i] == 0 || cellsinfo.cellsperproc[i] == 0 || cellsinfo.localcount.n == 0)
                        continue;
                if (MPI_Wait(&(cellcommdata->reqrecv[i]), &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication error", __FILE__, __LINE__);
        }
        /* some of the buffers can be deallocated */
        free(cellcommdata->sendcelloutdata);
        free(cellcommdata->reqsend);
        free(cellcommdata->reqrecv);
        return;
}
