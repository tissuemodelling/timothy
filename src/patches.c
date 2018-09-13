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
#include <math.h>
#include <limits.h>

#include "global.h"
#include "patches.h"
#include "utils.h"

/*! \file interp.c
 *  \brief contains grid to cellular data interpolation functions
 */

#define patchelement(p,f,i,j,k) (patches->data[p][patches->size[p].x*patches->size[p].y*patches->size[p].z*f+patches->size[p].y*patches->size[p].z*i+patches->size[p].z*j+k])

/*!
 * For each local cell we check its position in the grid.
 * If cell is located in the grid partition of other process
 * then the information about this cell should be send to
 * remote process by MPI communication. This is done by preparing
 * special patches and sending the information on all bunch of cells
 * rather than sendind it one by one.
 * This function identifies patches and allocate memory buffers for patches.
 */
void patches_alloc(systeminfo_t systeminfo, settings_t settings, patches_t *patches, cellsinfo_t *cellsinfo, grid_t *grid)
{

        int c,p;
        int643dv_t cellIdx;

        patches->data = (double **) calloc(systeminfo.size, sizeof(double *));
        patches->intersect = (int *) calloc(systeminfo.size, sizeof(int));
        patches->lowercorner = (int643dv_t *) calloc(systeminfo.size, sizeof(int643dv_t));
        patches->uppercorner = (int643dv_t *) calloc(systeminfo.size, sizeof(int643dv_t));
        patches->lowercornerR = (int643dv_t *) calloc(systeminfo.size, sizeof(int643dv_t));
        patches->uppercornerR = (int643dv_t *) calloc(systeminfo.size, sizeof(int643dv_t));
        patches->size = (int643dv_t *) calloc(systeminfo.size, sizeof(int643dv_t));

        for (p = 0; p < systeminfo.size; p++) {
                patches->intersect[p] = 0;
                patches->size[p].x = 0;
                patches->size[p].y = 0;
                patches->size[p].z = 0;
                patches->lowercorner[p].x = INT_MAX;
                patches->lowercorner[p].y = INT_MAX;
                if (cellsinfo->dimension == 3)
                        patches->lowercorner[p].z = INT_MAX;
                else
                        patches->lowercorner[p].z = 0;
                patches->uppercorner[p].x = INT_MIN;
                patches->uppercorner[p].y = INT_MIN;
                if (cellsinfo->dimension == 3)
                        patches->uppercorner[p].z = INT_MIN;
                else
                        patches->uppercorner[p].z = 0;
        }

        //#pragma omp parallel for default(none) private(p,c,cellIdx) shared(cells,gridResolution,lnc,MPIsize,gridStartIdx,gridEndIdx,lowercorner,uppercorner,intersect,sdim,lowerGridCorner)
        for (p = 0; p < systeminfo.size; p++) {

                for (c = 0; c < cellsinfo->localcount.n; c++) {

                        int ax, ay, az;

                        cellIdx.x = ((cellsinfo->cells[c].x - grid->lowercorner.x) / grid->resolution);
                        cellIdx.y = ((cellsinfo->cells[c].y - grid->lowercorner.y) / grid->resolution);
                        cellIdx.z = ((cellsinfo->cells[c].z - grid->lowercorner.z) / grid->resolution);

                        for (ax = 0; ax < 2; ax++)
                                for (ay = 0; ay < 2; ay++)
                                        for (az = 0; az < 2; az++) {

                                                if (cellIdx.x + ax >= grid->loweridx[p].x
                                                    && cellIdx.y + ay >= grid->loweridx[p].y
                                                    && cellIdx.z + az >= grid->loweridx[p].z
                                                    && cellIdx.x + ax <= grid->upperidx[p].x
                                                    && cellIdx.y + ay <= grid->upperidx[p].y
                                                    && cellIdx.z + az <= grid->upperidx[p].z) {
                                                        patches->lowercorner[p].x =
                                                                (patches->lowercorner[p].x >
                                                                 cellIdx.x + ax ? cellIdx.x +
                                                                 ax : patches->lowercorner[p].x);
                                                        patches->lowercorner[p].y =
                                                                (patches->lowercorner[p].y >
                                                                 cellIdx.y + ay ? cellIdx.y +
                                                                 ay : patches->lowercorner[p].y);
                                                        if (cellsinfo->dimension == 3)
                                                                patches->lowercorner[p].z =
                                                                        (patches->lowercorner[p].z >
                                                                         cellIdx.z + az ? cellIdx.z +
                                                                         az : patches->lowercorner[p].z);
                                                        patches->uppercorner[p].x =
                                                                (patches->uppercorner[p].x <
                                                                 cellIdx.x + ax ? cellIdx.x +
                                                                 ax : patches->uppercorner[p].x);
                                                        patches->uppercorner[p].y =
                                                                (patches->uppercorner[p].y <
                                                                 cellIdx.y + ay ? cellIdx.y +
                                                                 ay : patches->uppercorner[p].y);
                                                        if (cellsinfo->dimension == 3)
                                                                patches->uppercorner[p].z =
                                                                        (patches->uppercorner[p].z <
                                                                         cellIdx.z + az ? cellIdx.z +
                                                                         az : patches->uppercorner[p].z);

                                                        patches->intersect[p] = 1;
                                                }
                                        }
                }
        }

        for (p = 0; p < systeminfo.size; p++)
                if (patches->intersect[p]) {
                        patches->size[p].x = patches->uppercorner[p].x - patches->lowercorner[p].x + 1;
                        patches->size[p].y = patches->uppercorner[p].y - patches->lowercorner[p].y + 1;
                        if (cellsinfo->dimension == 3)
                                patches->size[p].z = patches->uppercorner[p].z - patches->lowercorner[p].z + 1;
                        else
                                patches->size[p].z = 1;
                        patches->data[p] =
                                (double *) calloc(patches->size[p].x * patches->size[p].y *
                                                  patches->size[p].z * settings.numberoffields, sizeof(double));
                }

        return;
}

/*!
 * For each local cell its density value is interpolated accross
 * neighbouring grid vertices with the use of Cloud-In-Cell method.
 * Computed values are stored in patches instead in field buffers.
 * No additional memory allocations are made here.
 */
void patches_cells2envinfo(systeminfo_t systeminfo, settings_t settings, patches_t *patches, cellsinfo_t *cellsinfo, grid_t *grid)
{

        int c, i, p, j, k, f;
        int643dv_t idx;
        double3dv_t d, t;
        int643dv_t cellIdx;
        double3dv_t cicCoord;

        for (p = 0; p < systeminfo.size; p++)
                for(f = 0; f < settings.numberoffields; f++)
                        for (i = 0; i < patches->size[p].x; i++)
                                for (j = 0; j < patches->size[p].y; j++)
                                        for (k = 0; k < patches->size[p].z; k++)
                                                patchelement(p, f, i, j, k) = 0.0;

        for (c = 0; c < cellsinfo->localcount.n; c++) {

                cellIdx.x = ((cellsinfo->cells[c].x - grid->lowercorner.x) / grid->resolution);
                cellIdx.y = ((cellsinfo->cells[c].y - grid->lowercorner.y) / grid->resolution);
                cellIdx.z = ((cellsinfo->cells[c].z - grid->lowercorner.z) / grid->resolution);

                for (p = 0; p < systeminfo.size; p++) {
                        int ax, ay, az;
                        for (ax = 0; ax < 2; ax++)
                                for (ay = 0; ay < 2; ay++)
                                        for (az = 0; az < 2; az++) {
                                                if (cellIdx.x + ax >= grid->loweridx[p].x
                                                    && cellIdx.y + ay >= grid->loweridx[p].y
                                                    && cellIdx.z + az >= grid->loweridx[p].z
                                                    && cellIdx.x + ax <= grid->upperidx[p].x
                                                    && cellIdx.y + ay <= grid->upperidx[p].y
                                                    && cellIdx.z + az <= grid->upperidx[p].z) {

                                                        idx.x = (cellIdx.x + ax) - patches->lowercorner[p].x;
                                                        idx.y = (cellIdx.y + ay) - patches->lowercorner[p].y;
                                                        idx.z = (cellIdx.z + az) - patches->lowercorner[p].z;

                                                        cicCoord.x = grid->lowercorner.x + cellIdx.x * grid->resolution;
                                                        cicCoord.y = grid->lowercorner.y + cellIdx.y * grid->resolution;
                                                        cicCoord.z = grid->lowercorner.z + cellIdx.z * grid->resolution;

                                                        d.x = (cellsinfo->cells[c].x - cicCoord.x) / grid->resolution;
                                                        d.y = (cellsinfo->cells[c].y - cicCoord.y) / grid->resolution;
                                                        d.z = (cellsinfo->cells[c].z - cicCoord.z) / grid->resolution;

                                                        t.x = 1.0 - d.x;
                                                        t.y = 1.0 - d.y;
                                                        t.z = 1.0 - d.z;

                                                        /*if(cells[c].ctype==1) { endothelial cell - production
                                                                patch(p,1,idx.x,idx.y,idx.z) +=
                                                                        1.0 * (ax * d.x + (1 - ax) * t.x) * (ay * d.y +
                                                                                                             (1 -
                                                                                                              ay) * t.y) *
                                                                        (az * d.z + (1 - az) * t.z);
                                                           } else */
                                                        if (cellsinfo->cells[c].phase != 5) { /* if not in necrotic phase */
                                                                if (cellsinfo->cells[c].phase == 0) { /* if in G0 phase - lower consumption */
                                                                        patchelement(p,0, idx.x, idx.y, idx.z) +=
                                                                                0.75 * (ax * d.x + (1 - ax) * t.x) * (ay * d.y +
                                                                                                                      (1 -
                                                                                                                       ay) * t.y) *
                                                                                (az * d.z + (1 - az) * t.z);
                                                                } else { /* if not in G0 phase - normal consumption */
                                                                        patchelement(p,0, idx.x, idx.y, idx.z) +=
                                                                                1.0 * (ax * d.x + (1 - ax) * t.x) * (ay * d.y +
                                                                                                                     (1 -
                                                                                                                      ay) * t.y) *
                                                                                (az * d.z + (1 - az) * t.z);
                                                                }
                                                        }

                                                }
                                        }
                }
        }
        return;
}

/*!
 *  Patches are being sent to receiving processes with non-blocking
 * communication scheme (Isend, Irecv).
 * Receiving patches are allocated here.
 * MPI_Request tables are allocated here.
 */
void patches_commcells2envinit(systeminfo_t systeminfo, settings_t settings, patches_t *patches)
{
        int p;

        patches->receiver = (int *) calloc(systeminfo.size, sizeof(int));
        patches->sender = (int *) calloc(systeminfo.size, sizeof(int));

        for (p = 0; p < systeminfo.size; p++)
                patches->receiver[p] = patches->intersect[p];

        MPI_Alltoall(patches->receiver, 1, MPI_INT, patches->sender, 1, MPI_INT,
                     MPI_COMM_WORLD);

        MPI_Alltoall(patches->lowercorner, sizeof(int643dv_t), MPI_BYTE,
                     patches->lowercornerR, sizeof(int643dv_t), MPI_BYTE,
                     MPI_COMM_WORLD);
        MPI_Alltoall(patches->uppercorner, sizeof(int643dv_t), MPI_BYTE,
                     patches->uppercornerR, sizeof(int643dv_t), MPI_BYTE,
                     MPI_COMM_WORLD);

        patches->recvdata = (double **) calloc(systeminfo.size, sizeof(double *));
        patches->reqsend = (MPI_Request *) malloc(sizeof(MPI_Request) * systeminfo.size);
        patches->reqrecv = (MPI_Request *) malloc(sizeof(MPI_Request) * systeminfo.size);

        for (p = 0; p < systeminfo.size; p++) {
                if (patches->receiver[p]) {
                        MPI_Isend(&(patches->data[p][0]),
                                  patches->size[p].x * patches->size[p].y * patches->size[p].z * settings.numberoffields,
                                  MPI_DOUBLE, p, systeminfo.rank, MPI_COMM_WORLD, &(patches->reqsend[p]));
                }
                if (patches->sender[p]) {
                        int recvSize;
                        recvSize =
                                (patches->uppercornerR[p].x - patches->lowercornerR[p].x +
                                 1) * (patches->uppercornerR[p].y - patches->lowercornerR[p].y +
                                       1) * (patches->uppercornerR[p].z - patches->lowercornerR[p].z +
                                             1);
                        if (!(patches->recvdata[p] = (double *) calloc(recvSize*settings.numberoffields, sizeof(double))))
                                terminate(systeminfo,"cannot allocate recvdata", __FILE__, __LINE__);
                        MPI_Irecv(&(patches->recvdata[p][0]), recvSize * settings.numberoffields, MPI_DOUBLE, p, p,
                                  MPI_COMM_WORLD, &(patches->reqrecv[p]));
                }
        }
        return;
}

/*!
 * Wait for communication to finish.
 * MPI_Request tables are deallocated here.
 */
int patches_commcells2envwait(systeminfo_t systeminfo, patches_t *patches)
{
        int p;
        MPI_Status status;

        for (p = 0; p < systeminfo.size; p++) {
                if (!patches->receiver)
                        continue;
                if (MPI_Wait(&(patches->reqsend[p]), &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication error commcells2envwait 1", __FILE__, __LINE__);
        }

        for (p = 0; p < systeminfo.size; p++) {
                if (!patches->sender[p])
                        continue;
                if (MPI_Wait(&(patches->reqrecv[p]), &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication error commcells2envwait 2", __FILE__, __LINE__);
        }

        free(patches->reqsend);
        free(patches->reqrecv);

        return 0;
}

/*!
 * Update local tissueField part with information from remote processes
 * received in patches. Receiveing patches are deallocated here.
 */
int patches_applycells2env(systeminfo_t systeminfo, settings_t settings, patches_t *patches, grid_t *grid,environment_t **environment)
{

        int i,p,f;
        int index;

        for(f=0; f<settings.numberoffields; f++)
                memset((*environment)[f].production,0,grid->localsize.x*grid->localsize.y*grid->localsize.z*sizeof(double));

        for (p = 0; p < systeminfo.size; p++) {
                int i, j, k;
                if (!patches->sender[p])
                        continue;
                for (i = patches->lowercornerR[p].x; i <= patches->uppercornerR[p].x; i++)
                        for (j = patches->lowercornerR[p].y; j <= patches->uppercornerR[p].y; j++)
                                for (k = patches->lowercornerR[p].z; k <= patches->uppercornerR[p].z; k++) {
                                        int643dv_t c, g, size;
                                        size.x = patches->uppercornerR[p].x - patches->lowercornerR[p].x + 1;
                                        size.y = patches->uppercornerR[p].y - patches->lowercornerR[p].y + 1;
                                        size.z = patches->uppercornerR[p].z - patches->lowercornerR[p].z + 1;
                                        c.x = i - patches->lowercornerR[p].x;
                                        c.y = j - patches->lowercornerR[p].y;
                                        c.z = k - patches->lowercornerR[p].z;
                                        if (i >= grid->loweridx[systeminfo.rank].x && i <= grid->upperidx[systeminfo.rank].x
                                            && j >= grid->loweridx[systeminfo.rank].y && j <= grid->upperidx[systeminfo.rank].y
                                            && k >= grid->loweridx[systeminfo.rank].z && k <= grid->upperidx[systeminfo.rank].z) {
                                                g.x = i - grid->loweridx[systeminfo.rank].x;
                                                g.y = j - grid->loweridx[systeminfo.rank].y;
                                                g.z = k - grid->loweridx[systeminfo.rank].z;
                                                index=grid->localsize.z * grid->localsize.y * g.x + grid->localsize.z * g.y + g.z;
                                                for(f=0; f<settings.numberoffields; f++)
                                                        (*environment)[f].production[index] +=
                                                                patches->recvdata[p][f*(size.x*size.y*size.z)+size.z * size.y * c.x + size.z * c.y + c.z];
                                        }
                                }
                free(patches->recvdata[p]);
        }

        free(patches->recvdata);

        for (p = 0; p < systeminfo.size; p++)
                if (patches->intersect[p])
                        free(patches->data[p]);
        free(patches->data);

        return 0;
}

/*!
 * Now each local cell should receive information about values
 * of global fields in the grid.
 * Field patches are filled with appropriate values.
 * Non-blocking communication is initiated.
 * This is executed after global fields are computed.
 * Field patches buffers are allocated here.
 * Sizes of the patches are the same as those from previous CIC communication.
 * Receiving field patches are also allocated here.
 * MPI_Request tables are allocated here.
 */
void patches_env2cellsinfo(systeminfo_t systeminfo, settings_t settings, patches_t *patches, grid_t *grid,environment_t **environment)
{

        int f; /* fields index */
        int p; /* process index */
        int643dv_t idx, g;
        int643dv_t size;
// UWAGA!!!!
        patches->commbuff = (double **) calloc(systeminfo.size, sizeof(double *));
        for (p = 0; p < systeminfo.size; p++) {
                int fieldsize;
                if (!patches->sender[p])
                        continue; /* continue to next process if current process do not overlap domain */
                size.x = patches->uppercornerR[p].x - patches->lowercornerR[p].x + 1;
                size.y = patches->uppercornerR[p].y - patches->lowercornerR[p].y + 1;
                size.z = patches->uppercornerR[p].z - patches->lowercornerR[p].z + 1;
                fieldsize = size.x * size.y * size.z;
                patches->commbuff[p] =
                        (double *) calloc((settings.numberoffields*4) * fieldsize, sizeof(double));
                /* fields */
                for (f = 0; f < settings.numberoffields; f++) {
                        int64_t i, j, k;
                        for (i = patches->lowercornerR[p].x; i <= patches->uppercornerR[p].x; i++)
                                for (j = patches->lowercornerR[p].y; j <= patches->uppercornerR[p].y; j++)
                                        for (k = patches->lowercornerR[p].z; k <= patches->uppercornerR[p].z; k++) {
                                                int index1,index2;
                                                idx.x = i - patches->lowercornerR[p].x;
                                                idx.y = j - patches->lowercornerR[p].y;
                                                idx.z = k - patches->lowercornerR[p].z;
                                                g.x = i - grid->loweridx[systeminfo.rank].x;
                                                g.y = j - grid->loweridx[systeminfo.rank].y;
                                                g.z = k - grid->loweridx[systeminfo.rank].z;
                                                index1 = f * fieldsize + size.z * size.y * idx.x + size.z * idx.y + idx.z;
                                                index2 = grid->localsize.z * grid->localsize.y * g.x + grid->localsize.z * g.y + g.z;
                                                patches->commbuff[p][index1] = (*environment)[f].data[index2];
                                        }
                }
                /* field gradients */
                for (f = 0; f < settings.numberoffields; f++) {
                        int64_t i, j, k;
                        for (i = patches->lowercornerR[p].x; i <= patches->uppercornerR[p].x; i++)
                                for (j = patches->lowercornerR[p].y; j <= patches->uppercornerR[p].y; j++)
                                        for (k = patches->lowercornerR[p].z; k <= patches->uppercornerR[p].z; k++) {
                                                int index1,index2;
                                                idx.x = i - patches->lowercornerR[p].x;
                                                idx.y = j - patches->lowercornerR[p].y;
                                                idx.z = k - patches->lowercornerR[p].z;
                                                g.x = i - grid->loweridx[systeminfo.rank].x;
                                                g.y = j - grid->loweridx[systeminfo.rank].y;
                                                g.z = k - grid->loweridx[systeminfo.rank].z;
                                                index1 = fieldsize*settings.numberoffields + 3* (f * fieldsize + size.z * size.y * idx.x + size.z * idx.y + idx.z);
                                                index2 = 3*(grid->localsize.z*grid->localsize.y*g.x + grid->localsize.z*g.y + g.z);
                                                patches->commbuff[p][index1] = (*environment)[f].gradient[index2];
                                                patches->commbuff[p][index1+1] = (*environment)[f].gradient[index2+1];
                                                patches->commbuff[p][index1+2] = (*environment)[f].gradient[index2+2];
                                        }
                }
        }
        return;
}

void patches_commenv2cellsinit(systeminfo_t systeminfo, settings_t settings, patches_t *patches){

        int p; /* process index */
        // UWAGA
        if (!(patches->buff = (double **) calloc(systeminfo.size, sizeof(double *))))
                terminate(systeminfo,"cannot allocate patches->buff", __FILE__, __LINE__);

        patches->reqsend = (MPI_Request *) malloc(sizeof(MPI_Request) * systeminfo.size);
        patches->reqrecv = (MPI_Request *) malloc(sizeof(MPI_Request) * systeminfo.size);

        for (p = 0; p < systeminfo.size; p++) {
                if (patches->sender[p]) {
                        int sendSize;
                        sendSize =  (patches->uppercornerR[p].x - patches->lowercornerR[p].x +
                                     1) * (patches->uppercornerR[p].y - patches->lowercornerR[p].y +
                                           1) * (patches->uppercornerR[p].z - patches->lowercornerR[p].z +
                                                 1) * (settings.numberoffields*4);
                        MPI_Isend(&(patches->commbuff[p][0]), sendSize, MPI_DOUBLE, p,
                                  systeminfo.rank, MPI_COMM_WORLD, &(patches->reqsend[p]));
                }
                if (patches->receiver) {
                        int recvSize;
                        recvSize = patches->size[p].x * patches->size[p].y * patches->size[p].z * (settings.numberoffields*4);
                        if (!(patches->buff[p] = (double *) calloc(recvSize, sizeof(double))))
                                terminate(systeminfo,"cannot allocate patches->buff[p]", __FILE__, __LINE__);
                        MPI_Irecv(&(patches->buff[p][0]), recvSize, MPI_DOUBLE, p, p,
                                  MPI_COMM_WORLD, &(patches->reqrecv[p]));
                }
        }
        return;
}

/*!
 * Wait for communication to finish.
 * MPI_Request tables are deallocated here.
 * Field patches buffers are deallocated here.
 */
void patches_commenv2cellswait(systeminfo_t systeminfo, patches_t *patches)
{
        int p;
        MPI_Status status;

        for (p = 0; p < systeminfo.size; p++) {
                if (!patches->sender[p])
                        continue;
                if (MPI_Wait(&patches->reqsend[p], &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication error commenv2cellswait 1", __FILE__, __LINE__);
        }

        for (p = 0; p < systeminfo.size; p++) {
                if (!patches->receiver)
                        continue;
                if (MPI_Wait(&patches->reqrecv[p], &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication error commenv2cellswait 2", __FILE__, __LINE__);
        }

        free(patches->reqsend);
        free(patches->reqrecv);

        for (p = 0; p < systeminfo.size; p++)
                free(patches->commbuff[p]);
        free(patches->commbuff);

        return;
}

/*!
 * Update local cell fields with information from
 * remote processes received in patches.
 * Receiveing field patches are deallocated here.
 */
void patches_applyenv2cells(systeminfo_t systeminfo, settings_t settings, patches_t *patches, cellsinfo_t *cellsinfo,grid_t *grid,cellenvdata_t ***cellenvdata)
{

        int p,c,f;
        int643dv_t idx;
        double3dv_t d, t;
        int643dv_t cellIdx;
        double3dv_t cicCoord;

        /* allocate cellenvdata */
        if(!((*cellenvdata)=(cellenvdata_t**)calloc(settings.numberoffields,sizeof(cellenvdata_t*))))
                terminate(systeminfo,"cannot allocate cellenvdata", __FILE__, __LINE__);
        for(f=0; f<settings.numberoffields; f++)
                if(!((*cellenvdata)[f]=(cellenvdata_t*)calloc(cellsinfo->localcount.n,sizeof(cellenvdata_t))))
                        terminate(systeminfo,"cannot allocate cellenvdata[f]", __FILE__, __LINE__);

        /* reset fields */
        for (f = 0; f < settings.numberoffields; f++)
                for (c = 0; c < cellsinfo->localcount.n; c++) {
                        (*cellenvdata)[f][c].value=0.0;
                        (*cellenvdata)[f][c].gx=0.0;
                        (*cellenvdata)[f][c].gy=0.0;
                        (*cellenvdata)[f][c].gz=0.0;
                }

        for (c = 0; c < cellsinfo->localcount.n; c++) { /* for every cell */
                cellIdx.x = ((cellsinfo->cells[c].x - grid->lowercorner.x) / grid->resolution);
                cellIdx.y = ((cellsinfo->cells[c].y - grid->lowercorner.y) / grid->resolution);
                cellIdx.z = ((cellsinfo->cells[c].z - grid->lowercorner.z) / grid->resolution);
                for (p = 0; p < systeminfo.size; p++) { /* for each process */
                        int ax, ay, az;
                        if (!patches->receiver)
                                continue; /* there is no patch from this process */
                        for (ax = 0; ax < 2; ax++)
                                for (ay = 0; ay < 2; ay++)
                                        for (az = 0; az < 2; az++) {
                                                if (cellIdx.x + ax >= grid->loweridx[p].x
                                                    && cellIdx.y + ay >= grid->loweridx[p].y
                                                    && cellIdx.z + az >= grid->loweridx[p].z
                                                    && cellIdx.x + ax <= grid->upperidx[p].x
                                                    && cellIdx.y + ay <= grid->upperidx[p].y
                                                    && cellIdx.z + az <= grid->upperidx[p].z) {

                                                        idx.x = (cellIdx.x + ax) - patches->lowercorner[p].x;
                                                        idx.y = (cellIdx.y + ay) - patches->lowercorner[p].y;
                                                        idx.z = (cellIdx.z + az) - patches->lowercorner[p].z;

                                                        cicCoord.x = grid->lowercorner.x + cellIdx.x * grid->resolution;
                                                        cicCoord.y = grid->lowercorner.y + cellIdx.y * grid->resolution;
                                                        cicCoord.z = grid->lowercorner.z + cellIdx.z * grid->resolution;

                                                        d.x = (cellsinfo->cells[c].x - cicCoord.x) / grid->resolution;
                                                        d.y = (cellsinfo->cells[c].y - cicCoord.y) / grid->resolution;
                                                        d.z = (cellsinfo->cells[c].z - cicCoord.z) / grid->resolution;

                                                        t.x = 1.0 - d.x;
                                                        t.y = 1.0 - d.y;
                                                        t.z = 1.0 - d.z;

                                                        /* interpolating back to cells */
                                                        /* scaling from mol/cm^3 to mol/cell */
                                                        for (f = 0; f < settings.numberoffields; f++) {
                                                                int index1 = f * patches->size[p].x * patches->size[p].y * patches->size[p].z + patches->size[p].y * patches->size[p].z * idx.x + patches->size[p].z * idx.y + idx.z;
                                                                (*cellenvdata)[f][c].value += patches->buff[p][index1] * (ax * d.x + (1 - ax) * t.x) * (ay * d.y + (1 - ay) * t.y) * (az * d.z + (1 - az) * t.z); //*cellVolume;
                                                        }
                                                        for (f = 0; f < settings.numberoffields; f++) {
                                                                int index1=settings.numberoffields * patches->size[p].x * patches->size[p].y * patches->size[p].z + 3* (f * patches->size[p].x * patches->size[p].y * patches->size[p].z + patches->size[p].y * patches->size[p].z * idx.x + patches->size[p].z * idx.y + idx.z);
                                                                (*cellenvdata)[f][c].gx += patches->buff[p][index1] * (ax * d.x + (1 - ax) * t.x) * (ay * d.y + (1 - ay) * t.y) * (az * d.z + (1 - az) * t.z);
                                                                (*cellenvdata)[f][c].gy += patches->buff[p][index1 + 1] * (ax * d.x + (1 - ax) * t.x) * (ay * d.y + (1 - ay) * t.y) * (az * d.z + (1 - az) * t.z);
                                                                (*cellenvdata)[f][c].gz += patches->buff[p][index1 + 2] * (ax * d.x + (1 - ax) * t.x) * (ay * d.y + (1 - ay) * t.y) * (az * d.z + (1 - az) * t.z);
                                                        }
                                                } // if
                                        }// az
                } // p
        } // c

        for (p = 0; p < systeminfo.size; p++)
                free(patches->buff[p]);
        free(patches->buff);

        return;
}

/*!
 * This function initializes data exchange between
 * processes required in cells-to-grid interpolation.
 * This function enables overlapping communication and
 * computations.
 */
void patches_cells2envinit(systeminfo_t systeminfo, settings_t settings, patches_t *patches, cellsinfo_t *cellsinfo, grid_t *grid)
{
        if (settings.numberoffields==0)
                return;
        patches_cells2envinfo(systeminfo,settings,patches,cellsinfo,grid);
        patches_commcells2envinit(systeminfo,settings,patches);
        return;
}

/*!
 * This function wait for patches communication to finish
 * in cells-to-grid interpolation.
 * This function enables overlapping communication and
 * computations.
 */
void patches_cells2envwait(systeminfo_t systeminfo, settings_t settings, patches_t *patches, grid_t *grid, environment_t **environment)
{
        if (settings.numberoffields==0)
                return;
        patches_commcells2envwait(systeminfo,patches);
        patches_applycells2env(systeminfo,settings,patches,grid,environment);
        return;
}

void patches_free(patches_t *patches) {
        free(patches->receiver);
        free(patches->sender);
        free(patches->intersect);
        free(patches->lowercorner);
        free(patches->uppercorner);
        free(patches->lowercornerR);
        free(patches->uppercornerR);
        free(patches->size);
        return;
}

/*!
 * This function is used to interpolate field data back to
 * cells. Overlapping of communication and computations is not
 * implemented here.
 * This function deallocates all important arrays used in interpolation.
 */
void patches_env2cellsinit(systeminfo_t systeminfo, settings_t settings, patches_t *patches, grid_t *grid, environment_t **environment)
{
        if (settings.numberoffields==0) return;
        patches_env2cellsinfo(systeminfo,settings,patches,grid,environment);
        patches_commenv2cellsinit(systeminfo,settings,patches);
        return;
}

void patches_env2cellswait(systeminfo_t systeminfo, settings_t settings, patches_t *patches,cellsinfo_t *cellsinfo, grid_t *grid,cellenvdata_t ***cellenvdata)
{
        if (settings.numberoffields==0) return;
        patches_commenv2cellswait(systeminfo,patches);
        patches_applyenv2cells(systeminfo,settings,patches,cellsinfo,grid,cellenvdata);
        return;
}
