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
#include <string.h>
#include <inttypes.h>
#include <mpi.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <float.h>
#include <time.h>

#if defined(__bg__) && defined(__bgq__)
#include <spi/include/kernel/process.h>
#include <spi/include/kernel/location.h>
#endif

#include "global.h"

/*! \file utils.c
 *  \brief contains various utility functions
 */

/*!
 * This function checks the endiannes of the systeminfo
 */
int checkendiannes(systeminfo_t *systeminfo)
{
        volatile uint32_t i = 0x01234567;
        /* return 0 for big endian, 1 for little endian. */
        systeminfo->endian = (*((uint8_t *) (&i))) == 0x67;
        return 0;
}

uint64_t n;
uint64_t g0phase;
uint64_t g1phase;
uint64_t sphase;
uint64_t g2phase;
uint64_t mphase;
uint64_t necroticphase;

/*!
 * This function swaps endiannes within a table of n elements
 * of a given size m (given in bytes).
 */
void swapnbyte(char *data, int n, int m)
{
        int i, j;
        char old_data[16];

        for (j = 0; j < n; j++) {
                memcpy(&old_data[0], &data[j * m], m);
                for (i = 0; i < m; i++)
                        data[j * m + i] = old_data[m - i - 1];
        }
}

void updateglobalcounts(cellsinfo_t* cellsinfo){
        uint64_t localcount[8];
        uint64_t globalcount[8];
        localcount[0]=cellsinfo->localcount.n;
        localcount[1]=cellsinfo->localcount.g0phase;
        localcount[2]=cellsinfo->localcount.g1phase;
        localcount[3]=cellsinfo->localcount.sphase;
        localcount[4]=cellsinfo->localcount.g2phase;
        localcount[5]=cellsinfo->localcount.mphase;
        localcount[6]=cellsinfo->localcount.necroticphase;
        localcount[7]=0;
        MPI_Allreduce(localcount,globalcount,8,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
        cellsinfo->globalcount.n=globalcount[0];
        cellsinfo->globalcount.g0phase=globalcount[1];
        cellsinfo->globalcount.g1phase=globalcount[2];
        cellsinfo->globalcount.sphase=globalcount[3];
        cellsinfo->globalcount.g2phase=globalcount[4];
        cellsinfo->globalcount.mphase=globalcount[5];
        cellsinfo->globalcount.necroticphase=globalcount[6];
        MPI_Allgather(&(cellsinfo->localcount.n),1,MPI_UINT64_T,cellsinfo->cellsperproc,1,MPI_UINT64_T,MPI_COMM_WORLD);
        return;
}

void terminate(systeminfo_t systeminfo, char *msg, char *file, int line) {
        if(systeminfo.rank==0) {
                fprintf(stderr,"error: %s (%s, line %d)\n",msg,file,line);
                fflush(stderr);
        }
        MPI_Abort(MPI_COMM_WORLD,997);
}

/*!
 * This function detects the amount of addressable memory available for each process.
 * Assumptions: each process allocates similar amount of data.
 * This functions is systeminfo dependent. Supported platforms: linux, AIX, Blue Gene/Q.
 */
size_t getMemoryPerProcess(int32_t lsize)
{
#if defined(__bg__) && defined(__bgq__)
        Personality_t pers;
        int32_t msize;
        Kernel_GetPersonality(&pers, sizeof(pers));
        msize = pers.DDR_Config.DDRSizeMB;
#else
#if defined(_AIX)
        int msize = sysconf(_SC_AIX_REALMEM) / 1024;
#else       /* aix */
        long psize = sysconf(_SC_PAGE_SIZE);
        long npages = sysconf(_SC_PHYS_PAGES);
        long msize = psize * npages / (1024 * 1024);
#endif        /* other then aix and bgq - linux or osx */
#endif        /* bgq */
        return msize / lsize;
}

/*!
 * This function detects number of processes on each node and assignes
 * local node ranks for each process.
 * Assumptions: number of processes per node is equal accross nodes.
 * This functions is systeminfo dependent. Supported platforms: linux, AIX, Blue Gene/Q.
 */
void getLocalRankAndSize(int rank, int size, int32_t * lrank,
                         int32_t * lsize)
{
        int i;
        int32_t r, s;
#if defined(__bg__) && defined(__bgq__)
        s = Kernel_ProcessCount();
        r = Kernel_MyTcoord();
#else
        int pnamelen;
        char pname[MPI_MAX_PROCESSOR_NAME];
        char pnametable[size][MPI_MAX_PROCESSOR_NAME];
        MPI_Get_processor_name(pname, &pnamelen);
        MPI_Allgather(pname, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, pnametable,
                      MPI_MAX_PROCESSOR_NAME, MPI_CHAR, MPI_COMM_WORLD);
        r = 0;
        s = 0;
        for (i = 0; i < size; i++)
                if (!strcmp(pnametable[rank], pnametable[i])) {
                        if (i < rank)
                                r++;
                        s++;
                }
#endif
        *lrank = r;
        *lsize = s;
}

void randomstreaminit(systeminfo_t *systeminfo,settings_t *settings)
{
        srand(time(NULL) + systeminfo->rank);
        //settings->rseed=time(NULL) + systeminfo->rank;
        return;
}

void printstatistics(systeminfo_t systeminfo,settings_t settings,cellsinfo_t cellsinfo,statistics_t* statistics)
{
        int p;
        statistics->mindist=DBL_MAX;
        statistics->minspeed=DBL_MAX;
        statistics->maxspeed=0.0;
        statistics->minsize=DBL_MAX;
        statistics->maxsize=0.0;
        statistics->mindens=DBL_MAX;
        statistics->maxdens=0.0;
        for (p = 0; p < cellsinfo.localcount.n; p++) {
                double speed;
                speed =  sqrt(cellsinfo.forces[p].x * cellsinfo.forces[p].x +
                              cellsinfo.forces[p].y * cellsinfo.forces[p].y +
                              cellsinfo.forces[p].z * cellsinfo.forces[p].z);
                statistics->minspeed = (speed < statistics->minspeed ? speed : statistics->minspeed);
                statistics->maxspeed = (speed > statistics->maxspeed ? speed : statistics->maxspeed);
                statistics->minsize = (cellsinfo.cells[p].size < statistics->minsize ? cellsinfo.cells[p].size : statistics->minsize);
                statistics->maxsize = (cellsinfo.cells[p].size > statistics->maxsize ? cellsinfo.cells[p].size : statistics->maxsize);
                statistics->mindist = (cellsinfo.cells[p].mindist < statistics->mindist ? cellsinfo.cells[p].mindist : statistics->mindist);
                statistics->mindens = (cellsinfo.cells[p].density < statistics->mindens ? cellsinfo.cells[p].density : statistics->mindens);
                statistics->maxdens = (cellsinfo.cells[p].density > statistics->maxdens ? cellsinfo.cells[p].density : statistics->maxdens);
        }

        MPI_Allreduce(MPI_IN_PLACE,&(statistics->minspeed), 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&(statistics->maxspeed), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&(statistics->minsize), 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&(statistics->maxsize), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&(statistics->mindist), 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

        if(systeminfo.rank==0) {
                printf("\n+++ simulation step %12d\n",settings.step);
                printf("%12s%10s%10s\n", "", "min", "max");
                printf("%12s%10.4lf%10.4lf\n", "size       ",statistics->minsize,statistics->maxsize);
                printf("%12s%10.4lf%10.4lf\n", "density    ",statistics->mindens,statistics->maxdens);
                printf("%12s%10.4lf%10.4lf\n", "speed      ",statistics->minspeed,statistics->maxspeed);
                printf("%12s%10.4lf%10s\n", "distance   ", statistics->mindist, "N/A");
        }

        return;
}
