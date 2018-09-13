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
#include <math.h>
#include <inttypes.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "global.h"
#include "inline.h"
#include "initialisation.h"
#include "utils.h"
#include "io.h"

/*! \file init.c
 *  \brief contains initialization functions
 */

void getsysteminfo(systeminfo_t* systeminfo) {
        checkendiannes(systeminfo);
        getLocalRankAndSize(systeminfo->rank, systeminfo->size, &(systeminfo->noderank), &(systeminfo->nodesize));
        systeminfo->memperproc = getMemoryPerProcess(systeminfo->nodesize);
        if (!POWER_OF_TWO(systeminfo->size))
                terminate(*systeminfo,"number of processes must be power of two", __FILE__, __LINE__);
        return;
}

void initialsettings(settings_t* settings){
        settings->maxcells=0;
        settings->numberofsteps=0;
        settings->secondsperstep=0;
        settings->numberofcelltypes=0;
        settings->numberoffields=0;
        settings->dimension=3;
        settings->restart=0;
        strcpy(settings->rstfilename,"restart.bin");
        strcpy(settings->outdir,"results");
        settings->visoutstep=0;
        settings->statoutstep=0;
        settings->rstoutstep=0;
        settings->maxspeed=0;
        settings->gfdt=0;
        settings->gfh=0;
        settings->simulationstart=0;
        return;
}

void initialcelltype(int numberofcelltypes,int numberoffields,celltype_t* celltype){
        int i,j;
        for(i=0; i<numberofcelltypes; i++) {
                sprintf(celltype->name,"celltype%d",i);
                celltype[i].g1=CELLTYPE_G1_DEFAULT;
                celltype[i].s=CELLTYPE_S_DEFAULT;
                celltype[i].g2=CELLTYPE_G2_DEFAULT;
                celltype[i].m=CELLTYPE_M_DEFAULT;
                celltype[i].v=CELLTYPE_V_DEFAULT;
                celltype[i].rd=CELLTYPE_RD_DEFAULT;
                celltype[i].criticaldensity=CELLTYPE_CDENS_DEFAULT;
                celltype[i].size=CELLTYPE_SIZE_DEFAULT;
                celltype[i].h=CELLTYPE_H_DEFAULT;
                celltype[i].h2=(celltype[i].h)*(celltype[i].h);
                celltype[i].h3=(celltype[i].h2)*(celltype[i].h);
                celltype[i].h4=(celltype[i].h3)*(celltype[i].h);
        }
        for(i=0; i<numberofcelltypes; i++) {
                for(j=0; j<numberoffields; j++) {
                        celltype[i].production[j]=CELLTYPE_PROD_DEFAULT;
                        celltype[i].consumption[j]=CELLTYPE_CONS_DEFAULT;
                        celltype[i].criticallevel1[j]=CELLTYPE_CL1_DEFAULT;
                        celltype[i].criticallevel2[j]=CELLTYPE_CL2_DEFAULT;
                }
        }
        return;
}

void initialfields(int numberoffields,environment_t* environment) {
        int i,j;
        for(i=0; i<numberoffields; i++) {
                sprintf(environment[i].name,"environment%d",i);
                environment[i].diffusioncoefficient=ENVIRONMENT_DC_DEFAULT;
                environment[i].boundarycondition=ENVIRONMENT_BC_DEFAULT;
                environment[i].initialconditionmean=ENVIRONMENT_ICMEAN_DEFAULT;
                environment[i].initialconditionvariance=ENVIRONMENT_ICVAR_DEFAULT;
                environment[i].lambdadelay=ENVIRONMENT_LAMBDA_DEFAULT;
        }
        return;
}

void initialisation(int argc, char **argv, systeminfo_t *systeminfo, settings_t* settings,celltype_t** celltype,environment_t** environment) {
        int i;
        int periods[3];
        int reorder;
        int maxlocalcells;
        struct stat st = {0};

        if (argc < 2 || argc >2) {
                if(systeminfo->rank==0) { printf("usage: timothy <parameter file>\n"); fflush(stdout); }
                MPI_Abort(MPI_COMM_WORLD,-1);
        }
        initialsettings(settings);
        readparamfile(argc,argv,*systeminfo,settings);
        settings->step=0;

        if (settings->numberofcelltypes<1)
                terminate(*systeminfo,"no cell types specified", __FILE__, __LINE__);

        if(!(*celltype=(celltype_t*)malloc((settings->numberofcelltypes)*sizeof(celltype_t))))
                terminate(*systeminfo,"cannot allocate celltype", __FILE__, __LINE__);
        if(!(*environment=(environment_t*)malloc((settings->numberoffields)*sizeof(environment_t))))
                terminate(*systeminfo,"cannot allocate environment", __FILE__, __LINE__);

        for(i=0; i<settings->numberofcelltypes; i++) {
                int size=settings->numberoffields*sizeof(float);
                (*celltype)[i].production=(float*)malloc(size);
                (*celltype)[i].consumption=(float*)malloc(size);
                (*celltype)[i].criticallevel1=(float*)malloc(size);
                (*celltype)[i].criticallevel2=(float*)malloc(size);
        }
        initialcelltype(settings->numberofcelltypes,settings->numberoffields,*celltype);
        readcellsfile(*systeminfo,settings,*celltype);
        initialfields(settings->numberoffields,*environment);
        readenvfile(*systeminfo,settings,*environment);

        randomstreaminit(systeminfo,settings);

        /* organizing processes in a Cartesian grid for global fields computations */
        MPI_Dims_create(systeminfo->size, settings->dimension, systeminfo->dim);
        periods[0] = 0;
        periods[1] = 0;
        periods[2] = 0;
        reorder = 0;
        MPI_Cart_create(MPI_COMM_WORLD, settings->dimension, systeminfo->dim, periods, reorder,
                        &(systeminfo->MPI_CART_COMM));

        if(!(systeminfo->coords = (int **) malloc(systeminfo->size * sizeof(int *))))
                terminate(*systeminfo,"cannot allocate systeminfo->coords", __FILE__, __LINE__);
        for (i = 0; i < systeminfo->size; i++) {
                if(!(systeminfo->coords[i] = (int *) malloc(3 * sizeof(int))))
                        terminate(*systeminfo,"cannot allocate systeminfo->coords[i]", __FILE__, __LINE__);
                MPI_Cart_coords(systeminfo->MPI_CART_COMM, i, settings->dimension, systeminfo->coords[i]);
        }

        maxlocalcells=settings->maxcells / systeminfo->size;
        if (systeminfo->rank < maxlocalcells % systeminfo->size)
                maxlocalcells++;

        settings->maxlocalcells=maxlocalcells;

        /* create output directories */
        if (stat("vtk", &st) == -1) {
                mkdir("vtk", 0755);
        }
        if (stat("rst", &st) == -1) {
                mkdir("rst", 0755);
        }

        /* define colormaps for PovRay outputs */
        /* defineColormaps(); */

        /* density critical levels (very important parameters) */
        /*      if (sdim == 3) {
                      densityCriticalLevel1 = 6 * h3 * sph_kernel(1.6 * csize); //1.8  //1.4 //1.75
                      densityCriticalLevel2 = 6 * h3 * sph_kernel(1.1 * csize); //1.1 //1.4
              }
              if (sdim == 2) {
                      densityCriticalLevel1 = 4 * h2 * sph_kernel(1.4 * csize); //1.4 //1.75
                      densityCriticalLevel2 = 4 * h2 * sph_kernel(1.15 * csize); //1.1 //1.4
              }*/


        return;
}

void initcount(cellcount_t *cellcount) {
        cellcount->n=0;
        cellcount->g0phase=0;
        cellcount->g1phase=0;
        cellcount->sphase=0;
        cellcount->g2phase=0;
        cellcount->mphase=0;
        cellcount->necroticphase=0;
        return;
}

void allocatecells(systeminfo_t systeminfo,settings_t settings,celltype_t *celltype,cellsinfo_t *cellsinfo) {
        int i;
        int maxlocalcells;

        maxlocalcells=settings.maxlocalcells;

        if(!(cellsinfo->cells=(celldata_t*)malloc(maxlocalcells*sizeof(celldata_t))))
                terminate(systeminfo,"cannot allocate cellsinfo->cells", __FILE__, __LINE__);

        if(!(cellsinfo->forces=(double3dv_t*)calloc(maxlocalcells,sizeof(double3dv_t))))
                terminate(systeminfo,"cannot allocate cellsinfo->forces", __FILE__, __LINE__);

        if(!(cellsinfo->cellsperproc=(uint64_t*)malloc(sizeof(uint64_t)*systeminfo.size)))
                terminate(systeminfo,"cannot allocate cellsinfo->cellsperproc", __FILE__, __LINE__);

        if(!(cellsinfo->localtypecount=(cellcount_t*)malloc(sizeof(cellcount_t)*settings.numberofcelltypes)))
                terminate(systeminfo,"cannot allocate cellsinfo->localtypecount", __FILE__, __LINE__);

        if(!(cellsinfo->globaltypecount=(cellcount_t*)malloc(sizeof(cellcount_t)*settings.numberofcelltypes)))
                terminate(systeminfo,"cannot allocate cellsinfo->globaltypecount", __FILE__, __LINE__);

        initcount(&(cellsinfo->localcount));
        initcount(&(cellsinfo->globalcount));
        for(i=0; i<settings.numberofcelltypes; i++) {
                initcount(&(cellsinfo->localtypecount[i]));
                initcount(&(cellsinfo->globaltypecount[i]));
        }

        for(i=0; i<systeminfo.size; i++)
                cellsinfo->cellsperproc[i]=0;

        cellsinfo->dimension=settings.dimension;

        if(settings.restart) {
                // read restart file
        } else {
                readcellpositions(systeminfo,settings,celltype,cellsinfo);
        }

        return;
}

void printinfo(systeminfo_t systeminfo) {
        if (systeminfo.rank == 0) {
                printf("\ntimothy, tissue modelling framework\n");
                printf("http://tissuemodelling.github.io/timothy\n");
                printf("version %s\n\n", VERSION);
                printf("number of MPI processes: %d\n",systeminfo.size);
                printf("processes per node: %d\n",systeminfo.nodesize);
                printf("threads per process: %d\n",systeminfo.nthreads);
                printf("systeminfo: ");
                if (systeminfo.endian)
                        printf("%s, little endian\n", CPUARCH);
                else
                        printf("%s, big endian\n", CPUARCH);
                printf("\n");
                fflush(stdout);
        }
        return;
}
