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
#include <mpi.h>
#include <omp.h>

#include "global.h"

#include "lb.h"
#include "io.h"
#include "cells.h"
#include "grid.h"
#include "step.h"
#include "utils.h"
#include "octree.h"
#include "exchange.h"
#include "initialisation.h"
#include "environment.h"

/*! \file main.c
 *  \brief contains the main simulation loop
 */

/*!
 * This function intializes MPI, calls Timothy initialization and allocation functions.
 * It also contains the main simulation loop where all important simulation steps are called.
 */

int main(int argc, char **argv)
{
        systeminfo_t systeminfo;
        settings_t settings;
        celltype_t *celltype;
        environment_t *environment;
        cellsinfo_t cellsinfo;
        grid_t grid;
        cellcommdata_t cellcommdata;
        fieldcommdata_t fieldcommdata;
        interpdata_t interpdata;
        statistics_t statistics;
        solverdata_t solverdata;
        solversettings_t solversettings;
        cellenvdata_t **cellenvdata;
        struct Zoltan_Struct *ztn;

        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &systeminfo.size);
        MPI_Comm_rank(MPI_COMM_WORLD, &systeminfo.rank);
        systeminfo.nthreads = omp_get_max_threads();

        getsysteminfo(&systeminfo);
        printinfo(systeminfo);
        initialisation(argc,argv,&systeminfo,&settings,&celltype,&environment);
        allocatecells(systeminfo,settings,celltype,&cellsinfo);
        allocategrid(systeminfo,settings,&grid);
        environment_allocate(systeminfo,settings,grid,&environment,&solverdata,&solversettings);
        environment_init(systeminfo,settings,grid,&environment);
        #ifdef HYPRE
        environment_initsystem(systeminfo,settings,&grid,&environment,&solverdata,&solversettings);
        #endif
        lbinit(argc,argv,MPI_COMM_WORLD,systeminfo,&ztn,&cellsinfo);

        settings.simulationstart=1;

        for (settings.step = 0; settings.step < settings.numberofsteps; settings.step++) {
                updateglobalcounts(&cellsinfo);
                lbexchange(systeminfo,ztn);
                octbuild(systeminfo,&cellsinfo,celltype);
                createexportlist(systeminfo,settings,cellsinfo,grid,ztn,celltype,&cellcommdata,&fieldcommdata);
                step_compute(systeminfo,settings,&cellsinfo,celltype,&grid,&environment,&cellcommdata,&interpdata,&cellenvdata,&solverdata,&solversettings);
                exchangecleanup(systeminfo,cellsinfo,&cellcommdata,&fieldcommdata);
                cellsupdate(systeminfo,settings,celltype,cellenvdata,&cellsinfo);
                if(!(settings.step%settings.statoutstep)) printstatistics(systeminfo,settings,cellsinfo,&statistics);
                if(!(settings.step%settings.visoutstep)) writevtk(systeminfo,settings,cellsinfo);
                octfree(&cellsinfo);
                step_clean(settings,&cellenvdata);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        //cellsdestroy();
        lbdestroy(&ztn);

        if (systeminfo.rank == 0)
                printf("\nEnd of simulation run.\n");

        MPI_Finalize();

        return 0;
}
