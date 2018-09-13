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

#include "global.h"
#include "potential.h"
#include "fieldgradient.h"
#include "exchange.h"
#include "patches.h"
#include "environment.h"

/*! \file compute.c
 *  \brief contains main computational function called in each time step of the simulation
 */

/*!
 * This function calls all important simulation steps (cellular dynamics and global fields computations).
 */
int step_compute(systeminfo_t systeminfo, settings_t settings, cellsinfo_t *cellsinfo, celltype_t* celltype, grid_t *grid, environment_t **environment,cellcommdata_t *cellcommdata,interpdata_t *interpdata,cellenvdata_t ***cellenvdata, solverdata_t *solverdata,solversettings_t *solversettings)
{
        int p;
        double sf;
        patches_t patches;
        /* 0. Initialization */
        patches_alloc(systeminfo,settings,&patches,cellsinfo,grid);
        patches_cells2envinit(systeminfo,settings,&patches,cellsinfo,grid);

        /* initiate asynchronous data transfers between processors */
        cellssendrecv(systeminfo,*cellsinfo,cellcommdata);

        /* 1. Compute potential for local cells */

        /* compute potential for local cells */
        computepotential(systeminfo,cellsinfo,celltype,*cellcommdata);

        /* 2. Solve global fields */

        if(settings.step>0) {
                patches_cells2envwait(systeminfo,settings,&patches,grid,environment);
                environment_compute(systeminfo,settings,grid,environment,solverdata,solversettings);
        }

        /* wait for data transfers to finish */
        cellswait(systeminfo,*cellsinfo,cellcommdata);

        /* 3. Compute potential for remote cells */
        computeremotepotential(systeminfo,cellsinfo,celltype,*cellcommdata);

        /* 4. Add chemotactic term to potential */

        /* add chemotaxis term to potential */

        /* 5. Compute gradient of the potential for local cells */
        /* initiate transfer of the density and potential data from remote cells */
        datasendrecv(systeminfo,*cellsinfo,cellcommdata);
        /* compute gradient of the potential for local cells */
        computegradient(systeminfo,cellsinfo,celltype,*cellcommdata);

        /* 6. Interpolate global fields and compute gradient */

        /* interpolate data */
        patches_env2cellsinit(systeminfo,settings,&patches,grid,environment);

        patches_env2cellswait(systeminfo,settings,&patches,cellsinfo,grid,cellenvdata);

        patches_free(&patches);
        /* compute gradient of global fields */
        fieldgradient(systeminfo,settings,environment,grid);

        /* 7. Compute gradient of the potential for remote cells */
        /* wait for density and potential data from remote cells */
        datawait(systeminfo,*cellsinfo,cellcommdata);

        /* compute gradient of the potential for remote cells */
        computeremotegradient(systeminfo,cellsinfo,celltype,*cellcommdata);

        return 0;
}

void step_clean(settings_t settings,cellenvdata_t ***cellenvdata) {
        int f;
        for(f=0; f<settings.numberoffields; f++) {
                free((*cellenvdata)[f]);
        }
        free((*cellenvdata));
}
