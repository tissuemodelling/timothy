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

#include <float.h>
#include <inttypes.h>
#include "_hypre_utilities.h"
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_krylov.h"
#include "global.h"
#include "utils.h"
#include "mpi.h"

void environment_allocate(systeminfo_t systeminfo,settings_t settings,grid_t grid,environment_t **environment,solverdata_t *solverdata,solversettings_t *solversettings) {
        int i;

        if(settings.numberoffields==0) return;

         #ifdef HYPRE
        solversettings->dt=(double*)calloc(settings.numberoffields,sizeof(double));
        solversettings->z=(double*)calloc(settings.numberoffields,sizeof(double));
        solversettings->vartypes=(HYPRE_SStructVariable*)calloc(settings.numberoffields,sizeof(HYPRE_SStructVariable));
        solversettings->stencil=(HYPRE_SStructStencil*)calloc(settings.numberoffields,sizeof(HYPRE_SStructStencil));
         #endif
        for(i=0; i<settings.numberoffields; i++)
                if(!((*environment)[i].data=(double*) calloc(grid.localsize.x*grid.localsize.y*grid.localsize.z,sizeof(double))))
                        terminate(systeminfo,"cannot allocate environment->data", __FILE__, __LINE__);
        for(i=0; i<settings.numberoffields; i++)
                if(!((*environment)[i].production=(double*) calloc(grid.localsize.x*grid.localsize.y*grid.localsize.z,sizeof(double))))
                        terminate(systeminfo,"cannot allocate environment->production", __FILE__, __LINE__);
        for(i=0; i<settings.numberoffields; i++)
                if(!((*environment)[i].gradient=(double*) calloc(grid.localsize.x*grid.localsize.y*grid.localsize.z*3,sizeof(double))))
                        terminate(systeminfo,"cannot allocate environment->gradient", __FILE__, __LINE__);
        return;
}

void environment_init(systeminfo_t systeminfo,settings_t settings,grid_t grid,environment_t **environment) {
        int i,j,k,f;

        if(settings.numberoffields==0) return;

        int yz=grid.localsize.y * grid.localsize.z;
        for(f=0; f<settings.numberoffields; f++) {
                for (i = 0; i < grid.localsize.x; i++)
                        for (j = 0; j < grid.localsize.y; j++)
                                for (k = 0; k < grid.localsize.z; k++) {
                                        (*environment)[f].data[yz * i + grid.localsize.z * j + k] =
                                                (*environment)[f].initialconditionmean + ((double)(rand())/RAND_MAX)*(*environment)[f].initialconditionvariance;
                                }
        }
        return;
}

/*!
 * This function sets boundary conditions for domain faces.
 */
void environment_setboundary(int coord, int boundary, solversettings_t *solversettings)
{
        if (coord == 0 && boundary == -1) {
                solversettings->bclower[0] = solversettings->lower[0];
                solversettings->bcupper[0] = solversettings->lower[0];
                solversettings->bclower[1] = solversettings->lower[1];
                solversettings->bcupper[1] = solversettings->upper[1];
                solversettings->bclower[2] = solversettings->lower[2];
                solversettings->bcupper[2] = solversettings->upper[2];
        }
        if (coord == 0 && boundary == 1) {
                solversettings->bclower[0] = solversettings->upper[0];
                solversettings->bcupper[0] = solversettings->upper[0];
                solversettings->bclower[1] = solversettings->lower[1];
                solversettings->bcupper[1] = solversettings->upper[1];
                solversettings->bclower[2] = solversettings->lower[2];
                solversettings->bcupper[2] = solversettings->upper[2];
        }
        if (coord == 1 && boundary == -1) {
                solversettings->bclower[0] = solversettings->lower[0];
                solversettings->bcupper[0] = solversettings->upper[0];
                solversettings->bclower[1] = solversettings->lower[1];
                solversettings->bcupper[1] = solversettings->lower[1];
                solversettings->bclower[2] = solversettings->lower[2];
                solversettings->bcupper[2] = solversettings->upper[2];
        }
        if (coord == 1 && boundary == 1) {
                solversettings->bclower[0] = solversettings->lower[0];
                solversettings->bcupper[0] = solversettings->upper[0];
                solversettings->bclower[1] = solversettings->upper[1];
                solversettings->bcupper[1] = solversettings->upper[1];
                solversettings->bclower[2] = solversettings->lower[2];
                solversettings->bcupper[2] = solversettings->upper[2];
        }
        if (coord == 2 && boundary == -1) {
                solversettings->bclower[0] = solversettings->lower[0];
                solversettings->bcupper[0] = solversettings->upper[0];
                solversettings->bclower[1] = solversettings->lower[1];
                solversettings->bcupper[1] = solversettings->upper[1];
                solversettings->bclower[2] = solversettings->lower[2];
                solversettings->bcupper[2] = solversettings->lower[2];
        }
        if (coord == 2 && boundary == 1) {
                solversettings->bclower[0] = solversettings->lower[0];
                solversettings->bcupper[0] = solversettings->upper[0];
                solversettings->bclower[1] = solversettings->lower[1];
                solversettings->bcupper[1] = solversettings->upper[1];
                solversettings->bclower[2] = solversettings->upper[2];
                solversettings->bcupper[2] = solversettings->upper[2];
        }
}

/*!
 * This function initializes grid, stencil and matrix for a given envical field.
 */
void environment_initsystem(systeminfo_t systeminfo,settings_t settings,grid_t *grid,environment_t **environment,solverdata_t *solverdata,solversettings_t *solversettings)
{
        int i, j, k, c;
        int entry;
        int var;

        if(settings.numberoffields==0) return;

        HYPRE_Int offsets[7][3] = {
                {0, 0, 0}, {-1, 0, 0}, {1, 0, 0}, {0, -1, 0},
                {0, 1, 0}, { 0, 0,-1}, {0, 0, 1}
        };

        if((settings.secondsperstep/settings.gfdt)>1) {
                int intdiv;
                float a1,a2;
                intdiv=(int)(settings.secondsperstep/settings.gfdt)+1;
                a1=settings.secondsperstep/intdiv;
                a2=settings.secondsperstep/(intdiv-1);
                if((settings.gfdt-a1)<(a2-settings.gfdt)) {
                        solversettings->numberofiters = intdiv;
                        for(var=0; var<settings.numberoffields; var++)
                                solversettings->dt[var]=a1;
                } else {
                        solversettings->numberofiters = intdiv-1;
                        for(var=0; var<settings.numberoffields; var++)
                                solversettings->dt[var]=a2;
                }
        } else {
                solversettings->numberofiters = 1;
                for(var=0; var<settings.numberoffields; var++)
                        solversettings->dt[var]=settings.secondsperstep;

        }

        solversettings->enviter=0;

        for(var=0; var<settings.numberoffields; var++)
                solversettings->vartypes[var]=HYPRE_SSTRUCT_VARIABLE_NODE;


        /* 1. INIT GRID */

        /* create an empty 3D grid object */
        HYPRE_SStructGridCreate(systeminfo.MPI_CART_COMM, 3, 1, &(solversettings->grid));

        /* set this process box */
        solversettings->lower[0] = grid->loweridx[systeminfo.rank].x;
        solversettings->lower[1] = grid->loweridx[systeminfo.rank].y;
        solversettings->lower[2] = grid->loweridx[systeminfo.rank].z;

        solversettings->upper[0] = grid->upperidx[systeminfo.rank].x;
        solversettings->upper[1] = grid->upperidx[systeminfo.rank].y;
        solversettings->upper[2] = grid->upperidx[systeminfo.rank].z;

        /* add a new box to the grid */
        HYPRE_SStructGridSetExtents(solversettings->grid, 0, solversettings->lower, solversettings->upper);

        HYPRE_SStructGridSetVariables(solversettings->grid, 0, settings.numberoffields, solversettings->vartypes);
        HYPRE_SStructGridAssemble(solversettings->grid);

        //  2. INIT STENCIL
        //HYPRE_SStructStencilCreate(3, 7, &solverdata->stencil);
        for(var = 0; var < settings.numberoffields; var++) {
                HYPRE_SStructStencilCreate(3, 7, &solversettings->stencil[var]);
                for(entry = 0; entry < 7; entry++)
                        HYPRE_SStructStencilSetEntry(solversettings->stencil[var], entry, offsets[entry],var);

        }

        // 3. SET UP THE GRAPH
        // assumption - all stencils are the same
        solversettings->envobjecttype = HYPRE_PARCSR;
        HYPRE_SStructGraphCreate(systeminfo.MPI_CART_COMM, solversettings->grid, &(solversettings->graph));
        HYPRE_SStructGraphSetObjectType(solversettings->graph, solversettings->envobjecttype);
        for(var=0; var<settings.numberoffields; var++)
                HYPRE_SStructGraphSetStencil(solversettings->graph, 0, var, solversettings->stencil[var]);
        HYPRE_SStructGraphAssemble(solversettings->graph);

        // 4. SET UP MATRIX
        long long nentries = 7;
        long long nvalues;
        double *values;
        HYPRE_Int stencil_indices[7];

        nvalues = nentries * grid->localsize.x * grid->localsize.y * grid->localsize.z;
        // create an empty matrix object
        HYPRE_SStructMatrixCreate(systeminfo.MPI_CART_COMM, solversettings->graph, &(solverdata->A));
        HYPRE_SStructMatrixSetObjectType(solverdata->A, solversettings->envobjecttype);
        // indicate that the matrix coefficients are ready to be set
        HYPRE_SStructMatrixInitialize(solverdata->A);

        values = calloc(nvalues, sizeof(double));

        for (j = 0; j < nentries; j++)
                stencil_indices[j] = j;

        for(var=0; var<settings.numberoffields; var++) {
                solversettings->dt[var] = settings.gfdt;
                solversettings->z[var] =
                        (*environment)[var].diffusioncoefficient * solversettings->dt[var] / (grid->resolution * grid->resolution);

                // set the standard stencil at each grid point,
                //  we will fix the boundaries later
                for (i = 0; i < nvalues; i += nentries) {
                        values[i] = 1 + 6.0 * solversettings->z[var] + (*environment)[var].lambdadelay * solversettings->dt[var];
                        for (j = 1; j < nentries; j++)
                                values[i + j] = -solversettings->z[var];
                }

                HYPRE_SStructMatrixSetBoxValues(solverdata->A, 0, solversettings->lower, solversettings->upper, var,
                                                nentries, stencil_indices, values);
        }

        free(values);
}

/*!
 * This function computes cell production/consumption function based on
 * the interpolated cell density field.
 */
void environment_pcfunction(systeminfo_t systeminfo,settings_t settings, int var, grid_t *grid, double *pc)
{
        int i, j, k;

        if (settings.step == 0)
                return;

        int idx = 0;
        for (k = 0; k < grid->localsize.z; k++)
                for (j = 0; j < grid->localsize.y; j++)
                        for (i = 0; i < grid->localsize.x; i++, idx++) {
                                if(i==5 && j==5 && k==5) pc[idx]=100; else
                                        pc[idx] = 0.0; //-fieldConsumption[var] * 0.1 * dt[var];
                                //envPC[idx] = -fieldConsumption[nch] * tissueField[gridSize.z * gridSize.y * i + gridSize.z * j + k] * dt[nch];
                                //   if(bvsim) envPC[idx]+=fieldProduction[nch] * vesselField[gridSize.z * gridSize.y * i + gridSize.z * j +k] * dt[nch] ;
                                //*(cellVolume/boxVolume);//*(1.0/cellVolume);//*dt[nch];//*dt[nch];
                        }
}

/*!
 * This function initializes boundary conditions for a given envical field.
 */
void environment_initboundary(systeminfo_t systeminfo,settings_t settings,grid_t *grid,environment_t **environment,solverdata_t *solverdata,solversettings_t *solversettings)
//struct state simstate,int settings.numberoffields,struct environment *nutrient,struct gridData grid)
{
        int i, j, k;
        int mi;
        int var;
        int nentries = 1;
        HYPRE_Int stencil_indices[1];
        long long nvalues = grid->localsize.x * grid->localsize.y * grid->localsize.z;
        double *values, *bvalues;
        double *pc;

        pc = (double *) calloc(nvalues, sizeof(double));
        values = calloc(nvalues, sizeof(double));
        bvalues = calloc(nvalues, sizeof(double));

        // 5. SETUP STRUCT VECTORS FOR B AND X

        // create an empty vector object
        HYPRE_SStructVectorCreate(systeminfo.MPI_CART_COMM, solversettings->grid, &(solverdata->b));
        HYPRE_SStructVectorCreate(systeminfo.MPI_CART_COMM, solversettings->grid, &(solverdata->x));

        // as with the matrix, set the appropriate object type for the vectors
        HYPRE_SStructVectorSetObjectType(solverdata->b,solversettings->envobjecttype);
        HYPRE_SStructVectorSetObjectType(solverdata->x,solversettings->envobjecttype);

        // indicate that the vector coefficients are ready to be set
        HYPRE_SStructVectorInitialize(solverdata->b);
        HYPRE_SStructVectorInitialize(solverdata->x);

        for(var=0; var<settings.numberoffields; var++) {

                environment_pcfunction(systeminfo,settings,var,grid,pc);

                // set the values
                mi = 0;
                for (k = solversettings->lower[2]; k <= solversettings->upper[2]; k++)
                        for (j = solversettings->lower[1]; j <= solversettings->upper[1]; j++)
                                for (i = solversettings->lower[0]; i <= solversettings->upper[0]; i++) {
                                        values[mi] = (*environment)[var].initialconditionmean;
                                        mi++;
                                }

                HYPRE_SStructVectorSetBoxValues(solverdata->b, 0, solversettings->lower, solversettings->upper, var, values);

                mi = 0;
                for (k = solversettings->lower[2]; k <= solversettings->upper[2]; k++)
                        for (j = solversettings->lower[1]; j <= solversettings->upper[1]; j++)
                                for (i = solversettings->lower[0]; i <= solversettings->upper[0]; i++) {
                                        values[mi] = (*environment)[var].initialconditionmean;
                                        mi++;
                                }

                HYPRE_SStructVectorSetBoxValues(solverdata->x, 0, solversettings->lower, solversettings->upper, var, values);

                // incorporate boundary conditions; Dirichlet on 6 faces

                for (i = 0; i < nvalues; i++)
                        values[i] = solversettings->z[var];
                for (i = 0; i < nvalues; i++)
                        bvalues[i] = solversettings->z[var] * (*environment)[var].boundarycondition;

                if (systeminfo.coords[systeminfo.rank][0] == 0) {
                        nvalues = nentries * grid->localsize.y * grid->localsize.z;
                        environment_setboundary(0, -1, solversettings);
                        stencil_indices[0] = 1;
                        HYPRE_SStructMatrixAddToBoxValues(solverdata->A, 0, solversettings->bclower, solversettings->bcupper, var, nentries, stencil_indices, values);
                        HYPRE_SStructVectorAddToBoxValues(solverdata->b, 0, solversettings->bclower, solversettings->bcupper, var, bvalues);
                }
                if (systeminfo.coords[systeminfo.rank][0] == systeminfo.dim[0] - 1) {
                        nvalues = nentries * grid->localsize.y * grid->localsize.z;
                        environment_setboundary(0, 1, solversettings);
                        stencil_indices[0] = 2;
                        HYPRE_SStructMatrixAddToBoxValues(solverdata->A, 0, solversettings->bclower, solversettings->bcupper, var, nentries, stencil_indices, values);
                        HYPRE_SStructVectorAddToBoxValues(solverdata->b, 0, solversettings->bclower, solversettings->bcupper, var, bvalues);
                }
                if (systeminfo.coords[systeminfo.rank][1] == 0) {
                        nvalues = nentries * grid->localsize.x * grid->localsize.z;
                        environment_setboundary(1, -1, solversettings);
                        stencil_indices[0] = 3;
                        HYPRE_SStructMatrixAddToBoxValues(solverdata->A, 0, solversettings->bclower, solversettings->bcupper, var, nentries, stencil_indices, values);
                        HYPRE_SStructVectorAddToBoxValues(solverdata->b, 0, solversettings->bclower, solversettings->bcupper, var, bvalues);
                }
                if (systeminfo.coords[systeminfo.rank][1] == systeminfo.dim[1] - 1) {
                        nvalues = nentries * grid->localsize.x * grid->localsize.z;
                        environment_setboundary(1, 1, solversettings);
                        stencil_indices[0] = 4;
                        HYPRE_SStructMatrixAddToBoxValues(solverdata->A, 0, solversettings->bclower, solversettings->bcupper, var, nentries, stencil_indices, values);
                        HYPRE_SStructVectorAddToBoxValues(solverdata->b, 0, solversettings->bclower, solversettings->bcupper, var, bvalues);
                }
                if (systeminfo.coords[systeminfo.rank][2] == 0) {
                        nvalues = nentries * grid->localsize.x * grid->localsize.y;
                        environment_setboundary(2, -1, solversettings);
                        stencil_indices[0] = 5;
                        HYPRE_SStructMatrixAddToBoxValues(solverdata->A, 0, solversettings->bclower, solversettings->bcupper, var, nentries, stencil_indices, values);
                        HYPRE_SStructVectorAddToBoxValues(solverdata->b, 0, solversettings->bclower, solversettings->bcupper, var, bvalues);
                }
                if (systeminfo.coords[systeminfo.rank][2] == systeminfo.dim[2] - 1) {
                        nvalues = nentries * grid->localsize.x * grid->localsize.y;
                        environment_setboundary(2, 1, solversettings);
                        stencil_indices[0] = 6;
                        HYPRE_SStructMatrixAddToBoxValues(solverdata->A, 0, solversettings->bclower, solversettings->bcupper, var, nentries, stencil_indices, values);
                        HYPRE_SStructVectorAddToBoxValues(solverdata->b, 0, solversettings->bclower, solversettings->bcupper, var, bvalues);
                }

                // add production consumption function to the right side
                HYPRE_SStructVectorAddToBoxValues(solverdata->b, 0, solversettings->lower, solversettings->upper, var, pc);
        }

        free(pc);
        free(values);
        free(bvalues);
        return;
}

/*!
 * This function initializes Hypre for solving a given envical field.
 */
void environment_initsolver(systeminfo_t systeminfo, solverdata_t *solverdata,solversettings_t *solversettings)
{

        HYPRE_SStructMatrixAssemble(solverdata->A);
        // This is a collective call finalizing the vector assembly.
        //   The vector is now ``ready to be used''
        HYPRE_SStructVectorAssemble(solverdata->b);
        HYPRE_SStructVectorAssemble(solverdata->x);

        HYPRE_SStructMatrixGetObject(solverdata->A, (void **) &(solverdata->parA));
        HYPRE_SStructVectorGetObject(solverdata->b, (void **) &(solverdata->parb));
        HYPRE_SStructVectorGetObject(solverdata->x, (void **) &(solverdata->parx));

        HYPRE_ParCSRPCGCreate(systeminfo.MPI_CART_COMM, &(solverdata->solver));
        HYPRE_ParCSRPCGSetTol(solverdata->solver, 1.0e-12);
        HYPRE_ParCSRPCGSetPrintLevel(solverdata->solver, 2);
        HYPRE_ParCSRPCGSetMaxIter(solverdata->solver, 50);

        HYPRE_BoomerAMGCreate(&(solversettings->precond));
        HYPRE_BoomerAMGSetMaxIter(solversettings->precond, 1);
        HYPRE_BoomerAMGSetTol(solversettings->precond, 0.0);
        HYPRE_BoomerAMGSetPrintLevel(solversettings->precond, 0);
        HYPRE_BoomerAMGSetCoarsenType(solversettings->precond, 6);
        HYPRE_BoomerAMGSetRelaxType(solversettings->precond, 6);
        HYPRE_BoomerAMGSetNumSweeps(solversettings->precond, 1);

        HYPRE_ParCSRPCGSetPrecond(solverdata->solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, solversettings->precond);
        HYPRE_ParCSRPCGSetup(solverdata->solver, solverdata->parA, solverdata->parb, solverdata->parx);

}

/*!
 * This is a driving function for solving next time step
 * of a given envical field.
 */
void environment_solve(systeminfo_t systeminfo,settings_t settings,grid_t *grid,environment_t **environment,solverdata_t *solverdata,solversettings_t *solversettings)
{
        int i, j, k;
        int idx;
        int var;
        double *values;
        int stepIter = 0;
        long long nvalues = grid->localsize.x * grid->localsize.y * grid->localsize.z;
        double *pc;
        if (systeminfo.rank == 0 ) {
                printf("Solving field.");
                fflush(stdout);
        }

        values = (double *) calloc(nvalues, sizeof(double));
        pc = (double *) calloc(nvalues, sizeof(double));

        while (stepIter < solversettings->numberofiters) {
                if (solversettings->enviter > 0) {
                        // update right hand side
                        for(var=0; var<settings.numberoffields; var++) {

                                environment_pcfunction(systeminfo,settings,var,grid,pc);

                                HYPRE_SStructVectorGetBoxValues(solverdata->x, 0, solversettings->lower, solversettings->upper, var, values);
                                HYPRE_SStructVectorSetBoxValues(solverdata->b, 0, solversettings->lower, solversettings->upper, var, values);

                                for (i = 0; i < nvalues; i++)
                                        values[i] = solversettings->z[var] * (*environment)[var].boundarycondition;

                                if (systeminfo.coords[systeminfo.rank][0] == 0) {
                                        environment_setboundary(0, -1, solversettings);
                                        HYPRE_SStructVectorAddToBoxValues(solverdata->b, 0, solversettings->bclower, solversettings->bcupper, var, values);
                                }
                                if (systeminfo.coords[systeminfo.rank][0] == systeminfo.dim[0] - 1) {
                                        environment_setboundary(0, 1, solversettings);
                                        HYPRE_SStructVectorAddToBoxValues(solverdata->b, 0, solversettings->bclower, solversettings->bcupper, var, values);
                                }
                                if (systeminfo.coords[systeminfo.rank][1] == 0) {
                                        environment_setboundary(1, -1, solversettings);
                                        HYPRE_SStructVectorAddToBoxValues(solverdata->b, 0, solversettings->bclower, solversettings->bcupper, var, values);
                                }
                                if (systeminfo.coords[systeminfo.rank][1] == systeminfo.dim[1] - 1) {
                                        environment_setboundary(1, 1, solversettings);
                                        HYPRE_SStructVectorAddToBoxValues(solverdata->b, 0, solversettings->bclower, solversettings->bcupper, var, values);
                                }
                                if (systeminfo.coords[systeminfo.rank][2] == 0) {
                                        environment_setboundary(2, -1, solversettings);
                                        HYPRE_SStructVectorAddToBoxValues(solverdata->b, 0, solversettings->bclower, solversettings->bcupper, var, values);
                                }
                                if (systeminfo.coords[systeminfo.rank][2] == systeminfo.dim[2] - 1) {
                                        environment_setboundary(2, 1, solversettings);
                                        HYPRE_SStructVectorAddToBoxValues(solverdata->b, 0, solversettings->bclower, solversettings->bcupper, var, values);
                                }
                                HYPRE_SStructVectorAddToBoxValues(solverdata->b, 0, solversettings->lower, solversettings->upper, var, pc);
                                HYPRE_SStructVectorAssemble(solverdata->b);
                                HYPRE_SStructVectorAssemble(solverdata->x);
                        }

                }

                HYPRE_ParCSRPCGSolve(solverdata->solver, solverdata->parA, solverdata->parb, solverdata->parx);

                (solversettings->enviter)++;
                stepIter++;
        }

        for(var=0; var<settings.numberoffields; var++) {
                HYPRE_SStructVectorGather(solverdata->x);
                HYPRE_SStructVectorGetBoxValues(solverdata->x, 0, solversettings->lower, solversettings->upper, var, values);
                idx = 0;
                for (k = 0; k < grid->localsize.z; k++)
                        for (j = 0; j < grid->localsize.y; j++)
                                for (i = 0; i < grid->localsize.x; i++, idx++) {
                                        //envField[nch][gridSize.y * gridSize.z * i + gridSize.z * j +
                                        //               k] = values[idx];
                                }

        }

        free(values);
        free(pc);
}

void environment_compute(systeminfo_t systeminfo,settings_t settings,grid_t *grid,environment_t **environment, solverdata_t *solverdata,solversettings_t *solversettings) {

        if(settings.numberoffields==0) return;

        environment_initboundary(systeminfo,settings,grid,environment,solverdata,solversettings);
        environment_initsolver(systeminfo,solverdata,solversettings);
        environment_solve(systeminfo,settings,grid,environment,solverdata,solversettings);
        return;
}
