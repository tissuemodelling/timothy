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

#include "fieldgradient.h"
#include "utils.h"
#include "cells.h"


void fieldgradient_allocate(systeminfo_t systeminfo,fieldgradientdata_t *fieldgradientdata,grid_t *grid)
{

        MPI_Cart_shift(systeminfo.MPI_CART_COMM,0,1,&fieldgradientdata->bx0,&fieldgradientdata->bx1);
        MPI_Cart_shift(systeminfo.MPI_CART_COMM,1,1,&fieldgradientdata->by0,&fieldgradientdata->by1);
        MPI_Cart_shift(systeminfo.MPI_CART_COMM,2,1,&fieldgradientdata->bz0,&fieldgradientdata->bz1);

        /* allocate send buffers */
        if(fieldgradientdata->bx0!=MPI_PROC_NULL) fieldgradientdata->sendx0=(double*)calloc(grid->localsize.y*grid->localsize.z,sizeof(double));
        if(fieldgradientdata->bx1!=MPI_PROC_NULL) fieldgradientdata->sendx1=(double*)calloc(grid->localsize.y*grid->localsize.z,sizeof(double));
        if(fieldgradientdata->by0!=MPI_PROC_NULL) fieldgradientdata->sendy0=(double*)calloc(grid->localsize.x*grid->localsize.z,sizeof(double));
        if(fieldgradientdata->by1!=MPI_PROC_NULL) fieldgradientdata->sendy1=(double*)calloc(grid->localsize.x*grid->localsize.z,sizeof(double));
        if(fieldgradientdata->bz0!=MPI_PROC_NULL) fieldgradientdata->sendz0=(double*)calloc(grid->localsize.y*grid->localsize.x,sizeof(double));
        if(fieldgradientdata->bz1!=MPI_PROC_NULL) fieldgradientdata->sendz1=(double*)calloc(grid->localsize.y*grid->localsize.x,sizeof(double));

        /* allocate receive buffers */
        if(fieldgradientdata->bx0!=MPI_PROC_NULL) fieldgradientdata->recvx0=(double*)calloc(grid->localsize.y*grid->localsize.z,sizeof(double));
        if(fieldgradientdata->bx1!=MPI_PROC_NULL) fieldgradientdata->recvx1=(double*)calloc(grid->localsize.y*grid->localsize.z,sizeof(double));
        if(fieldgradientdata->by0!=MPI_PROC_NULL) fieldgradientdata->recvy0=(double*)calloc(grid->localsize.x*grid->localsize.z,sizeof(double));
        if(fieldgradientdata->by1!=MPI_PROC_NULL) fieldgradientdata->recvy1=(double*)calloc(grid->localsize.x*grid->localsize.z,sizeof(double));
        if(fieldgradientdata->bz0!=MPI_PROC_NULL) fieldgradientdata->recvz0=(double*)calloc(grid->localsize.y*grid->localsize.x,sizeof(double));
        if(fieldgradientdata->bz1!=MPI_PROC_NULL) fieldgradientdata->recvz1=(double*)calloc(grid->localsize.y*grid->localsize.x,sizeof(double));

        return;
}

void fieldgradient_exchangeinit(systeminfo_t systeminfo, fieldgradientdata_t *fieldgradientdata, environment_t **environment, grid_t *grid, int f)
{
        int i,j,k;

        for(i=0; i<grid->localsize.x; i++)
                for(j=0; j<grid->localsize.y; j++)
                        for(k=0; k<grid->localsize.z; k++) {
                                double val;
                                val = (*environment)[f].data[grid->localsize.z*grid->localsize.y*i+grid->localsize.z*j+k];
                                if(i==0 && fieldgradientdata->bx0!=MPI_PROC_NULL) fieldgradientdata->sendx0[grid->localsize.z*j+k]=val;
                                if(i==grid->localsize.x-1 && fieldgradientdata->bx1!=MPI_PROC_NULL) fieldgradientdata->sendx1[grid->localsize.z*j+k]=val;
                                if(j==0 && fieldgradientdata->by0!=MPI_PROC_NULL) fieldgradientdata->sendy0[grid->localsize.z*i+k]=val;
                                if(j==grid->localsize.y-1 && fieldgradientdata->by1!=MPI_PROC_NULL) fieldgradientdata->sendy1[grid->localsize.z*i+k]=val;
                                if(k==0 && fieldgradientdata->bz0!=MPI_PROC_NULL) fieldgradientdata->sendz0[grid->localsize.y*i+j]=val;
                                if(k==grid->localsize.z-1 && fieldgradientdata->bz1!=MPI_PROC_NULL) fieldgradientdata->sendz1[grid->localsize.y*i+j]=val;
                        }

        if(fieldgradientdata->bx0!=MPI_PROC_NULL) {
                MPI_Isend(fieldgradientdata->sendx0,grid->localsize.y*grid->localsize.z,MPI_DOUBLE,fieldgradientdata->bx0,systeminfo.rank,systeminfo.MPI_CART_COMM,&(fieldgradientdata->reqsend)[0]);
                MPI_Irecv(fieldgradientdata->recvx0,grid->localsize.y*grid->localsize.z,MPI_DOUBLE,fieldgradientdata->bx0,systeminfo.size+fieldgradientdata->bx0,systeminfo.MPI_CART_COMM,&(fieldgradientdata->reqrecv)[0]);
        }
        if(fieldgradientdata->bx1!=MPI_PROC_NULL) {
                MPI_Isend(fieldgradientdata->sendx1,grid->localsize.y*grid->localsize.z,MPI_DOUBLE,fieldgradientdata->bx1,systeminfo.size+systeminfo.rank,systeminfo.MPI_CART_COMM,&(fieldgradientdata->reqsend)[1]);
                MPI_Irecv(fieldgradientdata->recvx1,grid->localsize.y*grid->localsize.z,MPI_DOUBLE,fieldgradientdata->bx1,fieldgradientdata->bx1,systeminfo.MPI_CART_COMM,&(fieldgradientdata->reqrecv)[1]);
        }
        if(fieldgradientdata->by0!=MPI_PROC_NULL) {
                MPI_Isend(fieldgradientdata->sendy0,grid->localsize.x*grid->localsize.z,MPI_DOUBLE,fieldgradientdata->by0,2*systeminfo.size+systeminfo.rank,systeminfo.MPI_CART_COMM,&(fieldgradientdata->reqsend)[2]);
                MPI_Irecv(fieldgradientdata->recvy0,grid->localsize.x*grid->localsize.z,MPI_DOUBLE,fieldgradientdata->by0,3*systeminfo.size+fieldgradientdata->by0,systeminfo.MPI_CART_COMM,&(fieldgradientdata->reqrecv)[2]);
        }
        if(fieldgradientdata->by1!=MPI_PROC_NULL) {
                MPI_Isend(fieldgradientdata->sendy1,grid->localsize.x*grid->localsize.z,MPI_DOUBLE,fieldgradientdata->by1,3*systeminfo.size+systeminfo.rank,systeminfo.MPI_CART_COMM,&(fieldgradientdata->reqsend)[3]);
                MPI_Irecv(fieldgradientdata->recvy1,grid->localsize.x*grid->localsize.z,MPI_DOUBLE,fieldgradientdata->by1,2*systeminfo.size+fieldgradientdata->by1,systeminfo.MPI_CART_COMM,&(fieldgradientdata->reqrecv)[3]);
        }
        if(fieldgradientdata->bz0!=MPI_PROC_NULL) {
                MPI_Isend(fieldgradientdata->sendz0,grid->localsize.y*grid->localsize.x,MPI_DOUBLE,fieldgradientdata->bz0,4*systeminfo.size+systeminfo.rank,systeminfo.MPI_CART_COMM,&(fieldgradientdata->reqsend)[4]);
                MPI_Irecv(fieldgradientdata->recvz0,grid->localsize.y*grid->localsize.x,MPI_DOUBLE,fieldgradientdata->bz0,5*systeminfo.size+fieldgradientdata->bz0,systeminfo.MPI_CART_COMM,&(fieldgradientdata->reqrecv)[4]);
        }
        if(fieldgradientdata->bz1!=MPI_PROC_NULL) {
                MPI_Isend(fieldgradientdata->sendz1,grid->localsize.y*grid->localsize.x,MPI_DOUBLE,fieldgradientdata->bz1,5*systeminfo.size+systeminfo.rank,systeminfo.MPI_CART_COMM,&(fieldgradientdata->reqsend)[5]);
                MPI_Irecv(fieldgradientdata->recvz1,grid->localsize.y*grid->localsize.x,MPI_DOUBLE,fieldgradientdata->bz1,4*systeminfo.size+fieldgradientdata->bz1,systeminfo.MPI_CART_COMM,&(fieldgradientdata->reqrecv)[5]);
        }

        return;
}

void fieldgradient_exchangewait(systeminfo_t systeminfo, fieldgradientdata_t *fieldgradientdata)
{
        MPI_Status status;
        if(fieldgradientdata->bx0!=MPI_PROC_NULL) {
                if (MPI_Wait(&fieldgradientdata->reqsend[0], &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication fieldgradientdata->reqsend[0]", __FILE__, __LINE__);
                if (MPI_Wait(&fieldgradientdata->reqrecv[0], &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication fieldgradientdata->reqrecv[0]", __FILE__, __LINE__);
        }
        if(fieldgradientdata->bx1!=MPI_PROC_NULL) {
                if (MPI_Wait(&fieldgradientdata->reqsend[1], &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication fieldgradientdata->reqsend[1]", __FILE__, __LINE__);
                if (MPI_Wait(&fieldgradientdata->reqrecv[1], &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication fieldgradientdata->reqrecv[1]", __FILE__, __LINE__);
        }
        if(fieldgradientdata->by0!=MPI_PROC_NULL) {
                if (MPI_Wait(&fieldgradientdata->reqsend[2], &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication fieldgradientdata->reqsend[2]", __FILE__, __LINE__);
                if (MPI_Wait(&fieldgradientdata->reqrecv[2], &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication fieldgradientdata->reqrecv[2]", __FILE__, __LINE__);
        }
        if(fieldgradientdata->by1!=MPI_PROC_NULL) {
                if (MPI_Wait(&fieldgradientdata->reqsend[3], &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication fieldgradientdata->reqsend[3]", __FILE__, __LINE__);
                if (MPI_Wait(&fieldgradientdata->reqrecv[3], &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication fieldgradientdata->reqrecv[3]", __FILE__, __LINE__);
        }
        if(fieldgradientdata->bz0!=MPI_PROC_NULL) {
                if (MPI_Wait(&fieldgradientdata->reqsend[4], &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication fieldgradientdata->reqsend[4]", __FILE__, __LINE__);
                if (MPI_Wait(&fieldgradientdata->reqrecv[4], &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication fieldgradientdata->reqrecv[4]", __FILE__, __LINE__);
        }
        if(fieldgradientdata->bz1!=MPI_PROC_NULL) {
                if (MPI_Wait(&fieldgradientdata->reqsend[5], &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication fieldgradientdata->reqsend[5]", __FILE__, __LINE__);
                if (MPI_Wait(&fieldgradientdata->reqrecv[5], &status) != MPI_SUCCESS)
                        terminate(systeminfo,"communication fieldgradientdata->reqrecv[5]", __FILE__, __LINE__);
        }

        return;
}


void fieldgradient_compute(systeminfo_t systeminfo, fieldgradientdata_t *fieldgradientdata, grid_t *grid, environment_t **environment, int f)
{
        int i,j,k;
        /* compute internal data */
        for(i=0; i<grid->localsize.x; i++)
                for(j=0; j<grid->localsize.y; j++)
                        for(k=0; k<grid->localsize.z; k++) {
                                /* x coord gradient */
                                if(i!=0 && i!=grid->localsize.x-1) {
                                        (*environment)[f].gradient[3*grid->localsize.z*grid->localsize.y*i+3*grid->localsize.z*j+3*k+0]=
                                                ((*environment)[f].data[grid->localsize.z*grid->localsize.y*(i+1)+grid->localsize.z*j+k]
                                                 -(*environment)[f].data[grid->localsize.z*grid->localsize.y*(i-1)+grid->localsize.z*j+k])/grid->resolution;

                                }
                                /* y coord gradient */
                                if(j!=0 && j!=grid->localsize.y-1) {
                                        (*environment)[f].gradient[3*grid->localsize.z*grid->localsize.y*i+3*grid->localsize.z*j+3*k+1]=
                                                ((*environment)[f].data[grid->localsize.z*grid->localsize.y*i+grid->localsize.z*(j+1)+k]
                                                 -(*environment)[f].data[grid->localsize.z*grid->localsize.y*i+grid->localsize.z*(j-1)+k])/grid->resolution;
                                }
                                /* z coord gradient */
                                if(k!=0 && k!=grid->localsize.z-1) {
                                        (*environment)[f].gradient[3*grid->localsize.z*grid->localsize.y*i+3*grid->localsize.z*j+3*k+2]=
                                                ((*environment)[f].data[grid->localsize.z*grid->localsize.y*i+grid->localsize.z*j+k+1]
                                                 -(*environment)[f].data[grid->localsize.z*grid->localsize.y*i+grid->localsize.z*j+k-1])/grid->resolution;

                                }
                        }
        /* wait for boundary data to arrive */
        fieldgradient_exchangewait(systeminfo,fieldgradientdata);
        /* update with boundary data */
        for(j=0; j<grid->localsize.y; j++)
                for(k=0; k<grid->localsize.z; k++) {
                        double x0,x1;
                        if(fieldgradientdata->bx0!=MPI_PROC_NULL) x0=fieldgradientdata->recvx0[grid->localsize.z*j+k];
                        else x0=(*environment)[f].data[grid->localsize.z*grid->localsize.y*0+grid->localsize.z*j+k];
                        if(fieldgradientdata->bx1!=MPI_PROC_NULL) x1=fieldgradientdata->recvx1[grid->localsize.z*j+k];
                        else x1=(*environment)[f].data[grid->localsize.z*grid->localsize.y*(grid->localsize.x-1)+grid->localsize.z*j+k];
                        (*environment)[f].gradient[3*grid->localsize.z*grid->localsize.y*0+3*grid->localsize.z*j+3*k+0]=
                                ((*environment)[f].data[grid->localsize.z*grid->localsize.y*1+grid->localsize.z*j+k]
                                 - x0)/grid->resolution;
                        (*environment)[f].gradient[3*grid->localsize.z*grid->localsize.y*(grid->localsize.x-1)+3*grid->localsize.z*j+3*k+0]=
                                (x1
                                 - (*environment)[f].data[grid->localsize.z*grid->localsize.y*(grid->localsize.x-2)+grid->localsize.z*j+k])/grid->resolution;
                }
        for(i=0; i<grid->localsize.x; i++)
                for(k=0; k<grid->localsize.z; k++) {
                        double y0,y1;
                        if(fieldgradientdata->by0!=MPI_PROC_NULL) y0=fieldgradientdata->recvy0[grid->localsize.z*i+k];
                        else y0=(*environment)[f].data[grid->localsize.z*grid->localsize.y*i+grid->localsize.z*0+k];
                        if(fieldgradientdata->by1!=MPI_PROC_NULL) y1=fieldgradientdata->recvy1[grid->localsize.z*i+k];
                        else y1=(*environment)[f].data[grid->localsize.z*grid->localsize.y*i+grid->localsize.z*(grid->localsize.y-1)+k];
                        (*environment)[f].gradient[3*grid->localsize.z*grid->localsize.y*i+3*grid->localsize.z*0+3*k+1]=
                                ((*environment)[f].data[grid->localsize.z*grid->localsize.y*i+grid->localsize.z*1+k]
                                 - y0)/grid->resolution;
                        (*environment)[f].gradient[3*grid->localsize.z*grid->localsize.y*i+3*grid->localsize.z*(grid->localsize.y-1)+3*k+1]=
                                (y1
                                 - (*environment)[f].data[grid->localsize.z*grid->localsize.y*i+grid->localsize.z*(grid->localsize.y-2)+k])/grid->resolution;
                }
        for(i=0; i<grid->localsize.x; i++)
                for(j=0; j<grid->localsize.y; j++) {
                        double z0,z1;
                        if(fieldgradientdata->bz0!=MPI_PROC_NULL) z0=fieldgradientdata->recvz0[grid->localsize.y*i+j];
                        else z0=(*environment)[f].data[grid->localsize.z*grid->localsize.y*i+grid->localsize.z*j+0];
                        if(fieldgradientdata->bz1!=MPI_PROC_NULL) z1=fieldgradientdata->recvz1[grid->localsize.y*i+j];
                        else z1=(*environment)[f].data[grid->localsize.z*grid->localsize.y*i+grid->localsize.z*j+(grid->localsize.z-1)];
                        (*environment)[f].gradient[3*grid->localsize.z*grid->localsize.y*i+3*grid->localsize.z*j+3*0+2]=
                                ((*environment)[f].data[grid->localsize.z*grid->localsize.y*i+grid->localsize.z*j+1]
                                 - z0)/grid->resolution;
                        (*environment)[f].gradient[3*grid->localsize.z*grid->localsize.y*i+3*grid->localsize.z*j+3*(grid->localsize.z-1)+2]=
                                (z1
                                 - (*environment)[f].data[grid->localsize.z*grid->localsize.y*i+grid->localsize.z*j+(grid->localsize.z-2)])/grid->resolution;
                }

        return;
}

void fieldgradient_free(fieldgradientdata_t *fieldgradientdata) {

        if(fieldgradientdata->bx0!=MPI_PROC_NULL) {
                free(fieldgradientdata->sendx0);
                free(fieldgradientdata->recvx0);
        }
        if(fieldgradientdata->bx1!=MPI_PROC_NULL) {
                free(fieldgradientdata->sendx1);
                free(fieldgradientdata->recvx1);
        }
        if(fieldgradientdata->by0!=MPI_PROC_NULL) {
                free(fieldgradientdata->sendy0);
                free(fieldgradientdata->recvy0);
        }
        if(fieldgradientdata->by1!=MPI_PROC_NULL) {
                free(fieldgradientdata->sendy1);
                free(fieldgradientdata->recvy1);
        }
        if(fieldgradientdata->bz0!=MPI_PROC_NULL) {
                free(fieldgradientdata->sendz0);
                free(fieldgradientdata->recvz0);
        }
        if(fieldgradientdata->bz1!=MPI_PROC_NULL) {
                free(fieldgradientdata->sendz1);
                free(fieldgradientdata->recvz1);
        }

        return;
}

void fieldgradient(systeminfo_t systeminfo, settings_t settings, environment_t **environment, grid_t *grid)
{
        int f;
        fieldgradientdata_t fieldgradientdata;
        if (settings.numberoffields==0) return;
        fieldgradient_allocate(systeminfo,&fieldgradientdata,grid);
        for(f=0; f<settings.numberoffields; f++) {
                fieldgradient_exchangeinit(systeminfo,&fieldgradientdata,environment,grid,f);
                fieldgradient_compute(systeminfo,&fieldgradientdata,grid,environment,f);
        }
        /*for(chf=0; chf<NCHEM; chf++) {
                if(chf==OXYG-NGLOB && !oxygen) continue;
                if(chf==GLUC-NGLOB && !glucose) continue;
                if(chf==HYDR-NGLOB && !hydrogenIon) continue;
                initFieldHaloExchange(systeminfo,&fieldgradientdata,chf);
                //waitFieldHaloExchange();
                computeFieldGradient(&fieldgradientdata,chf);

           }*/
        fieldgradient_free(&fieldgradientdata);
        return;
}
