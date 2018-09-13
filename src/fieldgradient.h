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

 typedef struct fieldgradientdata_t {
   int bx0,bx1;
   int by0,by1;
   int bz0,bz1;
   double *sendx0;
   double *sendx1;
   double *sendy0;
   double *sendy1;
   double *sendz0;
   double *sendz1;
   double *recvx0;
   double *recvx1;
   double *recvy0;
   double *recvy1;
   double *recvz0;
   double *recvz1;
   MPI_Request reqsend[6];
   MPI_Request reqrecv[6];
 } fieldgradientdata_t;

void fieldsInit();
void fieldsSolve(settings_t settings,cellsinfo_t *cellsinfo);
void fieldgradient_allocate(systeminfo_t systeminfo,fieldgradientdata_t *fieldgradientdata,grid_t *grid);
void fieldgradient_free(fieldgradientdata_t *fieldgradientdata);
void fieldgradient_exchangeinit(systeminfo_t systeminfo,fieldgradientdata_t *fieldgradientdata,environment_t **environment,grid_t *grid,int f);
void fieldgradient_compute(systeminfo_t systeminfo,fieldgradientdata_t *fieldgradientdata,grid_t *grid,environment_t **environment,int f);
void fieldgradient(systeminfo_t systeminfo, settings_t settings, environment_t **environment,grid_t *grid);
