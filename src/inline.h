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

/*! \file inline.h
 *  \brief contains various Timothy inline functions
 */
static inline int octnodeintersection(int idx,uint3dv_t minLocCode,uint3dv_t maxLocCode,cellsinfo_t cellsinfo)
{
  if( maxLocCode.x<=cellsinfo.octree[idx].xcode || minLocCode.x>=cellsinfo.octree[idx].xlimit ) return 0;
  if( maxLocCode.y<=cellsinfo.octree[idx].ycode || minLocCode.y>=cellsinfo.octree[idx].ylimit ) return 0;
  if( maxLocCode.z<=cellsinfo.octree[idx].zcode || minLocCode.z>=cellsinfo.octree[idx].zlimit ) return 0;
  return 1;
}

/*!
 * This function returns the value of the SPH kernel function for given two cells p1 and p2.
 * Kernel function as in:
 * - J.J.Monaghan,"Smoothed Particle Hydrodynamics",Annu.Rev.Astron.Astrophys.,1992,30:543-74
 * - V.Springel,"Smoothed Particle Hydrodynamics in Astrophysics",arXiv:1109.2219v1
 */
static inline float sphkernel(int dimension,double dist,double h,double h2,double h3)
{
  double u,c=1.0;
  if(dist<0.0) return 0.0;
  u=dist/h;
  if(dimension==2) c=40/(7*M_PI*h2);
  if(dimension==3) c=8/(M_PI*h3);
  if(u>=0.0 && u<=0.5) return c*(1-6*u*u+6*u*u*u);
  if(u>0.5 && u<=1.0) return c*(2*(1-u)*(1-u)*(1-u));
  if(u>1.0) return 0.0;
  return 0.0;
}

/*!
 * This function return the value of the SPH kernel function gradient.
 * Kernel function as in
 * - J.J.Monaghan,"Smoothed Particle Hydrodynamics",Annu.Rev.Astron.Astrophys.,1992,30:543-74
 * - V.Springel,"Smoothed Particle Hydrodynamics in Astrophysics",arXiv:1109.2219v1
 */
static inline int sphgradient(int dimension, celldata_t c1, celldata_t c2, double3dv_t* f,double h,double h3,double h4)
{
  double u=1.0,c=1.0,w=1.0;
  double dist;

  f->x=0.0;
  f->y=0.0;
  f->z=0.0;
  dist = sqrt((c1.x - c2.x) * (c1.x - c2.x) + (c1.y - c2.y) * (c1.y - c2.y) + (c1.z - c2.z) * (c1.z - c2.z));

  if(dist>=0.0 && dist<=h) {
    u=dist/h;
    if(dimension==2) c=(40*8)/(7*M_PI*h3);
    if(dimension==3) c=48/(M_PI*h4);

    if(u>=0.0 && u<=0.5) w=-2.0*u+3.0*u*u;
    if(u>0.5 && u<=1.0) w=-(1.0-u)*(1.0-u);
    if(u>1.0) w=0.0;

    f->x=w*c*(c2.x-c1.x)/dist;
    f->y=w*c*(c2.y-c1.y)/dist;
    if(dimension==3) f->z=w*c*(c2.z-c1.z)/dist;
  }
  return 0;
}
