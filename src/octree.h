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

 /* properties of the affine transformation */
 double3dv_t affShift;
 double affScale;

void octbuild(systeminfo_t systeminfo,cellsinfo_t *cellsinfo,celltype_t* celltype);
void octfree(cellsinfo_t *cellsinfo);
void octheapinit(systeminfo_t systeminfo,octheap_t *ttheap);
void octcomputebox(int64_t c,uint3dv_t *minLocCode,uint3dv_t *maxLocCode,cellsinfo_t cellsinfo,celltype_t* celltype);
int octlocateregion(uint3dv_t minLocCode,uint3dv_t maxLocCode,cellsinfo_t cellsinfo);
void octheappush(systeminfo_t systeminfo,octheap_t *ttheap,int idx);
int octheappop(octheap_t *ttheap);
void octheapfree(octheap_t *ttheap);
static inline int octnodeintersection(int idx,uint3dv_t minLocCode,uint3dv_t maxLocCode,cellsinfo_t cellsinfo);
void octcomputeboxr(int64_t c,uint3dv_t *minLocCode,uint3dv_t *maxLocCode,cellcommdata_t cellcommdata,celltype_t* celltype);
