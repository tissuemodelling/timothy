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
 
void createexportlist(systeminfo_t systeminfo,settings_t settings,cellsinfo_t cellsinfo,grid_t grid,struct Zoltan_Struct *ztn,celltype_t* celltype,cellcommdata_t *cellcommdata,fieldcommdata_t *fieldcommdata);
void exchangecleanup(systeminfo_t systeminfo,cellsinfo_t cellsinfo,cellcommdata_t *cellcommdata,fieldcommdata_t *fieldcommdata);
void cellssendrecv(systeminfo_t systeminfo, cellsinfo_t cellsinfo, cellcommdata_t *cellcommdata);
void cellswait(systeminfo_t systeminfo, cellsinfo_t cellsinfo, cellcommdata_t *cellcommdata);
void datasendrecv(systeminfo_t systeminfo, cellsinfo_t cellsinfo, cellcommdata_t *cellcommdata);
void datawait(systeminfo_t systeminfo, cellsinfo_t cellsinfo, cellcommdata_t *cellcommdata);
