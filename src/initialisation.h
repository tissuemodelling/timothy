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

void getsysteminfo(systeminfo_t* systeminfo);
void initialsettings(settings_t* settings);
void initialcelltype(int numberofcelltypes,int numberoffields,celltype_t* celltype);
void initialfields(int numberoffields,environment_t* environment);
void initialisation(int argc, char **argv, systeminfo_t *systeminfo, settings_t* settings,celltype_t** celltype,environment_t** environment);
void initcount(cellcount_t *cellcount);
void allocatecells(systeminfo_t systeminfo,settings_t settings,celltype_t *celltype,cellsinfo_t *cellsinfo);
void printinfo(systeminfo_t systeminfo);
