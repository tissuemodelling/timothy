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

int step_compute(systeminfo_t systeminfo, settings_t settings, cellsinfo_t *cellsinfo, celltype_t* celltype, grid_t *grid, environment_t **environment,cellcommdata_t *cellcommdata,interpdata_t *interpdata,cellenvdata_t ***cellenvdata,solverdata_t *solverdata, solversettings_t *solversettings);
void step_clean(settings_t settings,cellenvdata_t ***cellenvdata);
