/* **************************************************************************
 * This file is part of Timothy
 *
 * Copyright (c) 2014/15 Maciej Cytowski
 * Copyright (c) 2014/15 ICM, University of Warsaw, Poland
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

void patches_alloc(systeminfo_t systeminfo, settings_t settings, patches_t *patches, cellsinfo_t *cellsinfo, grid_t *grid);
void patches_free(patches_t *patches);
void patches_env2cellsinit(systeminfo_t systeminfo, settings_t settings, patches_t *patches, grid_t *grid, environment_t **environment);
void patches_env2cellswait(systeminfo_t systeminfo, settings_t settings, patches_t *patches,cellsinfo_t *cellsinfo, grid_t *grid,cellenvdata_t ***cellenvdata);
void patches_cells2envwait(systeminfo_t systeminfo, settings_t settings, patches_t *patches, grid_t *grid, environment_t **environment);
void patches_cells2envinit(systeminfo_t systeminfo, settings_t settings, patches_t *patches, cellsinfo_t *cellsinfo, grid_t *grid);
