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

void lbinit(int argc, char **argv, MPI_Comm Comm,systeminfo_t systeminfo,struct Zoltan_Struct **ztn,cellsinfo_t *cellsinfo);
void lbexchange(systeminfo_t systeminfo,struct Zoltan_Struct *ztn);
void lbdestroy(struct Zoltan_Struct **ztn);
