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

int checkendiannes(systeminfo_t *systeminfo);
void swapnbyte(char *data, int n, int m);
void updateglobalcounts(cellsinfo_t* cellsinfo);
void terminate(systeminfo_t systeminfo, char *msg, char *file, int line);
void stopRun(int ierr, char *name, char *file, int line);
size_t getMemoryPerProcess(int32_t lsize);
void getLocalRankAndSize(int rank, int size, int32_t * lrank, int32_t * lsize);
void randomstreaminit(systeminfo_t *systeminfo,settings_t *settings);
void printstatistics(systeminfo_t systeminfo,settings_t settings,cellsinfo_t cellsinfo,statistics_t* statistics);
