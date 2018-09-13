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
#include <string.h>
#include <mpi.h>
#include <inttypes.h>
#include <float.h>
#include <math.h>

#define _GNU_SOURCE
#include <fcntl.h>
#include <sys/stat.h>

#include "global.h"

#include "io.h"
#include "utils.h"
#include "cells.h"

/*! \file io.c
 *  \brief contains I/O functions
 */

#define INPUT_FILE_LINE_LENGTH 256

void readRstFile(int argc, char **argv);

void readcellpositions(systeminfo_t systeminfo,settings_t settings,celltype_t *celltype,cellsinfo_t *cellsinfo) {
        FILE *fhandle;
        char buf[128],bufx[16],bufy[16],bufz[16],bufn[16];
        int i;
        /* cell coordinate files are read only by rank 0 */
        if(systeminfo.rank==0) {
                for(i=0; i<settings.numberofcelltypes; i++) {
                        if(strlen(celltype[i].inputfile)==0) {
                                printf("filename for %s cells coordinates not specified\n",celltype[i].name);
                                continue;
                        } else {
                                fhandle = fopen(celltype[i].inputfile, "r");
                                if(fhandle==NULL) {
                                        printf("filename for %s cells coordinates not found\n",celltype[i].name);
                                        continue;
                                }
                                while(!feof(fhandle)) {
                                        char *ret;
                                        int ncoords;
                                        if(cellsinfo->localcount.n==settings.maxlocalcells) break;
                                        ret = fgets(buf, INPUT_FILE_LINE_LENGTH, fhandle);
                                        if (feof(fhandle))
                                                break;
                                        ncoords=sscanf(buf, "%s%s%s%s", bufx, bufy, bufz, bufn);

                                        if(ncoords==cellsinfo->dimension) {
                                                cellsinfo->cells[cellsinfo->localcount.n].x=atof(bufx);
                                                cellsinfo->cells[cellsinfo->localcount.n].y=atof(bufy);
                                                if(ncoords>2)
                                                        cellsinfo->cells[cellsinfo->localcount.n].z=atof(bufz);
                                                cellsinfo->cells[cellsinfo->localcount.n].ctype=i;
                                                cellsinfo->cells[cellsinfo->localcount.n].size=celltype[i].size;
                                                cellsinfo->cells[cellsinfo->localcount.n].density=0.0;
                                                cellsinfo->localcount.n+=1;
                                                cellsinfo->localcount.g0phase+=1;
                                                cellsinfo->localtypecount[i].n+=1;
                                                cellsinfo->localtypecount[i].g0phase+=1;
                                        }

                                        if(ncoords==cellsinfo->dimension+1) {
                                                int nrandom;
                                                cellsinfo->cells[cellsinfo->localcount.n].x=atof(bufx);
                                                cellsinfo->cells[cellsinfo->localcount.n].y=atof(bufy);
                                                if(cellsinfo->dimension==2) cellsinfo->cells[cellsinfo->localcount.n].z=0.0;
                                                if(cellsinfo->dimension==3) cellsinfo->cells[cellsinfo->localcount.n].z=atof(bufz);
                                                cellsinfo->cells[cellsinfo->localcount.n].ctype=i;
                                                cellsinfo->cells[cellsinfo->localcount.n].size=celltype[i].size;
                                                cellsinfo->cells[cellsinfo->localcount.n].density=0.0;
                                                cellsinfo->localcount.n+=1;
                                                cellsinfo->localcount.g0phase+=1;
                                                cellsinfo->localtypecount[i].n+=1;
                                                cellsinfo->localtypecount[i].g0phase+=1;
                                                if(cellsinfo->dimension==2) nrandom=atoi(bufz);
                                                if(cellsinfo->dimension==3) nrandom=atoi(bufn);
                                                cellsrandominit(nrandom,i,systeminfo,settings,celltype,cellsinfo);
                                        }

                                        if(ncoords<cellsinfo->dimension || ncoords>cellsinfo->dimension+1)
                                                printf("warning: missing cell coordinates in %s\n",celltype[i].inputfile);

                                }
                                fclose(fhandle);
                                if(cellsinfo->localcount.n==settings.maxlocalcells) {
                                        printf("warning: too many local cells, skipping rest of file %s\n",celltype[i].inputfile);
                                        break;
                                }
                        }
                }
        }
        return;
}

void readenvfile(systeminfo_t systeminfo,settings_t* settings,environment_t* environment) {
        char envfile[FNLEN];
        char buf[400], buf1[100], buf2[100], buf3[100];
        FILE *fhandle;
        int i,j;
        int nr;
        int firstnameset=0;

        if(settings->numberoffields==0) return;

        sprintf(envfile,"environment.inp");

        if (systeminfo.rank == 0) {
                printf("reading environment file: %s\n", envfile);
                fflush(stdout);
        }

        fhandle = fopen(envfile, "r");
        if (fhandle == NULL)
                terminate(systeminfo, "could not open environment file", __FILE__, __LINE__);

        for(i=0; i<settings->numberoffields; i++) {

                while(!feof(fhandle)) {
                        char *ret;
                        ret = fgets(buf, INPUT_FILE_LINE_LENGTH, fhandle);

                        if (sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
                                continue;

                        if (feof(fhandle))
                                break;

                        if (buf1[0] == '#')
                                continue;

                        if (strcmp(buf1,"ENVNAME")==0) {
                                if(firstnameset) i++;
                                else firstnameset=1;
                                strcpy(environment[i].name, buf2);
                                continue;
                        }

                        if (strcmp(buf1,"DC") == 0) { environment[i].diffusioncoefficient=atof(buf2); continue; }
                        if (strcmp(buf1,"BC") == 0) { environment[i].boundarycondition=atof(buf2); continue; }
                        if (strcmp(buf1,"ICMEAN") == 0) { environment[i].initialconditionmean=atof(buf2); continue; }
                        if (strcmp(buf1,"ICVAR") == 0) { environment[i].initialconditionvariance=atof(buf2); continue; }
                        if (strcmp(buf1,"LAMBDA") == 0) { environment[i].lambdadelay=atof(buf2); continue; }

                }
        }
        fclose(fhandle);
        return;
}

void readcellsfile(systeminfo_t systeminfo, settings_t* settings, celltype_t* celltype) {
        char cellsfile[FNLEN];
        char buf[400], buf1[100], buf2[100], buf3[100], buf4[100],buf5[100];
        FILE *fhandle;
        int i,k;
        int nr;
        int nbytes,bytesread;
        int firstnameset=0;

        sprintf(cellsfile,"cells.inp");

        if (systeminfo.rank == 0)
                printf("reading cells file: %s\n", cellsfile);
        fflush(stdout);

        fhandle = fopen(cellsfile, "r");
        if (fhandle == NULL)
                terminate(systeminfo, "could not open cell type file", __FILE__, __LINE__);

        for(i=0; i<settings->numberofcelltypes; i++) {

                while(!feof(fhandle)) {
                        char *ret;
                        ret = fgets(buf, INPUT_FILE_LINE_LENGTH, fhandle);

                        if (sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
                                continue;

                        if (feof(fhandle))
                                break;

                        if (buf1[0] == '#')
                                continue;

                        if (strcmp(buf1,"CNAME")==0) {
                                if(firstnameset) i++;
                                else firstnameset=1;
                                strcpy(celltype[i].name, buf2);
                                continue;
                        }

                        if (strcmp(buf1,"G1") == 0) { celltype[i].g1=atof(buf2); continue; }
                        if (strcmp(buf1,"S") == 0) { celltype[i].s=atof(buf2); continue; }
                        if (strcmp(buf1,"G2") == 0) { celltype[i].g2=atof(buf2); continue; }
                        if (strcmp(buf1,"M") == 0) { celltype[i].m=atof(buf2); continue; }
                        if (strcmp(buf1,"V") == 0) { celltype[i].v=atof(buf2); continue; }
                        if (strcmp(buf1,"RD") == 0) { celltype[i].rd=atof(buf2); continue; }
                        if (strcmp(buf1,"SIZE") == 0 ) { celltype[i].size=atof(buf2); continue;}
                        if (strcmp(buf1,"H") == 0 ) { celltype[i].h=atof(buf2); continue; }
                        if (strcmp(buf1,"INPUT") == 0) { strcpy(celltype[i].inputfile, buf2); continue; }
                        if (strcmp(buf1,"CDENS") == 0) { celltype[i].criticaldensity=atof(buf2); continue; }
                        if (strcmp(buf1,"RANDOMMOVE") == 0) { celltype[i].randommove=atof(buf2); continue; }
                        if (strcmp(buf1,"ENVCONS") == 0) {
                                bytesread=0;
                                sscanf(buf, "%s%n",buf4,&nbytes);
                                bytesread+=nbytes;
                                for(k=0; k<settings->numberoffields; k++) {
                                        int ret;
                                        ret=sscanf(buf+bytesread,"%s%n",buf5,&nbytes);
                                        if(ret<=0) break;
                                        bytesread+=nbytes;
                                        celltype[i].consumption[k] = atof(buf5);
                                }
                        }
                        if (strcmp(buf1,"ENVPROD") == 0) {
                                bytesread=0;
                                sscanf(buf, "%s%n",buf4,&nbytes);
                                bytesread+=nbytes;
                                for(k=0; k<settings->numberoffields; k++) {
                                        int ret;
                                        ret=sscanf(buf+bytesread,"%s%n",buf5,&nbytes);
                                        if(ret<=0) break;
                                        bytesread+=nbytes;
                                        celltype[i].production[k] = atof(buf5);
                                }
                        }
                        if (strcmp(buf1,"ENVCL1") == 0) {
                                bytesread=0;
                                sscanf(buf, "%s%n",buf4,&nbytes);
                                bytesread+=nbytes;
                                for(k=0; k<settings->numberoffields; k++) {
                                        int ret;
                                        ret=sscanf(buf+bytesread,"%s%n",buf5,&nbytes);
                                        if(ret<=0) break;
                                        bytesread+=nbytes;
                                        celltype[i].criticallevel1[k] = atof(buf5);
                                }
                        }
                        if (strcmp(buf1,"ENVCL2") == 0) {
                                bytesread=0;
                                sscanf(buf, "%s%n",buf4,&nbytes);
                                bytesread+=nbytes;
                                for(k=0; k<settings->numberoffields; k++) {
                                        int ret;
                                        ret=sscanf(buf+bytesread,"%s%n",buf5,&nbytes);
                                        if(ret<=0) break;
                                        bytesread+=nbytes;
                                        celltype[i].criticallevel2[k] = atof(buf5);
                                }
                        }


                }
        }
        fclose(fhandle);
        return;
}


void readparamfile(int argc, char **argv, systeminfo_t systeminfo, settings_t* settings)
{
        #define NPAR 15
        char paramfile[FNLEN];
        char buf[400], buf1[100], buf2[100], buf3[100];
        FILE *fhandle;
        int i;

        strcpy(settings->outdir, "results");

        if (strlen(argv[1]) >= FNLEN)
                terminate(systeminfo, "parameter file name too long", __FILE__, __LINE__);

        sprintf(paramfile, "%s", argv[1]);

        if (systeminfo.rank == 0)
                printf("reading parameters file: %s\n", paramfile);

        fflush(stdout);

        fhandle = fopen(paramfile, "r");
        if (fhandle == NULL)
                terminate(systeminfo, "could not open parameter file", __FILE__, __LINE__);

        /* look for parameters in the file */
        while (!feof(fhandle)) {

                char *ret;
                ret = fgets(buf, INPUT_FILE_LINE_LENGTH, fhandle);

                if (sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
                        continue;

                if (feof(fhandle))
                        break;

                if (buf1[0] == '#')
                        continue;

                if (strcmp(buf1,"RESTART") == 0) { settings->restart=atoi(buf2); continue; }
                if (strcmp(buf1,"MAXCELLS") == 0) { settings->maxcells=atol(buf2); continue; }

        }

        /* read restart file if given */
        //if (settings->restart == 1)
        //readRstFile(argc, argv);

        /* rewind the file */
        rewind(fhandle);

        /* convert some of the parameters to data */
        while (!feof(fhandle)) {

                char *ret;
                ret = fgets(buf, INPUT_FILE_LINE_LENGTH, fhandle);

                if (sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
                        continue;

                if (feof(fhandle))
                        break;

                if (buf1[0] == '#')
                        continue;

                if (strcmp(buf1,"GFDT") == 0) { settings->gfdt=atof(buf2); continue; }
                if (strcmp(buf1,"GFH") == 0) { settings->gfh=atof(buf2); continue; }
                if (strcmp(buf1,"MAXCELLS") == 0) { settings->maxcells=atol(buf2); continue; }
                if (strcmp(buf1,"NSTEPS") == 0) { settings->numberofsteps=atoi(buf2); continue; }
                if (strcmp(buf1,"SECPERSTEP") == 0) { settings->secondsperstep=atof(buf2); continue; }
                if (strcmp(buf1,"NCELLTYPES") == 0) { settings->numberofcelltypes=atoi(buf2); continue; }
                if (strcmp(buf1,"NFIELDS") == 0 ) { settings->numberoffields=atoi(buf2); continue; }
                if (strcmp(buf1,"DIMENSIONS") == 0 ) { settings->dimension=atoi(buf2); continue; }
                if (strcmp(buf1,"RESTART") == 0) { settings->restart=atoi(buf2); continue; }
                if (strcmp(buf1,"RSTFILE") == 0) { strcpy(settings->rstfilename,buf2); continue; }
                if (strcmp(buf1,"OUTDIR") == 0) { strcpy(settings->outdir,buf2); continue; }
                if (strcmp(buf1,"VISOUTSTEP") == 0 ) { settings->visoutstep=atoi(buf2); continue; }
                if (strcmp(buf1,"STATOUTSTEP") == 0 ) { settings->statoutstep=atoi(buf2); continue; }
                if (strcmp(buf1,"RSTOUTSTEP") == 0 ) { settings->rstoutstep=atoi(buf2); continue; }
                if (strcmp(buf1,"MAXSPEED") == 0) { settings->maxspeed=atof(buf2); continue; }

        }

        fclose(fhandle);
        return;
}


void writevtk(systeminfo_t systeminfo,settings_t settings,cellsinfo_t cellsinfo) {

        int i,j;
        MPI_File fhandle;
        float *floatvectorfield;
        float *floatscalarfield;
        int *intscalarfield;
        int64_t nprev=0;
        char fstname[256];
        char buffer[256];
        MPI_Offset offset,goffset;

        floatvectorfield = (float *) calloc(cellsinfo.localcount.n * 3, sizeof(float));
        floatscalarfield = (float *) calloc(cellsinfo.localcount.n, sizeof(float));
        intscalarfield = (int *) calloc(cellsinfo.localcount.n, sizeof(int));

        sprintf(fstname,"vtk/step%08d.vtk",settings.step);

        goffset = 0;
        MPI_File_open(MPI_COMM_WORLD, fstname, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fhandle);

        /* truncate the file */
        MPI_File_set_size(fhandle, 0);

        /* write vtk header */
        sprintf(buffer,"# vtk DataFile Version 2.0\nTimothy output\nBINARY\nDATASET UNSTRUCTURED_GRID");
        if (systeminfo.rank == 0)
                MPI_File_write(fhandle, &buffer, strlen(buffer), MPI_BYTE, MPI_STATUS_IGNORE);
        goffset += strlen(buffer);
        MPI_File_seek(fhandle, goffset, MPI_SEEK_SET);

        for (i = 0; i < systeminfo.rank; i++)
                nprev += cellsinfo.cellsperproc[i];

        /* write cell positions */
        memset(buffer, 0, 256);
        sprintf(buffer, "\nPOINTS %" PRId64 " float\n",cellsinfo.globalcount.n);
        if (systeminfo.rank == 0) MPI_File_write(fhandle, &buffer, strlen(buffer), MPI_BYTE, MPI_STATUS_IGNORE);
        goffset += strlen(buffer);
        MPI_File_seek(fhandle, goffset, MPI_SEEK_SET);
        offset = nprev * sizeof(float) * 3;
        MPI_File_seek(fhandle, goffset+offset, MPI_SEEK_SET);
        for (i = 0; i < cellsinfo.localcount.n; i++) {
                floatvectorfield[3 * i] = (float) (cellsinfo.cells[i].x);
                floatvectorfield[3 * i + 1] = (float) (cellsinfo.cells[i].y);
                floatvectorfield[3 * i + 2] = (float) (cellsinfo.cells[i].z);
        }
        if (systeminfo.endian) swapnbyte((char *) floatvectorfield, cellsinfo.localcount.n * 3, sizeof(float));
        MPI_File_write(fhandle, floatvectorfield, cellsinfo.localcount.n * 3, MPI_FLOAT, MPI_STATUS_IGNORE);
        goffset += cellsinfo.globalcount.n * sizeof(float) * 3;
        MPI_File_seek(fhandle, goffset, MPI_SEEK_SET);

        /* write cell types */
        memset(buffer, 0, 256);
        sprintf(buffer, "\nCELL_TYPES %" PRId64 "\n", cellsinfo.globalcount.n);
        if (systeminfo.rank == 0) MPI_File_write(fhandle, &buffer, strlen(buffer), MPI_BYTE, MPI_STATUS_IGNORE);
        goffset += strlen(buffer);
        MPI_File_seek(fhandle, goffset, MPI_SEEK_SET);
        offset = nprev * sizeof(int);
        MPI_File_seek(fhandle, goffset+offset, MPI_SEEK_SET);
        for (i = 0; i < cellsinfo.localcount.n; i++)
                intscalarfield[i] = 1;
        if (systeminfo.endian) swapnbyte((char *) intscalarfield, cellsinfo.localcount.n, sizeof(int));
        MPI_File_write(fhandle, intscalarfield, cellsinfo.localcount.n, MPI_INT, MPI_STATUS_IGNORE);
        goffset += cellsinfo.globalcount.n * sizeof(int);
        MPI_File_seek(fhandle, goffset, MPI_SEEK_SET);

        /* write point data */
        memset(buffer, 0, 256);
        sprintf(buffer, "\nPOINT_DATA %" PRId64, cellsinfo.globalcount.n);
        if (systeminfo.rank == 0) MPI_File_write(fhandle, &buffer, strlen(buffer), MPI_BYTE, MPI_STATUS_IGNORE);
        goffset += strlen(buffer);
        MPI_File_seek(fhandle, goffset, MPI_SEEK_SET);

        /* write ranks */
        memset(buffer, 0, 256);
        sprintf(buffer, "\nSCALARS rank integer 1\nLOOKUP_TABLE default\n");
        if (systeminfo.rank == 0) MPI_File_write(fhandle, &buffer, strlen(buffer), MPI_BYTE, MPI_STATUS_IGNORE);
        goffset += strlen(buffer);
        MPI_File_seek(fhandle, goffset, MPI_SEEK_SET);
        for (i = 0; i < cellsinfo.localcount.n; i++)
                intscalarfield[i] = systeminfo.rank;
        offset = nprev * sizeof(int);
        MPI_File_seek(fhandle, goffset+offset, MPI_SEEK_SET);
        if (systeminfo.endian) swapnbyte((char *) intscalarfield, cellsinfo.localcount.n, sizeof(int));
        MPI_File_write(fhandle, intscalarfield, cellsinfo.localcount.n, MPI_INT, MPI_STATUS_IGNORE);
        goffset += cellsinfo.globalcount.n * sizeof(int);
        MPI_File_seek(fhandle, goffset, MPI_SEEK_SET);

        /* write cell type */
        memset(buffer, 0, 256);
        sprintf(buffer, "\nSCALARS type integer 1\nLOOKUP_TABLE default\n");
        if (systeminfo.rank == 0) MPI_File_write(fhandle, &buffer, strlen(buffer), MPI_BYTE, MPI_STATUS_IGNORE);
        goffset += strlen(buffer);
        MPI_File_seek(fhandle, goffset, MPI_SEEK_SET);
        for (i = 0; i < cellsinfo.localcount.n; i++)
                intscalarfield[i] = cellsinfo.cells[i].ctype;
        offset = nprev * sizeof(int);
        MPI_File_seek(fhandle, goffset+offset, MPI_SEEK_SET);
        if (systeminfo.endian) swapnbyte((char *) intscalarfield, cellsinfo.localcount.n, sizeof(int));
        MPI_File_write(fhandle, intscalarfield, cellsinfo.localcount.n, MPI_INT, MPI_STATUS_IGNORE);
        goffset += cellsinfo.globalcount.n * sizeof(int);
        MPI_File_seek(fhandle, goffset, MPI_SEEK_SET);

        /* write density */
        memset(buffer, 0, 256);
        sprintf(buffer, "\nSCALARS density float 1\nLOOKUP_TABLE default\n");
        if (systeminfo.rank == 0) MPI_File_write(fhandle, &buffer, strlen(buffer), MPI_BYTE, MPI_STATUS_IGNORE);
        goffset += strlen(buffer);
        MPI_File_seek(fhandle, goffset, MPI_SEEK_SET);
        for (i = 0; i < cellsinfo.localcount.n; i++)
                floatscalarfield[i] = 1.0; //cellsinfo.cells[i].density;
        offset = nprev * sizeof(float);
        MPI_File_seek(fhandle, goffset+offset, MPI_SEEK_SET);
        if (systeminfo.endian) swapnbyte((char *) floatscalarfield, cellsinfo.localcount.n, sizeof(float));
        MPI_File_write(fhandle, floatscalarfield, cellsinfo.localcount.n, MPI_FLOAT, MPI_STATUS_IGNORE);
        goffset += cellsinfo.globalcount.n * sizeof(float);
        MPI_File_seek(fhandle, goffset, MPI_SEEK_SET);

        /* write forces */
        sprintf(buffer, "\nVECTORS force float\n");
        if (systeminfo.rank == 0) MPI_File_write(fhandle, &buffer, strlen(buffer), MPI_BYTE, MPI_STATUS_IGNORE);
        goffset += strlen(buffer);
        MPI_File_seek(fhandle, goffset, MPI_SEEK_SET);
        for (i = 0; i < cellsinfo.localcount.n; i++) {
                floatvectorfield[3 * i] = cellsinfo.forces[i].x;
                floatvectorfield[3 * i + 1] = cellsinfo.forces[i].y;
                floatvectorfield[3 * i + 2] = cellsinfo.forces[i].z;
        }
        offset = nprev * sizeof(float) * 3;
        MPI_File_seek(fhandle, goffset+offset, MPI_SEEK_SET);
        if (systeminfo.endian) swapnbyte((char *) floatvectorfield, cellsinfo.localcount.n * 3, sizeof(float));
        MPI_File_write(fhandle, floatvectorfield, cellsinfo.localcount.n * 3, MPI_FLOAT, MPI_STATUS_IGNORE);
        goffset += cellsinfo.localcount.n * 3 * sizeof(float);
        MPI_File_seek(fhandle, goffset, MPI_SEEK_SET);


        MPI_File_close(&fhandle);

        free(floatvectorfield);
        free(floatscalarfield);
        free(intscalarfield);

        return;
}





/*!
 * This function defines output for global fields.
 */
void ioDefineOutputGlobalFields()
{
        int f;
        /* output fields */
/*        for (f = 0; f < NFIELDS; f++) {
                //if(f==BVES && !bvsim) continue;
                //if(f==TEMP && !temperature) continue;
                //if(f==OXYG && !oxygen) continue;
                //if(f==GLUC && !glucose) continue;
                //if(f==HYDR && !hydrogenIon) continue;
   //                strcpy(nameOut[f], fieldName[f]);
                dimOut[f] = SCALAR;
                typeOut[f] = REAL;
   //                addrOut[f] = &fieldAddr[f];
                jumpOut[f] = sizeof(double);
        }
        // output gradient
        for(f=NFIELDS; f<NFIELDS+NCHEM; f++) {
                //if(f-NFIELDS==OXYG && !oxygen) continue;
                //if(f-NFIELDS==GLUC && !glucose) continue;
                //if(f-NFIELDS==HYDR && !hydrogenIon) continue;
   //                strcpy(nameOut[f], fieldName[NGLOB+f-NFIELDS]);
                sprintf(nameOut[f]+strlen(nameOut[f]),"Gradient");
                dimOut[f] = VECTOR;
                typeOut[f] = REAL;
                //              addrOut[f] = &gradAddr[f];
                jumpOut[f] = sizeof(double);
        }*/
}

/*!
 * This function prints global fields data in VisNow data format
 * http://visnow.icm.edu.pl
 */
/*void ioWriteFields(int step)
   {
        int j,f;
        MPI_File fh1, fh2;
        FILE *fh3;
        float3dv_t *floatVectorField;
        float *floatScalarField;
        int *integerScalarField;
        int64_t size;
        int bdim;
        int gsize[3];
        int bsize[3];
        int bstart[3];
        MPI_Datatype subarray1_t, subarray2_t, float3_t;

        if (!gfields)
                return;

        ioDefineOutputGlobalFields();

        bdim = 3;
        gsize[0] = gridI;
        gsize[1] = gridJ;
        gsize[2] = gridK;
        bsize[0] = gridSize.x;
        bsize[1] = gridSize.y;
        bsize[2] = gridSize.z;
        bstart[0] = gridStartIdx[MPIrank].x;
        bstart[1] = gridStartIdx[MPIrank].y;
        bstart[2] = gridStartIdx[MPIrank].z;

        MPI_Type_vector(1, 3, 0, MPI_FLOAT, &float3_t);
        MPI_Type_commit(&float3_t);

        MPI_Type_create_subarray(bdim, gsize, bsize, bstart, MPI_ORDER_C,
                                 float3_t, &subarray1_t);
        MPI_Type_commit(&subarray1_t);

        gsize[0] = gridI;
        gsize[1] = gridJ;
        gsize[2] = gridK;
        bsize[0] = gridSize.x;
        bsize[1] = gridSize.y;
        bsize[2] = gridSize.z;
        bstart[0] = gridStartIdx[MPIrank].x;
        bstart[1] = gridStartIdx[MPIrank].y;
        bstart[2] = gridStartIdx[MPIrank].z;

        MPI_Type_create_subarray(bdim, gsize, bsize, bstart, MPI_ORDER_C,
                                 MPI_FLOAT, &subarray2_t);
        MPI_Type_commit(&subarray2_t);

        for (f = 0; f < NFIELDS+NCHEM; f++) {

                if(f==BVES && !bvsim) continue;
                if(f==TEMP && !temperature) continue;
                if((f==OXYG || f-NFIELDS+NGLOB==OXYG) && !oxygen) continue;
                if((f==GLUC || f-NFIELDS+NGLOB==GLUC) && !glucose) continue;
                if((f==HYDR || f-NFIELDS+NGLOB==HYDR) && !hydrogenIon) continue;

                char fstname1[256];
                char fstname2[256];
                char fstname3[256];

                size = gridSize.x * gridSize.y * gridSize.z;

                floatVectorField = (float3dv_t *) malloc(size * sizeof(float3dv_t));
                floatScalarField = (float *) malloc(size * sizeof(float));
                integerScalarField = (int *) malloc(size * sizeof(int));

                sprintf(fstname1, "%s/%s%08dcoords.bin", outdir, nameOut[f], step);
                sprintf(fstname2, "%s/%s%08dvalues.bin", outdir, nameOut[f], step);
                sprintf(fstname3, "%s/%s%08d.vnf", outdir, nameOut[f], step);

                if(MPIrank==0) {
                        fh3=fopen(fstname3,"w");
                        fprintf(fh3,"#VisNow regular field\n");
                        fprintf(fh3,"field %s%08d, dim %ld %ld %ld, coords\n",nameOut[f],step,gridI,gridJ,gridK);
                        if(dimOut[f]==SCALAR) fprintf(fh3,"component %s float\n",nameOut[f]);
                        else fprintf(fh3,"component %s float, vector 3\n",nameOut[f]);
                        fprintf(fh3,"file=%s%08dcoords.bin binary ",nameOut[f],step);
                        if(!endian) fprintf(fh3,"big\n");
                        else fprintf(fh3,"little\n");
                        fprintf(fh3,"coords\n");
                        fprintf(fh3,"file=%s%08dvalues.bin binary ",nameOut[f],step);
                        if(!endian) fprintf(fh3,"big\n");
                        else fprintf(fh3,"little\n");
                        fprintf(fh3,"%s\n",nameOut[f]);
                        fclose(fh3);
                }

                MPI_File_open(MPI_COMM_WORLD, fstname1,
                              MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh1);
                MPI_File_set_view(fh1, 0, MPI_FLOAT, subarray1_t, "native",
                                  MPI_INFO_NULL);
                // truncate the first file
                MPI_File_set_size(fh1, 0);

                for (j = 0; j < size; j++) {
                        floatVectorField[j].x = (float) (gridBuffer[j].x);
                        floatVectorField[j].y = (float) (gridBuffer[j].y);
                        floatVectorField[j].z = (float) (gridBuffer[j].z);
                }
                if (!endian)
                        swapnbyte((char *) floatVectorField, size * 3, sizeof(float));
                MPI_File_write(fh1, floatVectorField, 3 * size, MPI_FLOAT,
                               MPI_STATUS_IGNORE);
                MPI_File_close(&fh1);

                MPI_File_open(MPI_COMM_WORLD, fstname2,
                              MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh2);
                if(dimOut[f]==SCALAR)
                        MPI_File_set_view(fh2, 0, MPI_FLOAT, subarray2_t, "native",
                                          MPI_INFO_NULL);
                else
                        MPI_File_set_view(fh2, 0, MPI_FLOAT, subarray1_t, "native",
                                          MPI_INFO_NULL);

                // truncate the second file
                MPI_File_set_size(fh2, 0);

                if(dimOut[f]==SCALAR) {
                        for (j = 0; j < size; j++)
                                floatScalarField[j] = fieldAddr[f][j];
                        if (!endian)
                                swapnbyte((char *) floatScalarField, size, sizeof(float));
                        MPI_File_write(fh2, floatScalarField, size, MPI_FLOAT,
                                       MPI_STATUS_IGNORE);
                } else {
                        for (j = 0; j < size; j++) {
                                floatVectorField[j].x = (float) (gradAddr[f-NFIELDS][j*3]);
                                floatVectorField[j].y = (float) (gradAddr[f-NFIELDS][j*3+1]);
                                floatVectorField[j].z = (float) (gradAddr[f-NFIELDS][j*3+2]);
                        }
                        if (!endian)
                                swapnbyte((char *) floatVectorField, size*3, sizeof(float));
                        MPI_File_write(fh2, floatVectorField, size*3, MPI_FLOAT,
                                       MPI_STATUS_IGNORE);
                }

                MPI_File_close(&fh2);

                free(floatVectorField);
                free(floatScalarField);
                free(integerScalarField);

        }
        MPI_Type_free(&subarray1_t);
        MPI_Type_free(&subarray2_t);
   }*/

/*!
 * This function redirects stdout to a given file.
 */
void switchStdOut(const char *newStream)
{
        fflush(stdout);
        dup2(fileno(stdout), fdSave);
        fflush(stdout);
        fdNew = open(newStream, O_WRONLY | O_CREAT | O_TRUNC, 0644);
        dup2(fdNew, fileno(stdout));
        close(fdNew);
}

/*!
 * This function brings back the stdout.
 */
void revertStdOut()
{
        fflush(stdout);
        dup2(fdSave, fileno(stdout));
        close(fdSave);
}

/*!
 * This function defines colormaps for PovRay output.
 */
void defineColormaps()
{
        int numberOfColormaps;

        /* on the ocassion - initialize the rotation angle */
        beta = 0.0;

        numberOfColormaps = 6;
        cmaps = (colormap *) malloc(numberOfColormaps * sizeof(colormap));

        /* medical colormap i=0 */
        strcpy(cmaps[0].name, "medical");
        cmaps[0].ncp = 4;
        cmaps[0].cp =
                (colormapPoint *) malloc(cmaps[0].ncp * sizeof(colormapPoint));
        cmaps[0].cp[0].position = 0.0;
        cmaps[0].cp[0].r = 0.0;
        cmaps[0].cp[0].g = 0.0;
        cmaps[0].cp[0].b = 0.0;
        cmaps[0].cp[1].position = 0.33;
        cmaps[0].cp[1].r = 122.0;
        cmaps[0].cp[1].g = 32.0;
        cmaps[0].cp[1].b = 32.0;
        cmaps[0].cp[2].position = 0.66;
        cmaps[0].cp[2].r = 255.0;
        cmaps[0].cp[2].g = 179.0;
        cmaps[0].cp[2].b = 77.0;
        cmaps[0].cp[3].position = 1.0;
        cmaps[0].cp[3].r = 255.0;
        cmaps[0].cp[3].g = 255.0;
        cmaps[0].cp[3].b = 255.0;

        /* rainbow colormap i=1 */
        strcpy(cmaps[1].name, "rainbow");
        cmaps[1].ncp = 5;
        cmaps[1].cp =
                (colormapPoint *) malloc(cmaps[1].ncp * sizeof(colormapPoint));
        cmaps[1].cp[0].position = 0.0;
        cmaps[1].cp[0].r = 0.0;
        cmaps[1].cp[0].g = 0.0;
        cmaps[1].cp[0].b = 255.0;
        cmaps[1].cp[1].position = 0.25;
        cmaps[1].cp[1].r = 0.0;
        cmaps[1].cp[1].g = 255.0;
        cmaps[1].cp[1].b = 255.0;
        cmaps[1].cp[2].position = 0.5;
        cmaps[1].cp[2].r = 0.0;
        cmaps[1].cp[2].g = 255.0;
        cmaps[1].cp[2].b = 0.0;
        cmaps[1].cp[3].position = 0.75;
        cmaps[1].cp[3].r = 255.0;
        cmaps[1].cp[3].g = 255.0;
        cmaps[1].cp[3].b = 0.0;
        cmaps[1].cp[4].position = 1.0;
        cmaps[1].cp[4].r = 255.0;
        cmaps[1].cp[4].g = 0.0;
        cmaps[1].cp[4].b = 0.0;

        /* blue red yellow */
        strcpy(cmaps[2].name, "bry");
        cmaps[2].ncp = 4;
        cmaps[2].cp =
                (colormapPoint *) malloc(cmaps[2].ncp * sizeof(colormapPoint));
        cmaps[2].cp[0].position = 0.0;
        cmaps[2].cp[0].r = 0.0;
        cmaps[2].cp[0].g = 0.0;
        cmaps[2].cp[0].b = 255.0;
        cmaps[2].cp[1].position = 0.33;
        cmaps[2].cp[1].r = 255.0;
        cmaps[2].cp[1].g = 0.0;
        cmaps[2].cp[1].b = 255.0;
        cmaps[2].cp[2].position = 0.67;
        cmaps[2].cp[2].r = 255.0;
        cmaps[2].cp[2].g = 0.0;
        cmaps[2].cp[2].b = 0.0;
        cmaps[2].cp[3].position = 1.0;
        cmaps[2].cp[3].r = 255.0;
        cmaps[2].cp[3].g = 255.0;
        cmaps[2].cp[3].b = 0.0;

        /* hot */
        strcpy(cmaps[3].name, "hot");
        cmaps[3].ncp = 5;
        cmaps[3].cp =
                (colormapPoint *) malloc(cmaps[3].ncp * sizeof(colormapPoint));
        cmaps[3].cp[0].position = 0.0;
        cmaps[3].cp[0].r = 107.0;
        cmaps[3].cp[0].g = 0.0;
        cmaps[3].cp[0].b = 0.0;
        cmaps[3].cp[1].position = 0.35;
        cmaps[3].cp[1].r = 255.0;
        cmaps[3].cp[1].g = 102.0;
        cmaps[3].cp[1].b = 28.0;
        cmaps[3].cp[2].position = 0.57;
        cmaps[3].cp[2].r = 250.0;
        cmaps[3].cp[2].g = 235.0;
        cmaps[3].cp[2].b = 128.0;
        cmaps[3].cp[3].position = 0.76;
        cmaps[3].cp[3].r = 232.0;
        cmaps[3].cp[3].g = 230.0;
        cmaps[3].cp[3].b = 230.0;
        cmaps[3].cp[4].position = 1.0;
        cmaps[3].cp[4].r = 156.0;
        cmaps[3].cp[4].g = 161.0;
        cmaps[3].cp[4].b = 255.0;

        /* hot1 */
        strcpy(cmaps[4].name, "hot1");
        cmaps[4].ncp = 5;
        cmaps[4].cp =
                (colormapPoint *) malloc(cmaps[4].ncp * sizeof(colormapPoint));
        cmaps[4].cp[0].position = 0.0;
        cmaps[4].cp[0].r = 128.0;
        cmaps[4].cp[0].g = 0.0;
        cmaps[4].cp[0].b = 0.0;
        cmaps[4].cp[1].position = 0.2;
        cmaps[4].cp[1].r = 255.0;
        cmaps[4].cp[1].g = 0.0;
        cmaps[4].cp[1].b = 0.0;
        cmaps[4].cp[2].position = 0.4;
        cmaps[4].cp[2].r = 255.0;
        cmaps[4].cp[2].g = 255.0;
        cmaps[4].cp[2].b = 0.0;
        cmaps[4].cp[3].position = 0.7;
        cmaps[4].cp[3].r = 255.0;
        cmaps[4].cp[3].g = 255.0;
        cmaps[4].cp[3].b = 255.0;
        cmaps[4].cp[4].position = 1.0;
        cmaps[4].cp[4].r = 128.0;
        cmaps[4].cp[4].g = 128.0;
        cmaps[4].cp[4].b = 255.0;

        /* my */
        strcpy(cmaps[5].name, "my");
        cmaps[5].ncp = 5;
        cmaps[5].cp =
                (colormapPoint *) malloc(cmaps[5].ncp * sizeof(colormapPoint));
        cmaps[5].cp[0].position = 0.0;
        cmaps[5].cp[0].r = 107.0;
        cmaps[5].cp[0].g = 0.0;
        cmaps[5].cp[0].b = 0.0;
        cmaps[5].cp[1].position = 0.35;
        cmaps[5].cp[1].r = 0.0;
        cmaps[5].cp[1].g = 100.0;
        cmaps[5].cp[1].b = 255.0;
        cmaps[5].cp[2].position = 0.57;
        cmaps[5].cp[2].r = 250.0;
        cmaps[5].cp[2].g = 235.0;
        cmaps[5].cp[2].b = 128.0;
        cmaps[5].cp[3].position = 0.76;
        cmaps[5].cp[3].r = 232.0;
        cmaps[5].cp[3].g = 230.0;
        cmaps[5].cp[3].b = 230.0;
        cmaps[5].cp[4].position = 1.0;
        cmaps[5].cp[4].r = 156.0;
        cmaps[5].cp[4].g = 161.0;
        cmaps[5].cp[4].b = 255.0;

}


/*!
 * This function prints PovRay output with cellular data.
 */
/*
   void ioWriteStepPovRay(int step, int type)
   {
        int c;
        int i;
        char fstname[256];
        MPI_File fh;
        MPI_Status status;
        MPI_Datatype subarray_t;
        float cr, cg, cb;
        char *const fmt =
                "sphere{ <%10.4lf,%10.4lf,%10.4lf>,%10.4lf texture { pigment { color rgb <%6.4f,%6.4f,%6.4f> } finish { phong 0.2 ambient .1 }}}         \n";
        char testBuffer[512];
        char *txtData;
        char *txtData_p;
        char *txtHeaderData;
        int numCharsPerCell;
        int headerLen;
        int headerSize;
        int error;
        int gdims;
        int gsize[1];
        int istart[1];
        int isize[1];
        int64_t printed = 0;
        int64_t tPrinted[MPIsize];
        double fmin, fmax, fepsilon;
        int cm;
        int cmReverse = 0;
        int cmShift = 0;
        double cmPad;
        double minCorner[3], maxCorner[3];
        double gMinCorner[3], gMaxCorner[3];

        double middlePointLocal[3];
        double middlePointGlobal[3];
        double lmass, gmass;

        // type: 0 - denisty, 1 - oxygen, 2 - phases, 3 - slice & phases

        minCorner[0] = DBL_MAX;
        minCorner[1] = DBL_MAX;
        minCorner[2] = DBL_MAX;
        maxCorner[0] = -DBL_MAX;
        maxCorner[1] = -DBL_MAX;
        maxCorner[2] = -DBL_MAX;

        for (i = 0; i < lnc; i++) {
                minCorner[0] =
                        (cells[i].x - cells[i].size <
                         minCorner[0] ? cells[i].x - cells[i].size : minCorner[0]);
                maxCorner[0] =
                        (cells[i].x + cells[i].size >
                         maxCorner[0] ? cells[i].x + cells[i].size : maxCorner[0]);
                minCorner[1] =
                        (cells[i].y - cells[i].size <
                         minCorner[1] ? cells[i].y - cells[i].size : minCorner[1]);
                maxCorner[1] =
                        (cells[i].y + cells[i].size >
                         maxCorner[1] ? cells[i].y + cells[i].size : maxCorner[1]);
                minCorner[2] =
                        (cells[i].z - cells[i].size <
                         minCorner[2] ? cells[i].z - cells[i].size : minCorner[2]);
                maxCorner[2] =
                        (cells[i].z + cells[i].size >
                         maxCorner[2] ? cells[i].z + cells[i].size : maxCorner[2]);
        }
        MPI_Allreduce(minCorner, gMinCorner, 3, MPI_DOUBLE, MPI_MIN,
                      MPI_COMM_WORLD);
        MPI_Allreduce(maxCorner, gMaxCorner, 3, MPI_DOUBLE, MPI_MAX,
                      MPI_COMM_WORLD);

        middlePointLocal[0] = 0.0;
        middlePointLocal[1] = 0.0;
        middlePointLocal[2] = 0.0;
        lmass = 0.0;
        gmass = 0.0;

        for (c = 0; c < lnc; c++) {
                middlePointLocal[0] += cells[c].size * cells[c].x;
                middlePointLocal[1] += cells[c].size * cells[c].y;
                middlePointLocal[2] += cells[c].size * cells[c].z;
                lmass += cells[c].size;
        }

        MPI_Allreduce(middlePointLocal, middlePointGlobal, 3, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&lmass, &gmass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        middlePointGlobal[0] /= gmass;
        middlePointGlobal[1] /= gmass;
        middlePointGlobal[2] /= gmass;

        middlePointGlobal[0] =
                gMinCorner[0] + (gMaxCorner[0] - gMinCorner[0]) / 2;
        middlePointGlobal[1] =
                gMinCorner[1] + (gMaxCorner[1] - gMinCorner[1]) / 2;
        middlePointGlobal[2] =
                gMinCorner[2] + (gMaxCorner[2] - gMinCorner[2]) / 2;

        // write data to test buffer
        cr = 0.0;
        cg = 0.0;
        cb = 0.0;
        cm = 0;
        numCharsPerCell =
                sprintf(testBuffer, fmt, cells[1].x, cells[1].y, cells[1].z,
                        cells[1].size, cr, cg, cb);
        if (!
            (txtData =
                     (char *) malloc(numCharsPerCell * (lnc + 16) * sizeof(char))))
                stopRun(106, "txtData", __FILE__, __LINE__);
        txtData_p = txtData;

        for (i = 0; i < MPIsize; i++)
                tPrinted[i] = 0;

        switch (type) {
        case 0:
                sprintf(fstname, "%s/step%08ddensity.pov", outdir, step);
                cm = 5;
                cmReverse = 1;
                break;
        case 1:
                sprintf(fstname, "%s/step%08doxygen.pov", outdir, step);
                MPI_Allreduce(&fieldMin[OXYG], &fmin, 1, MPI_DOUBLE, MPI_MIN,
                              MPI_COMM_WORLD);
                MPI_Allreduce(&fieldMax[OXYG], &fmax, 1, MPI_DOUBLE, MPI_MAX,
                              MPI_COMM_WORLD);
                fepsilon = (fmax - fmin) * 0.1;
                fmin -= fepsilon;
                fmax += fepsilon;
                cm = 0;
                break;
        case 2:
                sprintf(fstname, "%s/step%08dphases.pov", outdir, step);
                cm = 1;
                break;
        case 3:
                sprintf(fstname, "%s/step%08dslice.pov", outdir, step);
                cm = 1;
                break;
        case 4:
                sprintf(fstname, "%s/step%08doutside.pov", outdir, step);
                cm = 0;
                cmReverse = 1;
                cmShift = 1;
                cmPad = 0.15;
                break;
        }

        // open file
        error =
                MPI_File_open(MPI_COMM_WORLD, fstname,
                              MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
        if (error != MPI_SUCCESS)
                if (MPIrank == 0)
                        stopRun(113, NULL, __FILE__, __LINE__);

        MPI_File_set_size(fh, 0);

        // write header
        if (MPIrank == 0) {

                float cameraLocation[3];
                float lookAt[3];
                float lightSource1[3];
                float lightSource2[3];
                float a1, a2, a3, b1, b2, b3, dist;
                float ss, cc;
                float corner;

                a1 = (gMaxCorner[0] - gMinCorner[0]) / 2;
                a2 = (gMaxCorner[1] - gMinCorner[1]) / 2;
                a3 = (gMaxCorner[2] - gMinCorner[2]) / 2;

                b1 = a1 / (tan(30.0 * M_PI / 180.0));
                b2 = a2 / (tan(30.0 * M_PI / 180.0));
                b3 = a3 / (tan(30.0 * M_PI / 180.0));

                dist = (b1 >= b2 ? b1 : b2);
                dist = (dist >= b3 ? dist : b3);

                ss = sin(beta * M_PI / 180.0);
                cc = cos(beta * M_PI / 180.0);

                corner = sqrt(pow(a1, 2) + pow(a3, 2));

                cameraLocation[0]=middlePointGlobal[0]-2;
                cameraLocation[1]=middlePointGlobal[1]-0.8*dist;
                cameraLocation[2]=middlePointGlobal[2]-5;

                //cameraLocation[0] =
                //   -(middlePointGlobal[2] - corner - dist -
                //   middlePointGlobal[2]) * ss + middlePointGlobal[0];
                //   cameraLocation[1] = middlePointGlobal[1];
                //   cameraLocation[2] =
                //   (middlePointGlobal[2] - corner - dist -
                //   middlePointGlobal[2]) * cc + middlePointGlobal[2];

                lookAt[0] = middlePointGlobal[0]-2;
                lookAt[1] = middlePointGlobal[1];
                lookAt[2] = middlePointGlobal[2]-5;

                //   lookAt[0] = middlePointGlobal[0];
                //   lookAt[1] = middlePointGlobal[1];
                //   lookAt[2] = middlePointGlobal[2];

                lightSource1[0] = cameraLocation[0];
                lightSource1[1] = cameraLocation[1];
                lightSource1[2] = cameraLocation[2];

                lightSource2[0] = lightSource1[0];
                lightSource2[1] = lightSource1[1];
                lightSource2[2] = lightSource1[2];

                txtHeaderData = (char *) malloc(1024 * sizeof(char));

                headerLen =
                        sprintf(txtHeaderData,
                                "#include \"colors.inc\"\ncamera { location <%f,%f,%f> look_at <%f,%f,%f> angle 60}\nlight_source { <%f,%f,%f> color White }\nlight_source { <%f,%f,%f> color White }\nbackground { color White }\n",
                                cameraLocation[0], cameraLocation[1], cameraLocation[2],
                                lookAt[0], lookAt[1], lookAt[2], lightSource1[0],
                                lightSource1[1], lightSource1[2], lightSource2[0],
                                lightSource2[1], lightSource2[2]);

                headerSize = headerLen * sizeof(char);

                MPI_File_write(fh, txtHeaderData, headerLen, MPI_CHAR, &status);

                free(txtHeaderData);
        }

        MPI_Bcast(&headerSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

        for (c = 0; c < lnc; c++) {
                float color = 0.0;
                int jump;
                if(beta==360.0 && (cells[c].y<1.8 || cells[c].y>4.5)) continue;
                if (beta == 360.0
                    && (cells[c].z > middlePointGlobal[2] + 4.0 * csize
 || cells[c].z < middlePointGlobal[2] - 4.0 * csize))
                        continue;
                if (type == 0) {
                        color = ((cells[c].density) / (8*densityCriticalLevel2));                                                                                       // UWAGA!
                        if (cells[c].tumor)
                                color = 0.0;
                }                                                                                                   //2.5);
                if (type == 1)
                        color = ((cellFields[OXYG][c] - fmin) / (fmax - fmin));
                if (type == 2 || type == 3 || type == 4)
                        color = (((float) cells[c].phase) / 5.0);
                if (type == 2)
                        color = (((float) MPIrank) / 512.0);

                if (type == 0)
                        color = color / 4 + 0.75;

                if (cmReverse)
                        color = 1.0 - color;
                if (cmShift) {
                        color = color * (1 - cmPad);
                        if (color == 1 - cmPad)
                                color = 1.0;
                }

                if(color<0.0) color=0.0;

                if(cells[c].ctype==0) color=0.2;
                else color=0.0;

                for (i = 1; i < cmaps[cm].ncp; i++) {
                        float d, dr, dg, db;
                        d = cmaps[cm].cp[i].position - cmaps[cm].cp[i - 1].position;
                        dr = cmaps[cm].cp[i].r - cmaps[cm].cp[i - 1].r;
                        dg = cmaps[cm].cp[i].g - cmaps[cm].cp[i - 1].g;
                        db = cmaps[cm].cp[i].b - cmaps[cm].cp[i - 1].b;
                        if (color <= cmaps[cm].cp[i].position) {
                                cr = cmaps[cm].cp[i - 1].r +
                                     ((color - cmaps[cm].cp[i - 1].position) / d) * dr;
                                cg = cmaps[cm].cp[i - 1].g +
                                     ((color - cmaps[cm].cp[i - 1].position) / d) * dg;
                                cb = cmaps[cm].cp[i - 1].b +
                                     ((color - cmaps[cm].cp[i - 1].position) / d) * db;
                                break;
                        }

                }

                cr /= 255.0;
                cg /= 255.0;
                cb /= 255.0;

                jump =
                        sprintf(txtData_p, fmt, cells[c].x, cells[c].y, cells[c].z,
                                cells[c].size, cr, cg, cb);
                txtData_p += jump;
                printed += 1;

        }

        if (printed == 0) {
                float coords[3];
                float size = 0.0;
                float color[3];
                int jump;
                // artificial cell far away from the scene
                coords[0] = -512.0;
                coords[1] = -512.0;
                coords[2] = -512.0;
                color[0] = 0.0;
                color[1] = 0.0;
                color[2] = 0.0;

                printed = 1;

                jump =
                        sprintf(txtData_p, fmt, coords[0], coords[1], coords[2], size,
                                color[0], color[1], color[2]);
                txtData_p += jump;
        }

        gdims = 1;
        tPrinted[MPIrank] = printed;
        MPI_Allgather(&printed, 1, MPI_INT64_T, tPrinted, 1, MPI_INT64_T,
                      MPI_COMM_WORLD);

        gsize[0] = 0;
        istart[0] = 0;
        isize[0] = printed * numCharsPerCell;
        for (i = 0; i < MPIrank; i++)
                istart[0] += tPrinted[i] * numCharsPerCell;

        for (i = 0; i < MPIsize; i++)
                gsize[0] += tPrinted[i] * numCharsPerCell;

        MPI_Type_create_subarray(gdims, gsize, isize, istart, MPI_ORDER_C,
                                 MPI_CHAR, &subarray_t);
        MPI_Type_commit(&subarray_t);

        MPI_File_set_view(fh, headerSize, MPI_CHAR, subarray_t, "native",
                          MPI_INFO_NULL);

        MPI_File_write(fh, txtData, isize[0], MPI_CHAR, &status);

        free(txtData);

        MPI_File_close(&fh);
        MPI_Type_free(&subarray_t);

        if (beta < 360.0)
                beta += 1.0;
        else
                beta = 360.0;

   }
 */
