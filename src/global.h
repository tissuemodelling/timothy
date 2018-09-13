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

#include <zoltan.h>
#include "_hypre_utilities.h"
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_krylov.h"
#include <stdbool.h>

/*! \file global.h
 *  \brief contains the most important global variables, arrays and defines
 */

#define  VERSION   "1.01"

#define POWER_OF_TWO(x) !(x&(x-1))

/* architecture */
#if defined(__ia64) || defined(__itanium__) || defined(_M_IA64)
#define CPUARCH "Itanium"
#endif
#if defined(__powerpc__) || defined(__ppc__) || defined(__PPC__)
#if defined(__powerpc64__) || defined(__ppc64__) || defined(__PPC64__) || \
	defined(__64BIT__) || defined(_LP64) || defined(__LP64__)
#define CPUARCH "POWER64"
#else
#define CPUARCH "POWER32"
#endif
#endif
#if defined(__sparc)
#define CPUARCH "SPARC"
#endif
#if defined(__x86_64__) || defined(_M_X64)
#define CPUARCH "x86_64"
#elif defined(__i386) || defined(_M_IX86)
#define CPUARCH "x86_32"
#endif

#define MIN_CELLS_PER_PROC 16

#define BOXSIZEX 2048
#define BOXSIZEY 2048
#define BOXSIZEZ 2048

#define ZOLTAN_ID_TYPE int

#define N_LEVELS 30
#define ROOT_LEVEL N_LEVELS-1
#define MAXVALUE powf(2,ROOT_LEVEL)

#define CELLTYPE_G1_DEFAULT 11.0
#define CELLTYPE_S_DEFAULT 8.0
#define CELLTYPE_G2_DEFAULT 4.0
#define CELLTYPE_M_DEFAULT 1.0
#define CELLTYPE_V_DEFAULT 0.5
#define CELLTYPE_RD_DEFAULT 0.1
#define CELLTYPE_CDENS_DEFAULT 2
#define CELLTYPE_PROD_DEFAULT 0.0
#define CELLTYPE_CONS_DEFAULT 0.0
#define CELLTYPE_CL1_DEFAULT 100
#define CELLTYPE_CL2_DEFAULT 100
#define CELLTYPE_SIZE_DEFAULT 1.0
#define CELLTYPE_H_DEFAULT 1.5
#define NUMBER_OF_CELLENV_PAR 4

#define ENVIRONMENT_DC_DEFAULT 1.82e-5
#define ENVIRONMENT_BC_DEFAULT 0.1575e-6
#define ENVIRONMENT_ICMEAN_DEFAULT 0.1575e-6
#define ENVIRONMENT_ICVAR_DEFAULT 0.0
#define ENVIRONMENT_LAMBDA_DEFAULT 0.0

typedef struct cellcount_t{
  uint64_t n;
  uint64_t g0phase;
  uint64_t g1phase;
  uint64_t sphase;
  uint64_t g2phase;
  uint64_t mphase;
  uint64_t necroticphase;
} cellcount_t;

typedef struct celldata_t {
  int lifetime;         /* age of the cell */
  int phase;            /* actual phase of the cell (0=G0,1=G1,2=S,3=G2,4=M,5=Necrotic) */
  int age;              /* cell's age */
  int death;
  float phasetime;      /* actual phase time */
  float g1;
  float s;
  float g2;
  float m;
  float young;
  ZOLTAN_ID_TYPE gid;    /* global ID of the particle */
  double x;              /* x coordinate of the particle position */
  double y;              /* y coordinate of the particle position */
  double z;              /* z coordinate of the particle position */
  double size;           /* radius of the cell */
  double v;              /* particle potential */
  double density;        /* particle density */
	double mindist;
  int ctype;
} celldata_t;

typedef struct cellenvdata_t {
	double value;
	double gx;
	double gy;
	double gz;
} cellenvdata_t;

typedef struct double3dv_t {
  double x;
  double y;
  double z;
} double3dv_t;

typedef struct float3dv_t {
  float x;
  float y;
  float z;
} float3dv_t;

typedef struct int643dv_t {
	int64_t x;
	int64_t y;
	int64_t z;
} int643dv_t;

typedef struct uint3dv_t {
	unsigned int x;
  unsigned int y;
  unsigned int z;
} uint3dv_t;

typedef struct octnode_t {
  unsigned int xcode;
  unsigned int ycode;
  unsigned int zcode;
  unsigned int xlimit;
  unsigned int ylimit;
  unsigned int zlimit;
  unsigned int level;
  int64_t father;
  int64_t child[8];
  int data;
} octnode_t;

typedef struct cellsinfo_t{
	cellcount_t localcount;
	cellcount_t globalcount;
	cellcount_t *localtypecount;
	cellcount_t *globaltypecount;
	uint64_t *cellsperproc;
	celldata_t *cells;
	double3dv_t *forces;
	octnode_t *octree;
	int64_t octsize;
	int64_t octmaxsize;
	int dimension;
	double3dv_t *velocity;
} cellsinfo_t;

/* NEW */
typedef struct grid_t{
  int643dv_t globalsize;
  double3dv_t lowercorner,uppercorner;
  int643dv_t localsize;
  float resolution;
  int643dv_t *loweridx,*upperidx;
	double3dv_t *lowleftnear,*uprightfar;
  //double3dv_t *data;
} grid_t;

typedef struct expinterpdata_t {
	double x;
	double y;
	double z;
	int ctype;
} expinterpdata_t;

typedef struct interpdata_t {
	MPI_Request *reqsend;
	MPI_Request *reqrecv;
	int *recvcount;
	int *sendcount;
	expinterpdata_t *sendinterpdata;
	expinterpdata_t *recvinterpdata;
	int numexp;
	int numimp;
} interpdata_t;

/* NEW */

typedef struct explist_t {
        int64_t cell;
        int proc;
} explist_t;

typedef struct fieldindata_t {
	double x;
	double y;
	double z;
	int ctype;
} fieldindata_t;

typedef struct fieldoutdata_t {
	double *fdata;
} fieldoutdata_t;

typedef struct fieldcommdata_t {
	explist_t *explist;
	int explistmaxsize;
	MPI_Request *reqsend;
	MPI_Request *reqrecv;
	int64_t *sendoffset;
	int64_t *recvoffset;
	int *recvcount;
	int *sendcount;
	fieldindata_t *sendfieldindata;
	fieldindata_t *recvfieldindata;
	fieldoutdata_t *sendfieldoutdata;
	fieldoutdata_t *recvfieldoutdata;
	int numexp;
	int numimp;
} fieldcommdata_t;



typedef struct cellindata_t { /* this structure keeps cell data needed in potential computations */
  double x;
  double y;
  double z;
  double size;
  double young;
  int ctype;
} cellindata_t;

typedef struct celloutdata_t { /* this structure keeps additional cell data (potential & density) */
  double v;
  double density;
} celloutdata_t;

typedef struct cellcommdata_t {
	explist_t *explist;
	int explistmaxsize;
	MPI_Request *reqsend;
	MPI_Request *reqrecv;
	int64_t *sendoffset;
	int64_t *recvoffset;
	int *recvcount;
	int *sendcount;
	cellindata_t *sendcellindata;
	cellindata_t *recvcellindata;
	celloutdata_t *sendcelloutdata;
	celloutdata_t *recvcelloutdata;
	int numexp;
	int numimp;
} cellcommdata_t;

typedef struct settings_t {
 		int64_t maxcells;   /* maximal number of cells (set in parameter file) */
		int64_t maxlocalcells;
		int numberofsteps;
		float secondsperstep;
		int numberofcelltypes;
		int numberoffields;
		int dimension;
		int restart;
    char rstfilename[128];
		char outdir[128];
    int visoutstep;
    int statoutstep;
		int rstoutstep;
		float maxspeed;
		float gfdt;
		float gfh;
		int simulationstart;
		unsigned int rseed;
		int step;
} settings_t;

typedef struct systeminfo_t {
    int rank;                    /* MPI rank */
    int size;                    /* MPI size */
    int nthreads;
    int dim[3];
    int noderank;
    int nodesize;
    int memperproc;
    MPI_Comm MPI_CART_COMM;
    int **coords;
		int endian;
    //int restart;
} systeminfo_t;


typedef struct celltype_t {
	char name[128];
	float g1;
	float s;
	float g2;
	float m;
	float v;
	float rd;
	double size;
	double h;			/* neighbourhood of the cell */
	double h2;		/* h^2 */
	double h3;		/* h^3 */
	double h4;		/* h^4 */
	char inputfile[128];
	float criticaldensity;
	float *production;
	float *consumption;
	float *criticallevel1;
	float *criticallevel2;
	float randommove;
} celltype_t;

typedef struct patches_t {
	double **data;
	int *intersect;
	int *receiver;
	int *sender;
	double **recvdata;
	MPI_Request *reqsend;
	MPI_Request *reqrecv;
	int643dv_t *lowercorner;
	int643dv_t *uppercorner;
	int643dv_t *lowercornerR;
	int643dv_t *uppercornerR;
	int643dv_t *size;
	double **commbuff;
	double **buff;
} patches_t;

typedef struct environment_t {
	char name[128];
	double diffusioncoefficient;
	double boundarycondition;
	double initialconditionmean;
	double initialconditionvariance;
	double lambdadelay;
	double *data;
	double *gradient;
	double *production;
} environment_t;

#ifdef HYPRE
typedef struct solversettings_t {
	double *z;
	double *dt;
	int enviter;
	int numberofiters;
	int envobjecttype;
  HYPRE_SStructVariable *vartypes;
	HYPRE_SStructStencil *stencil;
	HYPRE_Solver precond;
	HYPRE_SStructGrid grid;
	HYPRE_SStructGraph graph;
	HYPRE_Int lower[3];
	HYPRE_Int upper[3];
	HYPRE_Int bclower[3];
	HYPRE_Int bcupper[3];
} solversettings_t;

typedef struct solverdata_t {
	HYPRE_SStructMatrix A;
	HYPRE_SStructVector b;
	HYPRE_SStructVector x;
	HYPRE_ParCSRMatrix parA;
	HYPRE_ParVector parb;
	HYPRE_ParVector parx;
	HYPRE_Solver solver;
} solverdata_t;
#endif

/* statistics */
typedef struct statistics_t {
  double minsize; /* Minimum size of cells */
  double maxsize; /* Maximum size of cells */
  double mindist; /* Minimum distance between neighbourhood cells */
  double maxspeed;  /* Maximum speed in the systeminfo */
  double minspeed;  /* Minimum speed in the systeminfo */
  double maxdens; /* Maximum density */
  double mindens; /* Minimum density */
  double densdev; /* Density deviation */
  double densavg; /* Average density */
} statistics_t;

typedef struct octheap_t {
  int size;
  int count;
  int *data;
} octheap_t;
