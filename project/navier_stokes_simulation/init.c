#include "helper.h"
#include "init.h"
#include "ns_definitions.h"
#include "boundary_val.h"
#include "parallel.h"

#include <string.h>
#include <mpi.h>

int read_parameters(const char * szFileName,
                    double *UI,                /* velocity x-direction */
                    double *VI,                /* velocity y-direction */
                    double *PI,                /* pressure */
                    double *GX,                /* gravitation x-direction */
                    double *GY,                /* gravitation y-direction */
                    double *t_end,             /* end time */
                    double *xlength,           /* length of the domain x-dir.*/
                    double *ylength,           /* length of the domain y-dir.*/
                    double *dt,                /* time step */
                    double *dx,                /* length of a cell x-dir. */
                    double *dy,                /* length of a cell y-dir. */
                    int  *imax,                /* number of cells x-direction*/
                    int  *jmax,                /* number of cells y-direction*/
                    double *alpha,             /* uppwind differencing factor*/
                    double *omg,               /* relaxation factor */
                    double *tau,               /* safety factor for time step*/
                    int  *itermax,             /* max. number of iterations  */
		                               	   	   /* for pressure per time step */
                    double *eps,               /* accuracy bound for pressure*/
                    double *dt_value,			/* time for output */
                    int *wl,			/* domain boundary conditions */
                    int *wr,			/* for right, left, top and */
                    int *wt,			/* bottom surfaces. */
                    int *wb,
                    int *iproc,					/* i-dim of processors */
                    int *jproc					/* j-dim of processors */
)
{
	/* Read domain boundary conditions*/
	READ_INT   ( szFileName, *wl );
	READ_INT   ( szFileName, *wr );
	READ_INT   ( szFileName, *wt );
	READ_INT   ( szFileName, *wb );

	/* Read all parameters */
	READ_DOUBLE( szFileName, *xlength );
	READ_DOUBLE( szFileName, *ylength );

	READ_DOUBLE( szFileName, *t_end );
	READ_DOUBLE( szFileName, *dt    );

//	READ_INT   ( szFileName, *imax );
//	READ_INT   ( szFileName, *jmax );

	READ_DOUBLE( szFileName, *omg   );
	READ_DOUBLE( szFileName, *eps   );
	READ_DOUBLE( szFileName, *tau   );
	READ_DOUBLE( szFileName, *alpha );

	READ_INT   ( szFileName, *itermax );
	READ_DOUBLE( szFileName, *dt_value );

	READ_DOUBLE( szFileName, *UI );
	READ_DOUBLE( szFileName, *VI );
	READ_DOUBLE( szFileName, *GX );
	READ_DOUBLE( szFileName, *GY );
	READ_DOUBLE( szFileName, *PI );

	READ_INT   ( szFileName, *iproc );
	READ_INT   ( szFileName, *jproc );

	*dx = *xlength / (double)(*imax);
	*dy = *ylength / (double)(*jmax);

	return 1;
}


/*
 * Initialize the U, V, and P grids to what UI, VU, and PI
 * are set to, respectively.
 */
void init_uvp(
	double UI,
	double VI,
	double PI,
	int imax,
	int jmax,
	double ***U,
	double ***V,
	double ***P
)
{
	/* Horizontal velocity */
	*U = matrix(0, imax+1, 0, jmax+1);
	init_matrix(*U, 0, imax+1, 0, jmax+1, UI);

	/* Vertical velocity */
	*V = matrix(0, imax+1, 0, jmax+1);
	init_matrix(*V, 0, imax+1, 0, jmax+1, VI);

	/* Pressure */
	*P = matrix(0, imax+1, 0, jmax+1);
	init_matrix(*P, 0, imax+1, 0, jmax+1, PI);
}

void init_flag(const char *problem, int xdim, int ydim, int imax, int jmax, int il, int ir, int jb, int jt, int ***flag)
{
	char *problemPBMFile = 0;
	int **obstacleFlag = 0;
	int i, j;

	/* initialize Flag as a matrix of integers */
	*flag = imatrix(0, xdim+1, 0, ydim+1);

	/* Retrieve full file name */
	problemPBMFile = malloc(strlen(problem) + 5);
	strcpy(problemPBMFile, problem);
	strcat(problemPBMFile, ".pbm");

	/* Input file containing fluid/obstacle data */
	obstacleFlag = read_pgm(problemPBMFile);

	/* Set the boundary values for flag */
	boundaryvalues_flag(imax, jmax, il, ir, jb, jt, flag, obstacleFlag);

	/* Set inner flags */
	for(i = 1; i <= xdim; i++)
	{
		for(j = 1; j <= ydim; j++)
		{
			/* Determine whether cell is a fluid or an obstacle */
			if(obstacleFlag[il+i-1][jb+j-1] == 1)
			{
				(*flag)[i][j] = C_F;  /* Fluid */
			}
			else
			{
				(*flag)[i][j] = C_B;  /* Obstacle/border */
			}

			/* Set flags according to neighboring fluids */
			if(obstacleFlag[il+i-2][jb+j-1] == 1)
			{
				(*flag)[i][j] |= B_W;
			}
			if(obstacleFlag[il+i][jb+j-1] == 1)
			{
				(*flag)[i][j] |= B_O;
			}
			if(obstacleFlag[il+i-1][jb+j] == 1)
			{
				(*flag)[i][j] |= B_N;
			}
			if(obstacleFlag[il+i-1][jb+j-2] == 1)
			{
				(*flag)[i][j] |= B_S;
			}
		}
	}

	free_imatrix(obstacleFlag, 0, xdim+1, 0, ydim+1);
	free(problemPBMFile);
}

void broadcast_parameters(
	int myrank,
	double *Re,
	double *UI,
	double *VI,
	double *PI,
	double *GX,
	double *GY,
	double *t_end,
	double *xlength,
	double *ylength,
	double *dt,
	double *dx,
	double *dy,
	int  *imax,
	int  *jmax,
	double *alpha,
	double *omg,
	double *tau,
	int  *itermax,
	double *eps,
	double *dt_value,
	int *wl,
	int *wr,
	int *wt,
	int *wb,
	int *iproc,
	int *jproc)
{
	/* Broadcast parameters to the remaining processes */
//	MPI_Bcast(Re, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(UI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(VI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(PI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(GX, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(GY, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(t_end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(xlength, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(ylength, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(dx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(dy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(imax, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(jmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(omg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(itermax, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(eps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(dt_value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(wl, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(wr, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(wt, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(wb, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(iproc, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(jproc, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if(myrank == 0)
	{
		/*printf("All simulation parameters have been broadcasted!\n\n");*/
	}
}

int read_args(int argc, char** argv, int *flag_Re, double* Re, double *viscosity, int* mc_id, int* imax, int* jmax, double *rho)
{
	if(argc == 6)
	{
		/*	flag_Re={1, or something else}, if 1 then argv[1] is reynold's number otherwise it is Viscosity*/
		*flag_Re = atoi(argv[1]);

		if (*flag_Re == 1)
		{
			*Re = atof(argv[2]);
			if(*Re <= 0)
			{
				Programm_Stop("Reynolds number must be greater than zero.");
				return 0;
			}
		}
		else
		{	
			*viscosity = atof(argv[2]);
			if(*viscosity <= 0)
			{
				Programm_Stop("Viscosity number must be greater than zero.");
				return 0;
			}
			/*	Re is calculated based on the channel width (height of step can also be cansidered?),
			 *	density(rho) specified in the data file and the initial flow velocity.
			 *	Re=100 corresponds to 0.0245 kg/(m·s)
			 */
			*Re = (2 * (*rho) * 1)/(*viscosity);
		}

		*mc_id = atoi(argv[3]);
		if(*mc_id < 0)
		{
			Programm_Stop("Monte Carlo ID must be positive.");
			return 0;
		}

		*imax = atoi(argv[4]);
		if(*imax < 0)
		{
			Programm_Stop("imax must be a positive integer.");
			return 0;
		}

		*jmax = atoi(argv[5]);
		if(*jmax < 0)
		{
			Programm_Stop("jmax must be a positive integer.");
			return 0;
		}
	}
	else
	{
		Programm_Stop("Please pass the correct number and type of arguments (Re, MC_ID, imax, jmax).");
		return 0;
	}

	return 1;
}
