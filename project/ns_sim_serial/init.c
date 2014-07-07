#include <omp.h>
#include <string.h>

#include "helper.h"
#include "init.h"
#include "ns_definitions.h"

int read_parameters(const char * szFileName,
                    double *UI,                /* velocity x-direction */
                    double *VI,                /* velocity y-direction */
                    double *PI,                /* pressure */
                    double *GX,                /* gravitation x-direction */
                    double *GY,                /* gravitation y-direction */
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
                    int *wb
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

	READ_DOUBLE( szFileName, *dt    );

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

	*dx = *xlength / (double)(*imax);
	*dy = *ylength / (double)(*jmax);

	return 1;
}

void read_special_BC(const char * szFileName,
		double *UI,	/* Read inlet velocity BC*/
		double *VI,
		double *delta_p
)
{
	READ_DOUBLE( szFileName, *UI );
	READ_DOUBLE( szFileName, *VI );

	if (strcmp(szFileName, "plane_shear_flow.dat") == 0)
	{
		READ_DOUBLE( szFileName, *delta_p );
	}

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

void init_flag(const char *problem, int imax, int jmax, int ***flag)
{
	char *problemPBMFile = 0;
	int **obstacleFlag = 0;
	int i, j;

	/* initialize Flag as a matrix of integers */
	*flag = imatrix(0, imax+1, 0, jmax+1);

	/* Set boundary flags */
	for(i = 0; i <= imax+1; i++)
	{
		(*flag)[i][0] 		= C_B;	/* Floor */
		(*flag)[i][jmax+1] 	= C_B;	/* Ceiling */
	}

	for(j = 0; j <= jmax+1; j++)
	{
		(*flag)[0][j] 		= C_B;	/* Left wall */
		(*flag)[imax+1][j] 	= C_B;	/* Right wall */
	}

	/* Retrieve full file name */
	problemPBMFile = malloc(strlen(problem) + 5);
	strcpy(problemPBMFile, problem);
	strcat(problemPBMFile, ".pbm");

	/* Input file containing fluid/obstacle data */
	obstacleFlag = read_pgm(problemPBMFile);

	/* Set inner flags */
	for(i = 1; i <= imax; i++)
	{
		for(j = 1; j <= jmax; j++)
		{
			/* Determine whether cell is a fluid or an obstacle */
			if(obstacleFlag[i][j] == 1)
			{
				(*flag)[i][j] = C_F;  /* Fluid */
			}
			else
			{
				(*flag)[i][j] = C_B;  /* Obstacle/border */
			}
		}
	}

	for(i = 1; i <= imax; i++)
	{
		for(j = 1; j <= jmax; j++)
		{
			/* Set flags according to neighbouring fluids */
			if(obstacleFlag[i-1][j] == 1)
			{
				(*flag)[i][j] |= B_W;
			}
			if(obstacleFlag[i+1][j] == 1)
			{
				(*flag)[i][j] |= B_O;
			}
			if(obstacleFlag[i][j+1] == 1)
			{
				(*flag)[i][j] |= B_N;
			}
			if(obstacleFlag[i][j-1] == 1)
			{
				(*flag)[i][j] |= B_S;
			}
		}
	}

	/* Specify the left and right boundaries where pressure BC is defined	*/
	/* We only specify pressure value on the left boundary	*/
	if (strcmp(problem, "plane_shear_flow") == 0)
	{
		for(j = 1; j <= jmax; j++)
		{
			(*flag)[0][j] = P_L;

		}
	}

	free_imatrix(obstacleFlag, 0, imax+1, 0, jmax+1);
	free(problemPBMFile);
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
				printf("Reynolds number must be greater than zero.");
				return 0;
			}
		}
		else
		{
			*viscosity = atof(argv[2]);
			if(*viscosity <= 0)
			{
				printf("Viscosity number must be greater than zero.");
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
			printf("Monte Carlo ID must be positive.");
			return 0;
		}

		*imax = atoi(argv[4]);
		if(*imax < 0)
		{
			printf("imax must be a positive integer.");
			return 0;
		}

		*jmax = atoi(argv[5]);
		if(*jmax < 0)
		{
			printf("jmax must be a positive integer.");
			return 0;
		}
	}
	else
	{
		printf("Please pass the correct number and type of arguments (Re, MC_ID, imax, jmax).");
		return 0;
	}

	return 1;
}

