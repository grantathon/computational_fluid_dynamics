#include "helper.h"
#include "init.h"
#include "ns_definitions.h"
#include <string.h>

int read_parameters(const char * szFileName,
					double *Re,                /* reynolds number   */
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

	READ_DOUBLE( szFileName, *Re    );
	READ_DOUBLE( szFileName, *t_end );
	READ_DOUBLE( szFileName, *dt    );

	READ_INT   ( szFileName, *imax );
	READ_INT   ( szFileName, *jmax );

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

	/* Set boundary flags */
	if(il == 1 && ir == imax)
	{
		for(j = 0; j <= ydim+1; j++)
		{
			(*flag)[0][j] 		= C_B;	/* Left wall */
			(*flag)[xdim+1][j] 	= C_B;	/* Right wall */
		}
	}
	else if(il != 1 && ir == imax)
	{
		for(j = 0; j <= ydim+1; j++)
		{
			(*flag)[xdim+1][j] 	= C_B;	/* Right wall */
		}

		for(j = 0; j <= ydim+1; j++)
		{
			if(obstacleFlag[il-1][jb+j-1] == 1)
			{
				(*flag)[0][j] 		= C_F;	/* Left wall */
			}
			else
			{
				(*flag)[0][j] 		= C_B;	/* Left wall */
			}
		}

//		for(j = 1; j <= jmax; j++)
		for(j = 0; j <= ydim+1; j++)
		{
			/* Set flags according to neighboring fluids */
			if(obstacleFlag[il-2][jb+j-1] == 1)
			{
				(*flag)[0][j] |= B_W;
			}
			if(obstacleFlag[il][jb+j-1] == 1)
			{
				(*flag)[0][j] |= B_O;
			}
			if(obstacleFlag[il-1][jb+j] == 1)
			{
				(*flag)[0][j] |= B_N;
			}
			if(obstacleFlag[il-1][jb+j-2] == 1)
			{
				(*flag)[0][j] |= B_S;
			}
		}
	}
	else if(il == 1 && ir != imax)
	{
		for(j = 0; j <= ydim+1; j++)
		{
			(*flag)[0][j] 		= C_B;	/* Left wall */
		}
	}
	else
	{

	}

//	for(i = 0; i <= imax+1; i++)
//	{
//		(*flag)[i][0] 		= C_B;	/* Floor */
//		(*flag)[i][jmax+1] 	= C_B;	/* Ceiling */
//	}
//	for(j = 0; j <= jmax+1; j++)
//	{
//		(*flag)[0][j] 		= C_B;	/* Left wall */
//		(*flag)[imax+1][j] 	= C_B;	/* Right wall */
//	}


	/* Set inner flags */
//	for(i = 1; i <= imax; i++)
	for(i = il; i <= ir; i++)
	{
//		for(j = 1; j <= jmax; j++)
		for(j = jb; j <= jt; j++)
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

//	for(i = 1; i <= imax; i++)
	for(i = il; i <= ir; i++)
	{
//		for(j = 1; j <= jmax; j++)
		for(j = jb; j <= jt; j++)
		{
			/* Set flags according to neighboring fluids */
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

	free_imatrix(obstacleFlag, 0, xdim+1, 0, ydim+1);
	free(problemPBMFile);
}


