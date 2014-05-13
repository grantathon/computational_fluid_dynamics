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

	*dx = *xlength / (double)(*imax);
	*dy = *ylength / (double)(*jmax);

	return 1;
}

void read_special_BC(const char * szFileName,
					double *xlength,           /* length of the domain x-dir.*/
                    double *ylength,           /* length of the domain y-dir.*/
                    double *dx,                /* length of a cell x-dir. */
                    double *dy,                /* length of a cell y-dir. */
                    int  *imax,                /* number of cells x-direction*/
                    int  *jmax,                /* number of cells y-direction*/
                    int *wl,			/* domain boundary conditions */
                    int *wr,			/* for right, left, top and */
                    int *wt,			/* bottom surfaces. */
                    int *wb,
                	double *u_wl_mov,
                	double *u_wr_mov,
                	double *u_wt_mov,
                	double *u_wb_mov,
                	double *u_wl_in,
                	double *u_wr_in,
                	double *u_wt_in,
                	double *u_wb_in,
                	int *in_prof_wl,
                	int *in_prof_wr,
                	int *in_prof_wt,
                	int *in_prof_wb,
                	double *U_max
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

	/* Read moving wall BCs*/
	READ_DOUBLE( szFileName, *u_wl_mov );
	READ_DOUBLE( szFileName, *u_wr_mov );
	READ_DOUBLE( szFileName, *u_wt_mov );
	READ_DOUBLE( szFileName, *u_wb_mov );

	/* Read string for inlet conditions */
	/*READ_INT( szFileName, inlets);*/
	READ_INT( szFileName, *in_prof_wl);
	READ_INT( szFileName, *in_prof_wr);
	READ_INT( szFileName, *in_prof_wt);
	READ_INT( szFileName, *in_prof_wb);

	/* Read inlet velocities */
	READ_DOUBLE( szFileName, *u_wl_in );
	READ_DOUBLE( szFileName, *u_wr_in );
	READ_DOUBLE( szFileName, *u_wt_in );
	READ_DOUBLE( szFileName, *u_wb_in );

	READ_DOUBLE( szFileName, *U_max );

	READ_DOUBLE( szFileName, *xlength );
	READ_DOUBLE( szFileName, *ylength );

	*dx = *xlength / (double)(*imax);
	*dy = *ylength / (double)(*jmax);
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
	/* TODO: Make sure in other implemented files that we don't
	 * 		treat border cells the same as obstacle cells */
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

			/* Check for neighboring fluid cells and set corresponding flags */
			if(obstacleFlag[i-1][j] == 1)  /* Left cell */
			{
				(*flag)[i][j] |= B_W;
			}
			if(obstacleFlag[i+1][j] == 1)  /* Right cell */
			{
				(*flag)[i][j] |= B_O;
			}
			if(obstacleFlag[i][j+1] == 1)  /* Top cell */
			{
				(*flag)[i][j] |= B_N;
			}
			if(obstacleFlag[i][j-1] == 1)  /* Bottom cell */
			{
				(*flag)[i][j] |= B_S;
			}
		}
	}

	free_imatrix(obstacleFlag, 0, imax+1, 0, jmax+1);
	free(problemPBMFile);
}



