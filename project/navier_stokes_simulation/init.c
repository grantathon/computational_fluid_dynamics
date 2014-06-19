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

	/* Set the boundary values for flag */
	// TODO: Put all of the boundary initializations in boundaryvalues_flag()
//	boundaryvalues_flag(imax, jmax, il, ir, jb, jt, flag, obstacleFlag);

	/* Set left and right boundary flags */
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
		/* Set right boundary flags */
		for(j = 0; j <= ydim+1; j++)
		{
			(*flag)[xdim+1][j] 	= C_B;	/* Right wall */
		}

		/* Set left boundary flags */
		for(j = 0; j <= ydim+1; j++)
		{
			if(obstacleFlag[il-1][jb+j-1] == 1)
			{
				(*flag)[0][j] 	= C_F;	/* Left wall */
			}
			else
			{
				(*flag)[0][j] 	= C_B;	/* Left wall */
			}

			/* Set flags according to neighboring fluids */
			if(obstacleFlag[il-2][jb+j-1] == 1)
			{
				(*flag)[0][j] |= B_W;
			}
			if(obstacleFlag[il][jb+j-1] == 1)
			{
				(*flag)[0][j] |= B_O;
			}
			if((jb+j) < jmax+1)
			{
				if(obstacleFlag[il-1][jb+j] == 1) // WARNING: May not exist!
				{
					(*flag)[0][j] |= B_N;
				}
			}
			if((jb+j-2) > 0)
			{
				if(obstacleFlag[il-1][jb+j-2] == 1) // WARNING: May not exist!
				{
					(*flag)[0][j] |= B_S;
				}
			}
		}
	}
	else if(il == 1 && ir != imax)
	{
		/* Set left boundary flags */
		for(j = 0; j <= ydim+1; j++)
		{
			(*flag)[0][j] 		= C_B;	/* Left wall */
		}

		/* Set right boundary flags */
		for(j = 0; j <= ydim+1; j++)
		{
			if(obstacleFlag[ir+1][jb+j-1] == 1)
			{
				(*flag)[xdim+1][j] 	= C_B;	/* Right wall */
			}
			else
			{
				(*flag)[xdim+1][j] 	= C_B;	/* Right wall */
			}

			/* Set flags according to neighboring fluids */
			if(obstacleFlag[ir][jb+j-1] == 1)
			{
				(*flag)[xdim+1][j] |= B_W;
			}
			if(obstacleFlag[ir+2][jb+j-1] == 1)
			{
				(*flag)[xdim+1][j] |= B_O;
			}
			if((jb+j) < jmax+1)
			{
				if(obstacleFlag[ir+1][jb+j] == 1) // WARNING: May not exist!
				{
					(*flag)[xdim+1][j] |= B_N;
				}
			}
			if((jb+j-2) > 0)
			{
				if(obstacleFlag[ir+1][jb+j-2] == 1 && jb != 1) // WARNING: May not exist!
				{
					(*flag)[xdim+1][j] |= B_S;
				}
			}
		}
	}
	else
	{
		/* Set left and right boundary flags */
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

			if(obstacleFlag[ir+1][jb+j-1] == 1)
			{
				(*flag)[xdim+1][j] 	= C_B;	/* Right wall */
			}
			else
			{
				(*flag)[xdim+1][j] 	= C_B;	/* Right wall */
			}

			/* Set flags according to neighboring fluids */
			if(obstacleFlag[il-2][jb+j-1] == 1)
			{
				(*flag)[0][j] |= B_W;
			}
			if(obstacleFlag[il][jb+j-1] == 1)
			{
				(*flag)[0][j] |= B_O;
			}
			if((jb+j) < jmax+1)
			{
				if(obstacleFlag[il-1][jb+j] == 1) // WARNING: May not exist!
				{
					(*flag)[0][j] |= B_N;
				}
			}
			if((jb+j-2) > 0)
			{
				if(obstacleFlag[il-1][jb+j-2] == 1 && jb != 1) // WARNING: May not exist!
				{
					(*flag)[0][j] |= B_S;
				}
			}

			/* Set flags according to neighboring fluids */
			if(obstacleFlag[ir][jb+j-1] == 1)
			{
				(*flag)[xdim+1][j] |= B_W;
			}
			if(obstacleFlag[ir+2][jb+j-1] == 1)
			{
				(*flag)[xdim+1][j] |= B_O;
			}
			if((jb+j) < jmax+1)
			{
				if(obstacleFlag[ir+1][jb+j] == 1) // WARNING: May not exist!
				{
					(*flag)[xdim+1][j] |= B_N;
				}
			}
			if((jb+j-2) > 0)
			{
				if(obstacleFlag[ir+1][jb+j-2] == 1) // WARNING: May not exist!
				{
					(*flag)[xdim+1][j] |= B_S;
				}
			}
		}
	}

	/* Set bottom and top boundary flags */
	if(jb == 1 && jt == jmax)
	{
		for(i = 0; i <= xdim+1; i++)
		{
			(*flag)[i][0] 		= C_B;	/* Floor */
			(*flag)[i][ydim+1] 	= C_B;	/* Ceiling */
		}
	}
	else if(jb != 1 && jt == jmax)
	{
		/* Set top boundary flags */
		for(i = 0; i <= xdim+1; i++)
		{
			(*flag)[i][ydim+1] 	= C_B;	/* Ceiling */
		}

		/* Set bottom boundary flags */
		for(i = 0; i <= xdim+1; i++)
		{
			if(obstacleFlag[il+i-1][jb-1] == 1)
			{
				(*flag)[i][0] 	= C_F;	/* Floor */
			}
			else
			{
				(*flag)[i][0] 	= C_B;	/* Floor */
			}

			/* Set flags according to neighboring fluids */
			if((il+i-2) > 0)
			{
				if(obstacleFlag[il+i-2][jb-1] == 1) // WARNING: May not exist!
				{
					(*flag)[i][0] |= B_W;
				}
			}
			if((il+i) < imax+1)
			{
				if(obstacleFlag[il+i][jb-1] == 1) // WARNING: May not exist!
				{
					(*flag)[i][0] |= B_O;
				}
			}
			if(obstacleFlag[il+i-1][jb] == 1)
			{
				(*flag)[i][0] |= B_N;
			}
			if(obstacleFlag[il+i-1][jb-2] == 1)
			{
				(*flag)[i][0] |= B_S;
			}
		}
	}
	else if(jb == 1 && jt != jmax)
	{
		/* Set bottom boundary flags */
		for(i = 0; i <= xdim+1; i++)
		{
			(*flag)[i][ydim+1] 		= C_B;	/* Floor */
		}

		/* Set top boundary flags */
		for(i = 0; i <= xdim+1; i++)
		{
			if(obstacleFlag[il+i-1][jt+1] == 1)
			{
				(*flag)[i][ydim+1] 	= C_B;	/* Ceiling */
			}
			else
			{
				(*flag)[i][ydim+1] 	= C_B;	/* Ceiling */
			}

			/* Set flags according to neighboring fluids */
			if((il+i-2) > 0)
			{
				if(obstacleFlag[il+i-2][jt+1] == 1) // WARNING: May not exist!
				{
					(*flag)[xdim+1][j] |= B_W;
				}
			}
			if((il+i) < imax+1)
			{
				if(obstacleFlag[il+i][jt+1] == 1) // WARNING: May not exist!
				{
					(*flag)[xdim+1][j] |= B_O;
				}
			}
			if(obstacleFlag[il+i-1][jt+2] == 1)
			{
				(*flag)[xdim+1][j] |= B_N;
			}
			if(obstacleFlag[il+i-1][jt] == 1)
			{
				(*flag)[xdim+1][j] |= B_S;
			}
		}
	}
	else
	{
		/* Set bottom and top boundary flags */
		for(i = 0; i <= xdim+1; i++)
		{
			if(obstacleFlag[il+i-1][jb-1] == 1)
			{
				(*flag)[i][0] 	= C_F;	/* Floor */
			}
			else
			{
				(*flag)[i][0] 	= C_B;	/* Floor */
			}

			if(obstacleFlag[il+i-1][jt+1] == 1)
			{
				(*flag)[i][ydim+1] 	= C_B;	/* Ceiling */
			}
			else
			{
				(*flag)[i][ydim+1] 	= C_B;	/* Ceiling */
			}

			/* Set flags according to neighboring fluids */
			if((il+i-2) > 0)
			{
				if(obstacleFlag[il+i-2][jb-1] == 1) // WARNING: May not exist!
				{
					(*flag)[i][0] |= B_W;
				}
			}
			if((il+i) < imax+1)
			{
				if(obstacleFlag[il+i][jb-1] == 1) // WARNING: May not exist!
				{
					(*flag)[i][0] |= B_O;
				}
			}
			if(obstacleFlag[il+i-1][jb] == 1)
			{
				(*flag)[i][0] |= B_N;
			}
			if(obstacleFlag[il+i-1][jb-2] == 1)
			{
				(*flag)[i][0] |= B_S;
			}

			/* Set flags according to neighboring fluids */
			if((il+i-2) > 0)
			{
				if(obstacleFlag[il+i-2][jt+1] == 1) // WARNING: May not exist!
				{
					(*flag)[xdim+1][j] |= B_W;
				}
			}
			if((il+i) < imax+1)
			{
				if(obstacleFlag[il+i][jt+1] == 1) // WARNING: May not exist!
				{
					(*flag)[xdim+1][j] |= B_O;
				}
			}
			if(obstacleFlag[il+i-1][jt+2] == 1)
			{
				(*flag)[xdim+1][j] |= B_N;
			}
			if(obstacleFlag[il+i-1][jt] == 1)
			{
				(*flag)[xdim+1][j] |= B_S;
			}
		}
	}

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


