#include <stdio.h>
#include "boundary_val.h"

void boundaryvalues(int imax, int jmax, double **U, double **V)
{
	int i, j;

	if(imax == jmax)
	{
		/* Iterate across the walls, ceiling, and floor */
		for(i = 1; i < imax+1; i++)
		{
			U[i][0] 		= -U[i][1];
			U[i][imax+1]	= 2.0-U[i][imax];

			U[0][i] 	= 0;
			U[imax][i] 	= 0;

			V[i][0] 	= 0;
			V[i][imax] 	= 0;

			V[0][i] 		= -V[1][i];
			V[imax+1][i] 	= -V[imax][i];
		}
	}
	else
	{
		/* Iterate across both walls */
		for(i = 1; i < imax+1; i++)
		{
			U[i][0] 		= -U[i][1];
			U[i][jmax+1] 	= 2.0-U[i][jmax];

			V[i][0] 	= 0;
			V[i][jmax] 	= 0;
		}

		/* Iterate across both the ceiling and floor */
		for(j = 1; j < jmax+1; j++)
		{
			U[0][j] 	= 0;
			U[imax][j] 	= 0;

			V[0][j]			= -V[1][j];
			V[imax+1][j]	= -V[imax][j];
		}
	}
}
