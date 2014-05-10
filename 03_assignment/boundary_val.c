#include <stdio.h>
#include "boundary_val.h"

void boundaryvalues(int imax, int jmax, int wl, int wr, int wt, int wb, double **U, double **V)
{
	int i, j;

	/* check the type of boudary conditions for each of the four walls;
	1 - no-slip
	2 - free-slip
	3 - outflow
	*/

	/* left wall, i.e. i = 0 */
	switch(wl)
	{
		/* no-slip bc */
		case 1:
			for(j = 1 ; j < jmax + 1 ; j++)
			{
				U[0][j] = 0;
				V[0][j] = -V[1][j];
			}
			break;
		/* free-slip bc */
		case 2:
			for(j = 1 ; j < jmax + 1 ; j++)
			{
				U[0][j] = 0;
				V[0][j] = V[1][j];
			}
			break;
		/* outflow bc */
		case 3:
			for(j = 1 ; j < jmax + 1 ; j++)
			{
				U[0][j] = U[1][j];
				V[0][j] = V[1][j];
			}
			break;
	}
	
	/* right wall, i.e. i = imax */
	switch(wr)
	{
		/* no-slip bc */
		case 1:
			for(j = 1 ; j < jmax + 1 ; j++)
			{
				U[imax][j] = 0;
				V[imax + 1][j] = -V[imax][j];
			}
			break;
		/* free-slip bc */
		case 2:
			for(j = 1 ; j < jmax + 1 ; j++)
			{
				U[imax][j] = 0;
				V[imax + 1][j] = V[imax][j];
			}
			break;
		/* outflow bc */
		case 3:
			for(j = 1 ; j < jmax + 1 ; j++)
			{
				U[imax][j] = U[imax - 1][j];
				V[imax + 1][j] = V[imax][j];
			}
			break;
	}

	/* top wall, i.e. j = jmax */
	switch(wt)
	{
		/* no-slip bc */
		case 1:
			for(i = 1 ; i < imax + 1 ; i++)
			{
				V[i][jmax] = 0;
				U[i][jmax + 1] = -U[i][jmax];
			}
			break;
		/* free-slip bc */
		case 2:
			for(i = 1 ; i < imax + 1 ; i++)
			{
				V[i][jmax] = 0;
				U[i][jmax + 1] = U[i][jmax];
			}
			break;
		/* outflow bc */
		case 3:
			for(i = 1 ; i < imax + 1 ; i++)
			{
				V[i][jmax] = V[i][jmax - 1];
				U[i][jmax + 1] = U[i][jmax];
			}
			break;
	}

	/* bottom wall, i.e. j = 0 */
	switch(wb)
	{
		/* no-slip bc */
		case 1:
			for(i = 1 ; i < imax + 1 ; i++)
			{
				V[i][0] = 0;
				U[i][0] = -U[i][1];
			}
			break;
		/* free-slip bc */
		case 2:
			for(i = 1 ; i < imax + 1 ; i++)
			{
				V[i][0] = 0;
				U[i][0] = U[i][1];
			}
			break;
		/* outflow bc */
		case 3:
			for(i = 1 ; i < imax + 1 ; i++)
			{
				V[i][0] = V[i][1];
				U[i][0] = U[i][1];
			}
			break;
	}
}

void spec_boundary_val(char *problem, int imax, int jmax, double **U, double **V)
{
	/*to be implemented */
}
