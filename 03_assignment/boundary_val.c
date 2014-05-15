#include <stdio.h>
#include <string.h>
#include "boundary_val.h"
#include "init.h"
#include "ns_definitions.h"
#include "helper.h"

/*
 *
 * I extented the current version of the boundary conditions, where I have rewritten the spec_boundary_val
 * funtion and also, I extended the boundaryvalues function, in order to incorporate the boundary conditions
 * for boundary cells of the obstacle(with respect to their flag, e.g. if flag = B_N, then implement 1.4)
 *
 */
/*
 * 1 = no slip; 2 = free-slip; 3 = outlet
 */


void boundaryvalues(int imax, int jmax, double **U, double **V,  int wl, int wr, int wt, int wb, int **Flag)
{
	int i, j;

	/* Iterate across both the ceiling and floor  */
	for(i = 1; i < imax+1; i++)
	{
		/* FLOOR BC*/
		if (wb == 1)	/* no-slip BC for floor*/
		{
			U[i][0] = -U[i][1];
			V[i][0]	= 0;
		}
		else if (wb == 2)	/* free-slip BC for floor*/
		{
			V[i][0]	= 0;
			U[i][0] = U[i][1];
		}

		/* CEILING BC */
		if (wt == 1)	/* no-slip BC for ceiling*/
		{
			U[i][jmax+1] 	= -U[i][jmax];
			V[i][jmax]		= 0;
		}
		else if (wt == 2)	/* free-slip BC for ceiling*/
		{
			V[i][jmax]		= 0;
			U[i][jmax+1] 	= U[i][jmax];

		}
	}

	/* Iterate across both the walls*/
	for(j = 1; j < jmax+1; j++)
	{
		/* LEFT wall BC*/
		if (wl == 1)	/* no-slip BC for left wall*/
		{
			V[0][j]	= -V[1][j];
			U[0][j]	= 0;
		}
		else if (wl == 2)	/* free-slip BC for left wall*/
		{
			U[0][j]	= 0;
			V[0][j]	= V[1][j];
		}
		else if (wl == 3)	/* outflow BC for left wall*/
		{
			U[0][j]	= U[1][j];
			V[0][j]	= V[1][j];
		}

		/* RIGHT wall BC*/
		if (wr == 1)	/* no-slip BC for right wall*/
		{
			V[imax+1][j] 	= -V[imax][j];
			U[imax][j]		= 0;
		}
		else if (wr == 2)	/* free-slip BC for right wall*/
		{
			U[imax][j]		= 0;
			V[imax+1][j] 	= V[imax][j];
		}
		else if (wr == 3)	/* outflow BC for right wall*/
		{
			U[imax][j] 		= U[imax-1][j];
			V[imax+1][j] 	= V[imax][j];
		}
	}

	/*
		loop through the Flag array and find each cells have the flag B_N, B_W, B_O, B_S
		we must define the protocol for B_N, B_W, B_O, B_S
	 */


	for(i = 1 ; i < imax + 1 ; i++)
	{
		for(j = 1 ; j < jmax + 1 ; j++)
		{
			/* take the normal bounary cells */
			/* always start with the East ;) */
			if(Flag[i][j] & B_O)
			{
				U[i][j] = 0;
				V[i][j - 1] = - V[i + 1][j - 1];
				V[i][j]	= -V[i + 1][j];
			}
			else if(Flag[i][j] & B_N)
			{
				V[i][j] = 0;
				U[i - 1][j] = - U[i - 1][j + 1];
				U[i][j]	= -U[i][j + 1];
			}
			else if(Flag[i][j] & B_W)
			{
				U[i - 1][j] = 0;
				V[i][j - 1] = - V[i - 1][j - 1];
				V[i][j]	= -V[i - 1][j];
			}
			else if(Flag[i][j] & B_S)
			{
				V[i][j - 1] = 0;
				U[i][j] = - U[i][j - 1];
				U[i - 1][j]	= - U[i - 1][j - 1];
			}

			/* take the corner cells */
			if(Flag[i][j] & B_NO)
			{
				U[i][j] = 0;
				V[i][j] = 0;
				U[i - 1][j]	= -U[i - 1][j + 1];
				V[i][j - 1]	= -V[i + 1][j - 1];
			}
			else if(Flag[i][j] & B_NW)
			{
				U[i - 1][j] = 0;
				V[i][j] = 0;
				U[i][j]	= -U[i][j + 1];
				V[i][j - 1]	= -V[i - 1][j - 1];
			}
			else if(Flag[i][j] & B_SO)
			{
				U[i][j] = 0;
				V[i][j - 1] = 0;
				U[i - 1][j]	= -U[i - 1][j - 1];
				V[i][j]	= -V[i + 1][j];
			}
			else if(Flag[i][j] & B_SW)
			{
				U[i - 1][j] = 0;
				V[i][j - 1] = 0;
				U[i][j]	= -U[i][j - 1];
				V[i][j]	= -V[i - 1][j];
			}
		}
	}
}


void spec_boundary_val(const char *problem, int imax, int jmax, double **U, double **V)
{
	int j;

	/* geometry details	*/
	double UI = 0;
	double VI = 0;
	double delta_p	= 0;
	char *problemDataFile = malloc(strlen(problem) + 5);

	strcpy(problemDataFile, problem);
	strcat(problemDataFile, ".dat");

	if(strcmp(problem, "plane_shear_flow") == 0)
	{
		read_special_BC(problemDataFile, &UI, &VI, &delta_p);

		for(j = 1; j < jmax+1; j++)
		{
			U[0][j] = UI;
			V[0][j] = 2 * VI - V[1][j];
		}
	}
	else if (strcmp(problem, "Karman_vortex") == 0)
	{
		/*
		for the Karman vortex street, we consider u = 1 and v = 0 at the left wall
		 */
		read_special_BC(problemDataFile, &UI, &VI, &delta_p);
		for(j = 1; j < jmax+1; j++)
		{
			U[0][j] = UI;
			V[0][j] = 2 * VI - V[1][j];
		}
	}
	else if(strcmp(problem, "flow_over_a_step") == 0)
	{
		/*	Inlet velocity BC for upper half of inlet only; lower half cells not used by any other function etc and can be removed LATER.*/
		read_special_BC(problemDataFile, &UI, &VI, &delta_p);
		for(j = 1; j <= jmax/2; j++)
		{
			U[0][j] = 0.0;
			V[0][j] = - V[1][j];
		}

		for(j = jmax/2 + 1; j < jmax + 1; j++)
		{
			U[0][j] = UI;
			V[0][j] = 2 * VI - V[1][j];
		}

	}

	free(problemDataFile);
}
