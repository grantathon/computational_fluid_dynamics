#include <stdio.h>
#include "boundary_val.h"
#include "init.h"
#include "ns_definitions.h"

/*
*
* I extented the current version of the boundary conditions, where I have rewritten the spec_boundary_val
* funtion and also, I extended the boundaryvalues function, in order to incorporate the boundary conditions
* for boundary cells of the obstacle(with respect to their flag, e.g. if flag = B_N, then implement 1.4)
*
*/


void boundaryvalues(int imax, int jmax, double **U, double **V,  int wl, int wr, int wt, int wb, int **Flag)
{
	int i, j;

	/* Iterate across both the ceiling and floor  */
	for(i = 1; i < imax+1; i++)
	{
		/* FLOOR BC*/
		if (wb == 1 || wb == 4 || wb == 5)	/* no-slip BC for floor*/
		{
			U[i][0] = -U[i][1];
			V[i][0]	= 0;
		}
		else if (wb == 2)	/* free-slip BC for floor*/
		{
			V[i][0]	= 0;
			U[i][0] = U[i][1];
		}
		else if (wb == 3) /*outflow BC for floor*/
		{
			U[i][0] = U[i][1];
			V[i][0]	= V[i][1];
		}

		/* CEILING BC */
		if (wt == 1 || wt == 4 || wt == 5)	/* no-slip BC for ceiling*/
		{
			U[i][jmax+1] 	= -U[i][jmax];
			V[i][jmax]		= 0;
		}
		else if (wt == 2)	/* free-slip BC for ceiling*/
		{
			V[i][jmax]		= 0;
			U[i][jmax+1] 	= U[i][jmax];

		}
		else if (wt == 3) /* (wb == 3) outflow BC for ceiling*/
		{
			U[i][jmax+1] 	= U[i][jmax];
			V[i][jmax] 		= V[i][jmax-1];
		}
	}

	/* Iterate across both the walls*/
	for(j = 1; j < jmax+1; j++)
	{
		/* LEFT wall BC*/
		if (wl == 1 || wl == 4 || wl == 5)	/* no-slip BC for left wall*/
		{
			V[0][j]	= -V[1][j];
			U[0][j]	= 0;
		}
		else if (wl == 2)	/* free-slip BC for left wall*/
		{
			U[0][j]	= 0;
			V[0][j]	= V[1][j];
		}
		else if (wl == 3) /* (wb == 3) outflow BC for left wall*/
		{
			U[0][j]	= U[1][j];
			V[0][j]	= V[1][j];
		}

		/* RIGHT wall BC*/
		if (wr == 1 || wr == 4 || wr == 5)	/* no-slip BC for right wall*/
		{
			V[imax+1][j] 	= -V[imax][j];
			U[imax][j]		= 0;
		}
		else if (wr == 2)	/* free-slip BC for right wall*/
		{
			U[imax][j]		= 0;
			V[imax+1][j] 	= V[imax][j];
		}
		else if (wr == 3) /* (wb == 3) outflow BC for right wall*/
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
			for(j = 1 ; j < imax + 1 ; j++)
			{
				/* take the normal bounary cells */
				/* always start with the East ;) */
				if(Flag[i][j] == B_O)
				{
					U[i][j] = 0;
					V[i][j - 1] = - V[i][j - 1];
					V[i][j]	= -V[i - 1][j];				
				} 
				else if(Flag[i][j] == B_N)
				{
					V[i][j] = 0;
					U[i - 1][j] = - U[i - 1][j + 1];
					U[i][j]	= -U[i][j + 1];			
				}
				else if(Flag[i][j] == B_W)
				{
					U[i - 1][j] = 0;
					V[i][j - 1] = - V[i][j - 1];
					V[i][j]	= -V[i - 1][j];			
				}
				else if(Flag[i][j] == B_S)
				{
					V[i][j-1] = 0;
					U[i][j] = - U[i][j + 1];
					U[i + 1][j]	= - U[i + 1][j + 1];				
				}

				/* take the corner cells */ 
				if(Flag[i][j] == B_NO)
				{
								
				}
				else if(Flag[i][j] == B_NW)
				{
								
				}
				else if(Flag[i][j] == B_SO)
				{
								
				} 
				else if(Flag[i][j] == B_SW)
				{
								
				}	
			}
		}
}


void spec_boundary_val(char *problem, int imax, int jmax, double **U, double **V)
{
	int i, j;	/* loop indices */

	/* Domain wall BC type */
	int wl = 0;
	int wr = 0;
	int wt = 0;
	int wb = 0;

	/* geometry details	*/
	double xlength 	= 0;
	double ylength 	= 0;
	double U_max 	= 0;
	double dx 		= 0;
	double dy 		= 0;
	double x 		= 0;	/* actual x-coordinate on the grid of a cell */
	double y 		= 0;

	/*	MOving wall velocities for each face of the domain	*/
	double u_wl_mov = 0.0;
	double u_wr_mov = 0.0;
	double u_wt_mov = 0.0;
	double u_wb_mov = 0.0;

	/* Velocities for the inlet BC; one for each face*/
	double u_wl_in = 0.0;
	double u_wr_in = 0.0;
	double u_wt_in = 0.0;
	double u_wb_in = 0.0;

	/* For inlet velocity BC this tells the type of inlet profile */
	int in_prof_wl = 0;
	int in_prof_wr = 0;
	int in_prof_wt = 0;
	int in_prof_wb = 0;

	/*read_special_BC(szFileName, &xlength, &ylength, &dx, &dy, &imax, &jmax, &wl, &wr, &wt, &wb, &u_wl_mov, &u_wr_mov, &u_wt_mov,
			&u_wb_mov, &u_wl_in, &u_wr_in, &u_wt_in, &u_wb_in, &in_prof_wl, &in_prof_wr, &in_prof_wt, &in_prof_wb,  &U_max);

	/* Floor */
	
	/* For moving wall */
	/*
	if (wb == 4)	
	{
		for(i = 1; i < imax+1; i++)
		{
			U[i][1] = 2 * u_wb_mov - U[i][0]; 
		}
	}
	else if (wb == 5 && in_prof_wb == 0)
	{
		for(i = 1; i < imax+1; i++)
		{
			V[i][0] = u_wb_in;
		}
	}
	else if (wb == 5 && in_prof_wb == 1)	
	{
		for(i = 1; i < imax+1; i++)
		{
			x = i * dx;
			V[i][0] = (4 * U_max * x * (1 - (x/xlength))) / xlength;
		}
	}
	*/

	/* Ceiling BC */
	/*
	if (wt == 4)	
	{
		for(i = 1; i < imax+1; i++)
		{
			U[i][jmax+1] = 2 * u_wt_mov - U[i][jmax];
		}
	}
	else if (wt == 5 && in_prof_wt == 0)	
	{
		for(i = 1; i < imax+1; i++)
		{
			V[i][jmax+1] = u_wt_in;
		}
	}
	else if (wt == 5 && in_prof_wt == 1)	
	{
		for(i = 1; i < imax+1; i++)
		{
			x = i * dx;
			V[i][jmax+1] = (4 * U_max * x * (1 - (x/xlength))) / xlength;
		}
	}
	*/

	/* Left wall*/
	/*
	if (wl == 4)	
	{
		for(j = 1; j < jmax+1; j++)
		{
			V[1][j]	= 2 * u_wl_mov - V[0][j];
		}
	}
	else if (wl == 5 && in_prof_wl == 0)	
	{
		for(j = 1; j < jmax+1; j++)
		{
			U[0][j] = u_wl_in;
		}
	}
	else if (wl == 5 && in_prof_wl == 1)	
	{
		for(j = 1; j < jmax+1; j++)
		{
			y = j * dy;
			U[0][j] = (4 * U_max * y * (1 - (y/ylength))) / ylength;
		}
	}
	*/

	/* Right wall */
	/*
	if (wr == 4)	
	{
		for(j = 1; j < jmax+1; j++)
		{
			V[imax+1][j] = 2 * u_wr_mov - V[imax][j];
		}
	}
	else if (wr == 5 && in_prof_wr == 0)	
	{
		for(j = 1; j < jmax+1; j++)
		{
			U[imax+1][j] = u_wr_in;
		}
	}
	else if (wr == 5 && in_prof_wr == 1)	
	{
		for(j = 1; j < jmax+1; j++)
		{
			y = j * dy;
			U[imax+1][j] = (4 * U_max * y * (1 - (y/ylength))) / ylength;
		}
	}*/

	if(strcmp(problem, "plane_shear_flow") == 0)
	{
		/*
		ask Ayman for this; practically, it is a parabolic profile of the form 4Umax*y(1-y)
		*/
		for(j = 1; j < jmax+1; j++)
		{
			y = j * dy;
			U[0][j] = (4 * U_max * y * (1 - (y/ylength))) / ylength;
		}
	}
	else if (strcmp(problem, "Karman_vortex") == 0)
	{
		/*
		for the Karman vortex street, we consider u = 1 and v = 0 at the left wall
		*/
		for(j = 1; j < jmax+1; j++)
		{
			U[0][j] = 1.0;
			V[0][j]	= 0.0;
		}
	}
	else if(strcmp(problem, "flow_over_a_step") == 0)
	{
		/*
		if I am not wrong, the upper half is u = 1, v = 0 and the lower half is u = 0 and 
		v = 0, because of the step; they say that IT MIGHT BE DONE IN init_uvp, but...:-j
		*/	
		for(j = 1; j <= jmax/2; j++)
		{
			U[0][j] = 0.0;
			V[0][j]	= 0.0;
		}

		for(j = jmax/2 + 1; j < jmax + 1; j++)
		{
			U[0][j] = 1.0;
			V[0][j]	= 0.0;
		}

	}
}
