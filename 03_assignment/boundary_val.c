#include <stdio.h>
#include "boundary_val.h"
#include "init.h"

void boundaryvalues(int imax, int jmax, double **U, double **V,  int wl, int wr, int wt, int wb)
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
}


void spec_boundary_val(const char *szFileName, int imax, int jmax, double **U, double **V)
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

	read_special_BC(szFileName, &xlength, &ylength, &dx, &dy, &imax, &jmax, &wl, &wr, &wt, &wb, &u_wl_mov, &u_wr_mov, &u_wt_mov,
			&u_wb_mov, &u_wl_in, &u_wr_in, &u_wt_in, &u_wb_in, &in_prof_wl, &in_prof_wr, &in_prof_wt, &in_prof_wb,  &U_max);

	/* Floor */
	if (wb == 4)	/* For moving wall */
	{
		for(i = 1; i < imax+1; i++)
		{
			U[i][1] = 2 * u_wb_mov - U[i][0]; /* U[i][jmax+1] 	= 2.0 - U[i][jmax];*/
		}
	}
	else if (wb == 5 && in_prof_wb == 0)	/* constant inlet BC */
	{
		for(i = 1; i < imax+1; i++)
		{
			V[i][0] = u_wb_in;
		}
	}
	else if (wb == 5 && in_prof_wb == 1)	/* parabolic inlet BC */
	{
		for(i = 1; i < imax+1; i++)
		{
			x = i * dx;
			V[i][0] = (4 * U_max * x * (1 - (x/xlength))) / xlength;
		}
	}

	/* Ceiling BC */
	if (wt == 4)	/* For moving wall */
	{
		for(i = 1; i < imax+1; i++)
		{
			U[i][jmax+1] = 2 * u_wt_mov - U[i][jmax];
		}
	}
	else if (wt == 5 && in_prof_wt == 0)	/* constant inlet BC */
	{
		for(i = 1; i < imax+1; i++)
		{
			V[i][jmax+1] = u_wt_in;
		}
	}
	else if (wt == 5 && in_prof_wt == 1)	/* parabolic inlet BC */
	{
		for(i = 1; i < imax+1; i++)
		{
			x = i * dx;
			V[i][jmax+1] = (4 * U_max * x * (1 - (x/xlength))) / xlength;
		}
	}

	/* Left wall*/
	if (wl == 4)	/* For moving wall */
	{
		for(j = 1; j < jmax+1; j++)
		{
			V[1][j]	= 2 * u_wl_mov - V[0][j];
		}
	}
	else if (wl == 5 && in_prof_wl == 0)	/* constant inlet BC */
	{
		for(j = 1; j < jmax+1; j++)
		{
			U[0][j] = u_wl_in;
		}
	}
	else if (wl == 5 && in_prof_wl == 1)	/* parabolic inlet BC */
	{
		for(j = 1; j < jmax+1; j++)
		{
			y = j * dy;
			U[0][j] = (4 * U_max * y * (1 - (y/ylength))) / ylength;
		}
	}

	/* Right wall */
	if (wr == 4)	/* For moving wall */
	{
		for(j = 1; j < jmax+1; j++)
		{
			V[imax+1][j] = 2 * u_wr_mov - V[imax][j];
		}
	}
	else if (wr == 5 && in_prof_wr == 0)	/* constant inlet BC */
	{
		for(j = 1; j < jmax+1; j++)
		{
			U[imax+1][j] = u_wr_in;
		}
	}
	else if (wr == 5 && in_prof_wr == 1)	/* parabolic inlet BC */
	{
		for(j = 1; j < jmax+1; j++)
		{
			y = j * dy;
			U[imax+1][j] = (4 * U_max * y * (1 - (y/ylength))) / ylength;
		}
	}
}
