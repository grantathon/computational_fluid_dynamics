#include <stdio.h>
#include <string.h>
#include "boundary_val.h"
#include "init.h"
#include "ns_definitions.h"
#include "helper.h"
#include <mpi.h>

void boundaryvalues(int imax, int jmax, double **U, double **V, int wl, int wr, int wt, int wb, int **Flag, int rank_l, int rank_r, int rank_b, int rank_t)
{
	int i, j;

	/* For floor wall */
	if (rank_b == MPI_PROC_NULL)
	{
		boundaryvalues_b(imax, jmax, U, V, wb);
	}

	/* For ceiling */
	if (rank_t == MPI_PROC_NULL)
	{
		boundaryvalues_t(imax, jmax, U, V, wt);
	}

	/* for left wall */
	if (rank_l == MPI_PROC_NULL)
	{
		boundaryvalues_l(imax, jmax, U, V, wl);
	}

	/* for right wall */
	if (rank_r == MPI_PROC_NULL)
	{
		boundaryvalues_r(imax, jmax, U, V, wr);
	}

	/*
		loop through the Flag array and find each cell having flags B_N, B_W, B_O, B_S
	 */
	for(i = 1 ; i < imax + 1 ; i++)
	{
		for(j = 1 ; j < jmax + 1 ; j++)
		{
			/* check the boundary cells */
			if(Flag[i][j] == B_O)
			{
				U[i][j] = 0;
				V[i][j - 1] = - V[i + 1][j - 1];
				V[i][j]	= -V[i + 1][j];
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
				V[i][j - 1] = - V[i - 1][j - 1];
				V[i][j]	= -V[i - 1][j];
			}
			else if(Flag[i][j] == B_S)
			{
				V[i][j - 1] = 0;
				U[i][j] = - U[i][j - 1];
				U[i - 1][j]	= - U[i - 1][j - 1];
			}
			else if(Flag[i][j] == B_NO)
			{
				U[i][j] = 0;
				V[i][j] = 0;
				U[i - 1][j]	= -U[i - 1][j + 1];
				V[i][j - 1]	= -V[i + 1][j - 1];
			}
			else if(Flag[i][j] == B_NW)
			{
				U[i - 1][j] = 0;
				V[i][j] = 0;
				U[i][j]	= -U[i][j + 1];
				V[i][j - 1]	= -V[i - 1][j - 1];
			}
			else if(Flag[i][j] == B_SO)
			{
				U[i][j] = 0;
				V[i][j - 1] = 0;
				U[i - 1][j]	= -U[i - 1][j - 1];
				V[i][j]	= -V[i + 1][j];
			}
			else if(Flag[i][j] == B_SW)
			{
				U[i - 1][j] = 0;
				V[i][j - 1] = 0;
				U[i][j]	= -U[i][j - 1];
				V[i][j]	= -V[i - 1][j];
			}
		}
	}
}

void boundaryvalues_flag(int imax, int jmax, int il, int ir, int jb, int jt, double **flag, double **obstacle)
{

}

void boundaryvalues_t(int imax, int jmax, double **U, double **V, int wt)
{
	int i;
	/* Iterate across the ceiling */
	for(i = 1; i < imax+1; i++)
	{
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
		else if (wt == 3) /* (wb == 3) outflow BC for ceiling*/
		{
			U[i][jmax+1] 	= U[i][jmax];
			V[i][jmax] 		= V[i][jmax-1];
		}
	}
}

void boundaryvalues_b(int imax, int jmax, double **U, double **V, int wb)
{
	int i;
	/* Iterate across the floor  */
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
		else if (wb == 3)	/* outflow BC for floor*/
		{
			V[i][0]	= V[i][1];
			U[i][0] = U[i][1];
		}
	}
}

void boundaryvalues_l(int imax, int jmax, double **U, double **V, int wl)
{
	int j;
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
	}
}

void boundaryvalues_r(int imax, int jmax, double **U, double **V, int wr)
{
	int j;
	for(j = 1; j < jmax+1; j++)
	{
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
}

void spec_boundary_val(const char *problem, int imax, int jmax, int y_dim, double **U, double **V, double UI, double VI, int omg_j, int jb, int jt, int rank_l)
{
	int j;

	/*	Inlet velocity BC for upper half of inlet only; lower half cells not used by any other function etc and can be removed LATER.*/
	if (rank_l == MPI_PROC_NULL)
	{
		/*	subdomain completely inside step */
		if (jb < (jmax/2) && jt <= (jmax/2))
		{
			for(j = 1; j <= y_dim; j++)
			{
				U[0][j] = 0.0;
				V[0][j] = - V[1][j];
			}
		}
		/*	subdomain partially inside step */
		else if (jb <= (jmax/2) && jt > (jmax/2))
		{
			for(j = 1; j <= (y_dim/2 - jb); j++)
			{
				U[0][j] = 0.0;
				V[0][j] = - V[1][j];
			}

			for(j = (y_dim/2 - jb + 1); j < y_dim + 1; j++)
			{
				U[0][j] = UI;
				V[0][j] = 2 * VI - V[1][j];
			}
		}
		/*	subdomain fully out of step */
		else if (jb > (jmax/2) && jt > (jmax/2))
		{
			for(j = 1; j < y_dim; j++)
			{
				U[0][j] = UI;
				V[0][j] = 2 * VI - V[1][j];
			}
		}
	}
}

