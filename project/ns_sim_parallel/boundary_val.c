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

void boundaryvalues_flag(int imax, int jmax, int il, int ir, int jb, int jt, int ***flag, int **obstacle)
{
	int i, j;
	int xdim = ir - il + 1;
	int ydim = jt - jb + 1;

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
			if(obstacle[il-1][jb+j-1] == 1)
			{
				(*flag)[0][j] 	= C_F;	/* Left wall */
			}
			else
			{
				(*flag)[0][j] 	= C_B;	/* Left wall */
			}

			/* Set flags according to neighboring fluids */
			if(obstacle[il-2][jb+j-1] == 1)
			{
				(*flag)[0][j] |= B_W;
			}
			if(obstacle[il][jb+j-1] == 1)
			{
				(*flag)[0][j] |= B_O;
			}
			if((jb+j) < jmax+1)
			{
				if(obstacle[il-1][jb+j] == 1) // WARNING: May not exist!
				{
					(*flag)[0][j] |= B_N;
				}
			}
			if((jb+j-2) > 0)
			{
				if(obstacle[il-1][jb+j-2] == 1) // WARNING: May not exist!
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
			if(obstacle[ir+1][jb+j-1] == 1)
			{
				(*flag)[xdim+1][j] 	= C_B;	/* Right wall */
			}
			else
			{
				(*flag)[xdim+1][j] 	= C_B;	/* Right wall */
			}

			/* Set flags according to neighboring fluids */
			if(obstacle[ir][jb+j-1] == 1)
			{
				(*flag)[xdim+1][j] |= B_W;
			}
			if(obstacle[ir+2][jb+j-1] == 1)
			{
				(*flag)[xdim+1][j] |= B_O;
			}
			if((jb+j) < jmax+1)
			{
				if(obstacle[ir+1][jb+j] == 1) // WARNING: May not exist!
				{
					(*flag)[xdim+1][j] |= B_N;
				}
			}
			if((jb+j-2) > 0)
			{
				if(obstacle[ir+1][jb+j-2] == 1 && jb != 1) // WARNING: May not exist!
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
			if(obstacle[il-1][jb+j-1] == 1)
			{
				(*flag)[0][j] 		= C_F;	/* Left wall */
			}
			else
			{
				(*flag)[0][j] 		= C_B;	/* Left wall */
			}

			if(obstacle[ir+1][jb+j-1] == 1)
			{
				(*flag)[xdim+1][j] 	= C_B;	/* Right wall */
			}
			else
			{
				(*flag)[xdim+1][j] 	= C_B;	/* Right wall */
			}

			/* Set flags according to neighboring fluids */
			if(obstacle[il-2][jb+j-1] == 1)
			{
				(*flag)[0][j] |= B_W;
			}
			if(obstacle[il][jb+j-1] == 1)
			{
				(*flag)[0][j] |= B_O;
			}
			if((jb+j) < jmax+1)
			{
				if(obstacle[il-1][jb+j] == 1) // WARNING: May not exist!
				{
					(*flag)[0][j] |= B_N;
				}
			}
			if((jb+j-2) > 0)
			{
				if(obstacle[il-1][jb+j-2] == 1 && jb != 1) // WARNING: May not exist!
				{
					(*flag)[0][j] |= B_S;
				}
			}

			/* Set flags according to neighboring fluids */
			if(obstacle[ir][jb+j-1] == 1)
			{
				(*flag)[xdim+1][j] |= B_W;
			}
			if(obstacle[ir+2][jb+j-1] == 1)
			{
				(*flag)[xdim+1][j] |= B_O;
			}
			if((jb+j) < jmax+1)
			{
				if(obstacle[ir+1][jb+j] == 1) // WARNING: May not exist!
				{
					(*flag)[xdim+1][j] |= B_N;
				}
			}
			if((jb+j-2) > 0)
			{
				if(obstacle[ir+1][jb+j-2] == 1) // WARNING: May not exist!
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
			if(obstacle[il+i-1][jb-1] == 1)
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
				if(obstacle[il+i-2][jb-1] == 1) // WARNING: May not exist!
				{
					(*flag)[i][0] |= B_W;
				}
			}
			if((il+i) < imax+1)
			{
				if(obstacle[il+i][jb-1] == 1) // WARNING: May not exist!
				{
					(*flag)[i][0] |= B_O;
				}
			}
			if(obstacle[il+i-1][jb] == 1)
			{
				(*flag)[i][0] |= B_N;
			}
			if(obstacle[il+i-1][jb-2] == 1)
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
			if(obstacle[il+i-1][jt+1] == 1)
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
				if(obstacle[il+i-2][jt+1] == 1) // WARNING: May not exist!
				{
					(*flag)[xdim+1][j] |= B_W;
				}
			}
			if((il+i) < imax+1)
			{
				if(obstacle[il+i][jt+1] == 1) // WARNING: May not exist!
				{
					(*flag)[xdim+1][j] |= B_O;
				}
			}
			if(obstacle[il+i-1][jt+2] == 1)
			{
				(*flag)[xdim+1][j] |= B_N;
			}
			if(obstacle[il+i-1][jt] == 1)
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
			if(obstacle[il+i-1][jb-1] == 1)
			{
				(*flag)[i][0] 	= C_F;	/* Floor */
			}
			else
			{
				(*flag)[i][0] 	= C_B;	/* Floor */
			}

			if(obstacle[il+i-1][jt+1] == 1)
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
				if(obstacle[il+i-2][jb-1] == 1) // WARNING: May not exist!
				{
					(*flag)[i][0] |= B_W;
				}
			}
			if((il+i) < imax+1)
			{
				if(obstacle[il+i][jb-1] == 1) // WARNING: May not exist!
				{
					(*flag)[i][0] |= B_O;
				}
			}
			if(obstacle[il+i-1][jb] == 1)
			{
				(*flag)[i][0] |= B_N;
			}
			if(obstacle[il+i-1][jb-2] == 1)
			{
				(*flag)[i][0] |= B_S;
			}

			/* Set flags according to neighboring fluids */
			if((il+i-2) > 0)
			{
				if(obstacle[il+i-2][jt+1] == 1) // WARNING: May not exist!
				{
					(*flag)[xdim+1][j] |= B_W;
				}
			}
			if((il+i) < imax+1)
			{
				if(obstacle[il+i][jt+1] == 1) // WARNING: May not exist!
				{
					(*flag)[xdim+1][j] |= B_O;
				}
			}
			if(obstacle[il+i-1][jt+2] == 1)
			{
				(*flag)[xdim+1][j] |= B_N;
			}
			if(obstacle[il+i-1][jt] == 1)
			{
				(*flag)[xdim+1][j] |= B_S;
			}
		}
	}
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

void spec_boundary_val(int imax, int jmax, int y_dim, double **U, double **V, double UI, double VI, int omg_j, int jb, int jt, int rank_l)
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

