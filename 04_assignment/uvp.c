#include "math.h"
#include "uvp.h"
#include "parallel.h"
#include <mpi.h>

double matrix_abs_max(double **A, int imax, int jmax)
{
	int i, j;
	double max_val = 0;

	/* Iterate only over the inner cells */
	for(i = 1; i < imax+1; i++)
	{
		for(j = 1; j < jmax+1; j++)
		{
			if(fabs(A[i][j]) > max_val)
			{
				/* Only grab the maximum absolute value */
				max_val = fabs(A[i][j]);
			}
		}
	}

	return max_val;
}

void calculate_fg(
  double Re,
  double GX,
  double GY,
  double alpha,
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  int il,
  int ir,
  int jb,
  int jt,
  int rank_l,
  int rank_r,
  int rank_b,
  int rank_t
)
{
	int i, j;
	double du_x_2, du_y_2, dv_x_2, dv_y_2, duv_x, duv_y, du2_x, dv2_y;

	/* Iterate only over the inner cells of F */
	for(i = 1; i <= imax-1; i++)
	{
		for(j = 1; j <= jmax; j++)
		{
			/* Setup 1st and 2nd order derivatives for U and V */
			du_x_2 = (U[i-1][j] - 2*U[i][j] + U[i+1][j])/pow(dx, 2);
			du_y_2 = (U[i][j-1] - 2*U[i][j] + U[i][j+1])/pow(dy, 2);

			duv_y = (1/dy)*((((V[i][j] + V[i+1][j])*(U[i][j] + U[i][j+1]))/4)
						- (((V[i][j-1] + V[i+1][j-1])*(U[i][j-1] + U[i][j]))/4))
					+ (alpha/dy)*(((fabs(V[i][j] + V[i+1][j])*(U[i][j] - U[i][j+1]))/4)
							- ((fabs(V[i][j-1] + V[i+1][j-1])*(U[i][j-1] - U[i][j]))/4));

			du2_x = (1/dx)*(pow((U[i][j] + U[i+1][j])/2, 2) - pow((U[i-1][j] + U[i][j])/2, 2))
					+ (alpha/dx)*(((fabs(U[i][j] + U[i+1][j])*(U[i][j] - U[i+1][j]))/4)
							- ((fabs(U[i-1][j] + U[i][j])*(U[i-1][j] - U[i][j]))/4));

			F[i][j] = U[i][j] + dt*(((du_x_2 + du_y_2)/Re) - du2_x - duv_y + GX);
		}
	}

	/* Iterate only over the inner cells of G */
	for(i = 1; i <= imax; i++)
	{
		for(j = 1; j <= jmax-1; j++)
		{
			/* Setup 1st and 2nd order derivatives for U and V */
			dv_x_2 = (V[i-1][j] - 2*V[i][j] + V[i+1][j])/pow(dx, 2);
			dv_y_2 = (V[i][j-1] - 2*V[i][j] + V[i][j+1])/pow(dy, 2);

			duv_x = (1/dx)*((((U[i][j] + U[i][j+1])*(V[i][j] + V[i+1][j]))/4)
						- (((U[i-1][j] + U[i-1][j+1])*(V[i-1][j] + V[i][j]))/4))
					+ (alpha/dx)*(((fabs(U[i][j] + U[i][j+1])*(V[i][j] - V[i+1][j]))/4)
							- ((fabs(U[i-1][j] + U[i-1][j+1])*(V[i-1][j] - V[i][j]))/4));

			dv2_y = (1/dy)*(pow((V[i][j] + V[i][j+1])/2, 2) - pow((V[i][j-1] + V[i][j])/2, 2))
					+ (alpha/dy)*(((fabs(V[i][j] + V[i][j+1])*(V[i][j] - V[i][j+1]))/4)
							- ((fabs(V[i][j-1] + V[i][j])*(V[i][j-1] - V[i][j]))/4));

			G[i][j] = V[i][j] + dt*(((dv_x_2 + dv_y_2)/Re) - duv_x - dv2_y + GY);
		}
	}

	/* Set boundary values */
	for(i = 1; i < imax+1; i++)
	{
		G[i][0] = V[i][0];
		G[i][jmax] = V[i][jmax];
	}

	for(j = 1; j < jmax+1; j++)
	{
		F[0][j] = U[0][j];
		F[imax][j] = U[imax][j];
	}
}


/**
 * This operation computes the right hand side of the pressure poisson equation.
 * The right hand side is computed according to the formula
 *
 * @f$ rs = \frac{1}{\delta t} \left( \frac{F^{(n)}_{i,j}-F^{(n)}_{i-1,j}}{\delta x} + \frac{G^{(n)}_{i,j}-G^{(n)}_{i,j-1}}{\delta y} \right)  @f$
 *
 */
void calculate_rs(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **F,
  double **G,
  double **RS
)
{
	int i, j;

	/* Iterate only over the inner cells */
	for(i = 1; i <= imax; i++)
	{
		for(j = 1; j <= jmax; j++)
		{
			RS[i][j] = ( (F[i][j] - F[i-1][j])/dx + (G[i][j] - G[i][j-1])/dy)/dt;
		}
	}
}


/**
 * Determines the maximal time step size. The time step size is restricted
 * accordin to the CFL theorem. So the final time step size formula is given
 * by
 *
 * @f$ {\delta t} := \tau \, \min\left( \frac{Re}{2}\left(\frac{1}{{\delta x}^2} + \frac{1}{{\delta y}^2}\right)^{-1},  \frac{{\delta x}}{|u_{max}|},\frac{{\delta y}}{|v_{max}|} \right) @f$
 *
 */
void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V
)
{
	/* Only update dt if tau is positive */
	if(tau > 0)
	{
		double loc_uv_max[2] = {0.0};
		double glob_uv_max[2] = {0.0};

		/* Find the max of both velocity matrices for the local domain */
		loc_uv_max[0] = matrix_abs_max(U, imax, jmax);
		loc_uv_max[1] = matrix_abs_max(V, imax, jmax);

		/* Communicate with the other processes to find the global max values of U and V */
		MPI_Allreduce(loc_uv_max, glob_uv_max, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

		*dt = tau*fmin((Re/2)*(1/((1/pow(dx, 2)) + (1/pow(dy, 2)))), fmin(dx/glob_uv_max[0], dy/glob_uv_max[1]));
	}
}


/**
 * Calculates the new velocity values according to the formula
 *
 * @f$ u_{i,j}^{(n+1)}  =  F_{i,j}^{(n)} - \frac{\delta t}{\delta x} (p_{i+1,j}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 * @f$ v_{i,j}^{(n+1)}  =  G_{i,j}^{(n)} - \frac{\delta t}{\delta y} (p_{i,j+1}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 *
 * As always the index range is
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 * @image html calculate_uv.jpg
 */
void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P,
  int il,
  int ir,
  int jb,
  int jt,
  int rank_l,
  int rank_r,
  int rank_b,
  int rank_t
)
{
	double *bufSend = 0;
	double *bufRecv = 0;
	MPI_Status status;
	int chunk = 0;
	int i, j;

	/* Iterate only over the inner cells */
	/* for U_ij */
	for(i = 1; i <= imax-1; i++)
	{
		for(j = 1; j <= jmax; j++)
		{
			U[i][j] = F[i][j] - (dt/dx) * (P[i+1][j] - P[i][j]);
		}
	}

	/* for V_ij */
	for(i = 1; i <= imax; i++)
	{
		for(j = 1; j <= jmax-1; j++)
		{
			V[i][j] = G[i][j] - (dt/dy) * (P[i][j+1] - P[i][j]);
		}
	}

	/* Communicate between processes regarding velocities */
	uv_comm(U, V, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, &status, chunk);
}
