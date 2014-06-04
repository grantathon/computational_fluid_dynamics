#include "sor.h"
#include "parallel.h"
#include <math.h>
#include <mpi.h>

void sor(
  double omg,
  double dx,
  double dy,
  double **P,
  double **RS,
  double *res,
  int il,
  int ir,
  int jb,
  int jt,
  int rank_l,
  int rank_r,
  int rank_b,
  int rank_t,
  int imax,
  int jmax
)
{
	int i, j;
	double rloc = 0.0;
	double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));
	double *bufSend = 0;
	double *bufRecv = 0;
	MPI_Status status;
	int chunk = 0;
	int x_dim = ir - il + 1;
	int y_dim = jt - jb + 1;

	/* Set left & right global domain boundaries according to Neumann boundary conditions */
	if(rank_l == MPI_PROC_NULL && rank_r != MPI_PROC_NULL)  /* Only receive/send data from/to right */
	{
		for(j = 1; j <= y_dim; j++) {
			P[0][j] = P[1][j];
		}
	}
	else if(rank_l != MPI_PROC_NULL && rank_r == MPI_PROC_NULL)  /* Only send/receive data to/from left */
	{
		for(j = 1; j <= y_dim; j++) {
			P[x_dim+1][j] = P[x_dim][j];
		}
	}
	else if(rank_l == MPI_PROC_NULL && rank_r == MPI_PROC_NULL)  /* No bordering processes */
	{
		for(j = 1; j <= y_dim; j++) {
			P[0][j] = P[1][j];
			P[x_dim+1][j] = P[x_dim][j];
		}
	}

	/* Set top & bottom global domain boundaries according to Neumann boundary conditions */
	if(rank_t == MPI_PROC_NULL && rank_b != MPI_PROC_NULL)  /* Only receive/send data from/to bottom */
	{
		/* Implement Neumann conditions on global domain boundaries */
		for(i = 1; i <= x_dim; i++) {
			P[i][y_dim+1] = P[i][y_dim];
		}
	}
	else if(rank_t != MPI_PROC_NULL && rank_b == MPI_PROC_NULL)  /* Only send/receive data to/from top */
	{
		/* Implement Neumann conditions on global domain boundaries */
		for(i = 1; i <= x_dim; i++) {
			P[i][0] = P[i][1];
		}
	}
	else if(rank_t == MPI_PROC_NULL && rank_b == MPI_PROC_NULL)  /* No bordering processes */
	{
		/* Implement Neumann conditions on global domain boundaries */
		for(i = 1; i <= x_dim; i++) {
			P[i][0] = P[i][1];
			P[i][y_dim+1] = P[i][y_dim];
		}
	}

	/* SOR iteration */
	for(i = 1; i <= (ir - il + 1); i++) {
		for(j = 1; j <= (jt - jb + 1); j++) {
			P[i][j] = (1.0-omg)*P[i][j]
				  + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
		}
	}

	/* Communicate between processes regarding pressure boundaries */
	pressure_comm(P, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, &status, chunk);

	/* compute the residual */
	for(i = 1; i <= (ir - il + 1); i++) {
		for(j = 1; j <= (jt - jb + 1); j++) {
			rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
				  ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
		}
	}

	/* Sum the squares of all local residuals then square root that sum for global residual */
	MPI_Allreduce(&rloc, res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	*res = sqrt((*res)/(imax*jmax));

}
