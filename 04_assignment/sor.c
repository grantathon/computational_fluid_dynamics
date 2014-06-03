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

	/* set boundary values */
	for(i = 1; i <= (ir - il + 1); i++) {
		P[i][0] = P[i][1];
		P[i][jt - jb + 2] = P[i][jt - jb + 1];
	}
	for(j = 1; j <= (jt - jb + 1); j++) {
		P[0][j] = P[1][j];
		P[ir - il + 2][j] = P[ir - il + 1][j];
	}
}
