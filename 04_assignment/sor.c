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
  int rank_t
)
{
  int i,j;
  double rloc;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));
  double *bufSend = 0;
  double *bufRecv = 0;
  MPI_Status status;
  int chunk = 0;

  /* SOR iteration */
  for(i = il; i <= ir; i++) {
    for(j = jb; j <= jt; j++) {
      P[i][j] = (1.0-omg)*P[i][j]
              + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }

  /* Communicate between processes regarding pressure boundaries */
  pressure_comm(P, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, &status, chunk);

  /* compute the residual */
  rloc = 0;
  for(i = il; i <= ir; i++) {
    for(j = jb; j <= jt; j++) {
      rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
              ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }

  rloc = rloc/(ir*jt);
  rloc = sqrt(rloc);
  *res = rloc;	/* set residual */

  /* set boundary values */
  /* TODO: Do we treat process boundaries the same as domain boundaries? */
  for(i = il; i <= ir; i++) {
    P[i][0] = P[i][1];
    P[i][jt+1] = P[i][jt];
  }
  for(j = jb; j <= jt; j++) {
    P[0][j] = P[1][j];
    P[ir+1][j] = P[ir][j];
  }
}
