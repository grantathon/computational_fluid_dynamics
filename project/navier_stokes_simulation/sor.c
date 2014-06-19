#include "sor.h"
#include "ns_definitions.h"
#include "parallel.h"
#include <math.h>
#include <mpi.h>

void sor(double omg, double dx, double dy, double **P, double **RS, double *res, int il, int ir, int jb, int jt, int rank_l, int rank_r, int rank_b, int rank_t, int imax, int jmax, int **Flag)
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

  /* new local variable, that is used to count the total number of fluid cells ; used for the normalization
  of the residual */
  int local_fluid_cells_count = 0;
  int fluid_cells_count = 0;

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
    for(i = 1; i <= x_dim; i++) {
      P[i][y_dim+1] = P[i][y_dim];
    }
  }
  else if(rank_t != MPI_PROC_NULL && rank_b == MPI_PROC_NULL)  /* Only send/receive data to/from top */
  {
    for(i = 1; i <= x_dim; i++) {
      P[i][0] = P[i][1];
    }
  }
  else if(rank_t == MPI_PROC_NULL && rank_b == MPI_PROC_NULL)  /* No bordering processes */
  {
    for(i = 1; i <= x_dim; i++) {
      P[i][0] = P[i][1];
      P[i][y_dim+1] = P[i][y_dim];
    }
  }

  /* boundary conditions for the pressure at the boundary stripe*/
  for(i = 1 ; i < imax + 1 ; i++)
  {
      for(j = 1 ; j < jmax + 1 ; j++)
      {
        if(Flag[i][j] == B_O)
        {
          P[i][j] = P[i + 1][j];
        } 
        else if(Flag[i][j] == B_N)
        {
          P[i][j] = P[i][j + 1];    
        }
        else if(Flag[i][j] == B_W)
        {
          P[i][j] = P[i - 1][j];    
        }
        else if(Flag[i][j] == B_S)
        {
          P[i][j] = P[i][j - 1];       
        }
        else if(Flag[i][j] == B_NO)
        {
          P[i][j] = (P[i][j + 1] + P[i + 1][j])/2;
        }
        else if(Flag[i][j] == B_NW)
        {
          P[i][j] = (P[i][j + 1] + P[i - 1][j])/2;
        }
        else if(Flag[i][j] == B_SO)
        {
          P[i][j] = (P[i][j - 1] + P[i + 1][j])/2;
        } 
        else if(Flag[i][j] == B_SW)
        {
          P[i][j] = (P[i][j - 1] + P[i - 1][j])/2;    
        }
      }
  }

  /* SOR iteration */
  for(i = 1; i < imax + 1; i++)
  {
    for(j = 1; j < jmax + 1; j++)
    {
      P[i][j] = (1.0-omg)*P[i][j]
              + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);  
    }
  }

  /* Communicate between processes regarding pressure boundaries */
  pressure_comm(P, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, &status, chunk);

  /* compute the residual */
  /* now, the computation is restricted only for the fluid cells */
  rloc = 0;
  for(i = 0; i <= imax; i++) 
  {
    for(j = 0; j <= jmax; j++) 
    {
      if(Flag[i][j] & C_F)
      {
        rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
                ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
        local_fluid_cells_count++;
      }
    }
  }

  /* Sum the squares of all local residuals then square root that sum for global residual */
  MPI_Allreduce(&rloc, res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  /* Sum all local fluid cell counts*/
  MPI_Allreduce(&local_fluid_cells_count, &fluid_cells_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  *res = sqrt((*res)/fluid_cells_count);
}
