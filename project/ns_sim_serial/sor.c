#include "sor.h"
#include "ns_definitions.h"

#include <omp.h>
#include <math.h>

void sor(
  double omg,
  double dx,
  double dy,
  int    imax,
  int    jmax,
  double **P,
  double **RS,
  double *res,
  int **Flag
) {
  int i,j;
  double rloc;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));

  /* new variable, that is used to count the total number of fluid cells ; used for the normalization
  of the residual */
  int fluid_cells_count = 0;

  /* SOR iteration */
  for(i = 1; i < imax + 1; i++)
  {
    for(j = 1; j < jmax + 1; j++)
    {
      P[i][j] = (1.0-omg)*P[i][j]
              + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);  
    }
  }

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
        fluid_cells_count++;
      }
    }
  }

  /* the residual is now normalized by the total number of fluid cells */
  rloc = rloc/(fluid_cells_count);
  rloc = sqrt(rloc);

  /* set residual */
  *res = rloc;

  /* set boundary values */
  /* pressure BC based on left boundary pressure value and pressure gradient */
  for(i = 0; i <= imax + 1; i++) {
    P[i][0] = P[i][1];
    P[i][jmax+1] = P[i][jmax];
  }
  for(j = 0; j <= jmax + 1; j++)
  {
	  P[0][j] = P[1][j];
	  P[imax+1][j] = P[imax][j];
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
}
