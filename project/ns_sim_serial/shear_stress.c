#include "shear_stress.h"
#include "ns_definitions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void shear_stress_calc(
  double *x_loc,
  double dx,
  double dy,
  int    imax,
  double **U,
  double **V
)
{
  int i;
  double dv_dx, du_dy, stress, stress_ahead;

  du_dy = (U[imax][2] - U[imax][0])/(2*dy);
  dv_dx = (V[imax + 1][1] - V[imax - 1][1])/(2*dx);
  stress_ahead = du_dy + dv_dx;

  for(i = imax-1; i > 1; i--)
  {
    du_dy = (U[i][2] - U[i][0])/(2*dy);
    dv_dx = (V[i + 1][1] - V[i - 1][1])/(2*dx);
    stress = du_dy + dv_dx;

    if (stress_ahead > 0 && stress < 0)
    {
      *x_loc = i * dx;
      break;
    }

    stress_ahead = stress;

    printf("%f \t %f \t %f \t %f \n ", (i * dx),  stress, du_dy, dv_dx);
  }
}