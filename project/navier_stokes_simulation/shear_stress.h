#ifndef __SHEAR_STRESS_H_
#define __SHEAR_STRESS_H_

void shear_stress_calc(
  double *x_loc,
  double dx,
  double dy,
  int il,
  int    imax,
  double **U,
  double **V
);


#endif
