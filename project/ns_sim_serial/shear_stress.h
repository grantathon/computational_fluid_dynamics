#ifndef __SHEAR_STRESS_H_
#define __SHEAR_STRESS_H_

void separation_point_shear_stress(
  double *x_loc,
  double dx,
  double dy,
  int    imax,
  double **U,
  double **V
);

void separation_point_U(
  double *x_loc,
  double dx,
  double dy,
  int    imax,
  double **U,
  double **V
);


#endif
