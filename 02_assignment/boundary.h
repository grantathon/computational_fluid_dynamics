#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#define FLUID 0
#define NO_SLIP 1
#define MOVING_WALL 2

/** handles the boundaries in our simulation setup */
void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity,int xlength);

#endif

