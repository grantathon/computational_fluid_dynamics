#ifndef _COLLISION_H_
#define _COLLISION_H_

#include "computeCellValues.h"
/*#include <stdlib.h>*/

#define FLUID 0
#define NO_SLIP 1
#define MOVING_WALL 2

/** computes the post-collision distribution functions according to the BGK update rule and
 *  stores the results again at the same position.
 */
void computePostCollisionDistributions(double *currentCell, const double * const tau, const double * const feq);


/** carries out the whole local collision process. Computes density and velocity and
 *  equilibrium distributions. Carries out BGK update.
 */
void doCollision(double *collideField,int *flagField,const double * const tau,int xlength);
#endif

