#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
	int imax,
	int jmax,
	double **U,
	double **V,
	int rank_l,
	int rank_r,
	int rank_b,
	int rank_t
);

#endif
