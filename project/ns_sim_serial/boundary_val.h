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
		int wl,
		int wr,
		int wt,
		int wb,
		int **Flag
);

/* Additional BCs implemented here.*/
void spec_boundary_val(
		int imax,
		int jmax,
		double **U,
		double **V,
		double UI,
		double VI
);


#endif
