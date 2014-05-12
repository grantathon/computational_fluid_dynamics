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
		int wb
);

/* Additional BCs implemented here.*/
void spec_boundary_val(
		const char *problem,
		int imax,
		int jmax,
		double **U,
		double **V
);


#endif
