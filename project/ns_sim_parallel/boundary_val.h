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
		int **Flag,
		int rank_l,
		int rank_r,
		int rank_b,
		int rank_t
);

void boundaryvalues_flag(int imax, int jmax, int il, int ir, int jb, int jt, int ***flag, int **obstacle);

void boundaryvalues_t(int imax, int jmax, double **U, double **V, int wt);

void boundaryvalues_b(int imax, int jmax, double **U, double **V, int wb);

void boundaryvalues_l(int imax, int jmax, double **U, double **V, int wl);

void boundaryvalues_r(int imax, int jmax, double **U, double **V, int wr);


/* Additional BCs implemented here.*/
void spec_boundary_val(
		int imax,
		int jmax,
		int y_dim,
		double **U,
		double **V,
		double UI,
		double VI,
		int omg_j,
		int jb,
		int jt,
		int rank_l
);


#endif
