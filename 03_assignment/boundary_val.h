#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/**
 * The boundary values of the problem are set(in the case of no-slip, free-slip and outflow bc).
 */
void boundaryvalues(
	int imax, 
	int jmax, 
	int wl, 
	int wr, 
	int wt, 
	int wb, 
	double **U, 
	double **V);

/**
 * The boundary values of the problem are set(in the case of inflow bc).
 */
void spec_boundary_val(
	char *problem,
	int imax,
	int jmax,
	double **U,
	double **V);

#endif
