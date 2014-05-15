#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"
#include "string.h"
#include <stdio.h>

/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */
int main(int argc, char *argv[])
{
	/* Input file with user parameters */
	/*const char *szFileName = "cavity100.dat";*/
	const char *problem = argv[1];
	char *problemDataFile = malloc(strlen(problem) + 5);
	char *problemOutput = malloc(strlen(problem) + 7);

	int readParamError = 0;

	/* Geometry data */
	double xlength 	= 0;
	double ylength 	= 0;
	int imax 		= 0;
	int jmax 		= 0;
	double dx 		= 0;
	double dy 		= 0;

	/* Domain boundary Condition variables */
	int wl = 0;
	int wr = 0;
	int wt = 0;
	int wb = 0;

	/* Time-stepping data */
	double t 		= 0;
	double t_end 	= 0;
	double dt 		= 0;
	double tau 		= 0;
	double dt_value = 0;
	int n 			= 0;

	/* Pressure iteration data */
	int itermax 	= 0;
	int it 			= 0;
	double res 		= 0;
	double eps 		= 0;
	double omg 		= 0;
	double alpha 	= 0;

	/* Problem dependent quantities */
	double Re = 0;
	double GX = 0;
	double GY = 0;
	double UI = 0;
	double VI = 0;
	double PI = 0;

	/* Resulting system quantities */
	double **U 	= 0;
	double **V 	= 0;
	double **P 	= 0;
	double **RS = 0;
	double **F 	= 0;
	double **G 	= 0;
	int **Flag 	= 0;

	/* User must enter correct problem strings */
	if(strcmp(problem, "Karman_vortex") && strcmp(problem, "plane_shear_flow") && strcmp(problem, "flow_over_a_step"))
	{
		printf("Run program using command 'sim [problem]', "
				"where problem is 'Karman_vortex', 'plane_shear_flow', or 'flow_over_a_step'\n");
		return 0;
	}
	strcpy(problemDataFile, problem);
	strcat(problemDataFile, ".dat");

	/* Append name for output vtk file */
	strcpy(problemOutput, problem);
	strcat(problemOutput, "_Output");

	/* Read parameters from DAT file, store locally, and check for potential error */
	readParamError = read_parameters(problemDataFile, &Re, &UI, &VI, &PI, &GX, &GY,
			&t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau,
			&itermax, &eps, &dt_value, &wl, &wr, &wt, &wb);
	if(readParamError != 1)
	{
		printf("ERROR: Input parameters potentially corrupt!");
		return 0;
	}

	/* Initialize matrices for velocity, pressure, rhs, etc. */
	init_uvp(UI, VI, PI, imax, jmax, &U, &V, &P);
	RS = matrix(0, imax+1, 0, jmax+1);
	F = matrix(0, imax+1, 0, jmax+1);
	G = matrix(0, imax+1, 0, jmax+1);

	init_flag(problem, imax, jmax, &Flag);

	/* Begin the time iteration process */
	while(t < t_end)
	{
		calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V);
		boundaryvalues(imax, jmax, U, V, wl, wr, wt, wb, Flag);

		spec_boundary_val(problem, imax, jmax, U, V);

		calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G, Flag);
		calculate_rs(dt, dx, dy, imax, jmax, F, G, RS);

		it = 0;
		do
		{
			sor(omg, dx, dy, imax, jmax, P, RS, &res, Flag);
			it++;
		}
		while( it < itermax && res > eps);
		printf("n=%u, res=%f, it=%u ", n, res, it);

		calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P, Flag);

		/* Visualize U, V, and P */
		write_vtkFile(problemOutput, n, xlength, ylength, imax, jmax, dx, dy, U, V, P);

		n++;
		t += dt;
		printf("t=%f, dt=%f\n", t, dt);
	}

	/* Visualize U, V, and P */
	write_vtkFile(problemOutput, n, xlength, ylength, imax, jmax, dx, dy, U, V, P);

	/* Print end value of U[imax/2][jmax/2], i.e., at center of the domain	*/
	printf("\nEnd value of U[imax/2][jmax/2]= %f \n", U[imax/2][(jmax/2]);

	/* Deallocate heap memory */
	free_matrix(U, 0, imax+1, 0, jmax+1);
	free_matrix(V, 0, imax+1, 0, jmax+1);
	free_matrix(P, 0, imax+1, 0, jmax+1);
	free_matrix(RS, 0, imax+1, 0, jmax+1);
	free_matrix(F, 0, imax+1, 0, jmax+1);
	free_matrix(G, 0, imax+1, 0, jmax+1);
	free_imatrix(Flag, 0, imax+1, 0, jmax+1);

	return -1;
}
