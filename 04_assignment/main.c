#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"
#include "parallel.h"
#include "string.h"
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

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
int main(int argn, char** args)
{
	/* Input file with user parameters */
	const char *szFileName = "cavity_lab4.dat";
	char *szProblem = malloc(strlen(szFileName) + 5);
	int readParamError = 0;

	/* Geometry data */
	double xlength = 0;
	double ylength = 0;
	int imax = 0;
	int jmax = 0;
	double dx = 0;
	double dy = 0;

	/* Time-stepping data */
	double t = 0;
	double t_end = 0;
	double dt = 0;
	double tau = 0;
	double dt_value = 0;
	double visual_n = 1.0;
	int n = 0;

	/* Pressure iteration data */
	int itermax = 0;
	int it = 0;
	double res = 0;
	double eps = 0;
	double omg = 0;
	double alpha = 0;

	/* Problem dependent quantities */
	double Re = 0;
	double GX = 0;
	double GY = 0;
	double UI = 0;
	double VI = 0;
	double PI = 0;

	/* Resulting system quantities */
	double **U = 0;
	double **V = 0;
	double **P = 0;
	double **RS = 0;
	double **F = 0;
	double **G = 0;

	/* MPI variables */
	int iproc = 0;
	int jproc = 0;
	int myrank = 0;
	int il = 0;
	int ir = 0;
	int jb = 0;
	int jt = 0;
	int rank_l = 0;
	int rank_r = 0;
	int rank_b = 0;
	int rank_t = 0;
	int omg_i = 0;
	int omg_j = 0;
	int num_proc = 0;
//	const char *np_flag = "-np ";

//	char *np_command = malloc(strlen(np_flag) + 3);

	/* Initialize MPI and begin parallel processes */
	MPI_Init(&argn, &args);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	int x_dim, y_dim;

	strcpy(szProblem, szFileName);
	sprintf(szProblem, "%i", myrank);

	/* Perform necessary initializations in the main process */
	if(myrank == 0)
	{
		/* Read parameters from DAT file, store locally, and check for potential error */
		readParamError = read_parameters(szFileName, &Re, &UI, &VI, &PI, &GX, &GY,
				&t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau,
				&itermax, &eps, &dt_value, &iproc, &jproc);

		if(readParamError != 1)
		{
			Programm_Stop("Input parameters potentially corrupt!");
			return -1;
		}
	}

//	sprintf(np_command, "-np %d", iproc*jproc);
//	printf("np_command = %s\n", np_command);
//	printf("args[0] = %s\n", args[0]);
//	printf("args[1] = %s\n", args[1]);

	/* Broadcast parameters to the remaining processes */
	MPI_Bcast(&Re, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&UI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&VI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&PI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&GX, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&GY, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&t_end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&xlength, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ylength, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&dx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&dy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&imax, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&jmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&omg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&itermax, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&eps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&dt_value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&iproc, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&jproc, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if(myrank == 0)
	{
		printf("All simulation parameters have been broadcasted!\n\n");
	}

	/* Initialize process dependent variables */
	init_parallel(iproc, jproc, imax, jmax, &myrank, &il, &ir, &jb, &jt, &rank_l,
            &rank_r, &rank_b, &rank_t, &omg_i, &omg_j, num_proc);

	x_dim = ir - il + 1;
	y_dim = jt - jb + 1;

	/* Initialize matrices for velocity, pressure, rhs, etc. */
	init_uvp(UI, VI, PI, x_dim, y_dim, &U, &V, &P);
	F = matrix(0, x_dim+1, 0, y_dim+1);
	G = matrix(0, x_dim+1, 0, y_dim+1);
	RS = matrix(0, x_dim+1, 0, y_dim+1);

	/* Begin the time iteration process */
	while(t < t_end)
	{
		calculate_dt(Re, tau, &dt, dx, dy, x_dim, y_dim, U, V);
		boundaryvalues(x_dim, y_dim, U, V, rank_l, rank_r, rank_b, rank_t);
		calculate_fg(Re, GX, GY, alpha, dt, dx, dy, x_dim, y_dim, U, V, F, G, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t);
		calculate_rs(dt, dx, dy, x_dim, y_dim, F, G, RS);

		it = 0;
		do
		{
			sor(omg, dx, dy, P, RS, &res, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, imax, jmax);
			it++;
		}
		while(it < itermax && res > eps);

		calculate_uv(dt, dx, dy, x_dim, y_dim, U, V, F, G, P, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t);

		/* Visualize U, V, and P depending on dt_value */
		if(((t / dt_value) >= visual_n) || (t == dt))
		{
			write_vtkFile(szProblem, visual_n, xlength, ylength, x_dim, y_dim, dx, dy, U, V, P);
			visual_n++;
		}

		n++;
		t += dt;

		// output only for master rank
		if (myrank == 0)
		{
			printf("res=%f, it=%u, t=%f, dt=%f\n", res, it, t, dt);
		}
		//write_vtkFile(szProblem, visual_n, xlength, ylength, x_dim, y_dim, dx, dy, U, V, P);
	}

	/* Deallocate heap memory */
	free_matrix(U, 0, x_dim+1, 0, y_dim+1);
	free_matrix(V, 0, x_dim+1, 0, y_dim+1);
	free_matrix(P, 0, x_dim+1, 0, y_dim+1);
	free_matrix(F, 0, x_dim+1, 0, y_dim+1);
	free_matrix(G, 0, x_dim+1, 0, y_dim+1);
	free_matrix(RS, 0, x_dim+1, 0, y_dim+1);

	/* Finalize MPI */
	Programm_Stop("End of simulation.");

	return 0;
}
