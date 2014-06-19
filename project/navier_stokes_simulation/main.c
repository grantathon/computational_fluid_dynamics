#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"
#include "string.h"
#include "parallel.h"

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
int main(int argc, char *argv[])
{
	/* Input file with user parameters */
	const char *problem = "flow_over_a_step";
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
	double start_time = 0.0;
	double end_time = 0.0;

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
	double *output_values = 0;
	double **U 	= 0;
	double **V 	= 0;
	double **P 	= 0;
	double **RS = 0;
	double **F 	= 0;
	double **G 	= 0;
	int **flag 	= 0;

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
	int comm_size = 0;
	int x_dim, y_dim;

	/* Initialize MPI and begin parallel processes */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	/* Setup retrieval of data configuration file */
	strcpy(problemDataFile, problem);
	strcat(problemDataFile, ".dat");

	/* Append name for output vtk file */
	strcpy(problemOutput, problem);
	strcat(problemOutput, "_Output");

	/* Perform necessary initializations in the main process */
	if(myrank == 0)
	{
		/* Start timer */
		start_time = MPI_Wtime();

		/* Read parameters from DAT file, store locally, and check for potential error */
		readParamError = read_parameters(problemDataFile, &Re, &UI, &VI, &PI, &GX, &GY,
				&t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau,
				&itermax, &eps, &dt_value, &wl, &wr, &wt, &wb, &iproc, &jproc);

		if(readParamError != 1)
		{
			Programm_Stop("Input parameters potentially corrupt!");
			return 0;
		}

		MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
		if(comm_size != (iproc*jproc))
		{
			printf("The number of processes entered via 'mpirun -np %u' does not equal the size specified by the input file (iproc * jproc = %u)\n", comm_size, iproc*jproc);
			return 0;
		}
	}

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

	/* Init flags indicating cell type and neighbors cell types */
	init_flag(problem, x_dim, y_dim, imax, jmax, il, ir, jb, jt, &flag);

	/* Begin the time iteration process */
	printf("Begin the main computation...\n");
	while(t < t_end)
	{
		boundaryvalues(imax, jmax, U, V, wl, wr, wt, wb, flag, rank_l, rank_r, rank_b, rank_t);
		//boundaryvalues(imax, jmax, U, V, wl, wr, wt, wb, flag);

		spec_boundary_val(problem, imax, jmax, U, V, UI, VI, rank_l);
		//spec_boundary_val(problem, imax, jmax, U, V, UI, VI);

		calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G, flag);
		//calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G, flag);
		calculate_rs(dt, dx, dy, imax, jmax, F, G, RS);

		it = 0;
		do
		{
			sor(omg, dx, dy, P, RS, &res, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, imax, jmax, flag);
			//sor(omg, dx, dy, imax, jmax, P, RS, &res, flag, Pw, delta_p);
			it++;
		}
		while( it < itermax && res > eps);

		// Output info to user
		printf("n=%u, res=%f, it=%u, t=%f, dt=%f\n", n, res, it, t, dt);

		calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P, flag, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t);
		//calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P, flag);
		calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V);
		//calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V);
//        write_vtkFile(problemOutput, n, xlength, ylength, imax, jmax, dx, dy, U, V, P);

		n++;
		t += dt;
	}

	/* Visualize U, V, and P */
	write_vtkFile(problemOutput, n, xlength, ylength, imax, jmax, dx, dy, U, V, P);

	/* TODO: Write output array to shared file */
	output_values = malloc(4*sizeof(double));
	output_values[0] = 0;  	// Re
	output_values[1] = 0;	// U
	output_values[2] = 0;	// V
	output_values[3] = 0;	// P

	/* Deallocate heap memory */
	free_matrix(U, 0, x_dim+1, 0, y_dim+1);
	free_matrix(V, 0, x_dim+1, 0, y_dim+1);
	free_matrix(P, 0, x_dim+1, 0, y_dim+1);
	free_matrix(F, 0, x_dim+1, 0, y_dim+1);
	free_matrix(G, 0, x_dim+1, 0, y_dim+1);
	free_matrix(RS, 0, x_dim+1, 0, y_dim+1);
	free_imatrix(flag, 0, x_dim+1, 0, y_dim+1);

	if(myrank == 0)
	{
		/* End timer */
		end_time = MPI_Wtime();
		printf("Program time: %f\n", (end_time - start_time));
	}

	/* Finalize MPI */
	Programm_Stop("End of simulation.");

	return 1;
}
