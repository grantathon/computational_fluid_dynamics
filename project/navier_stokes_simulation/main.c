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

int main(int argc, char *argv[])
{
	/* Input file with user parameters */
	const char *problem = "flow_over_a_step";
	char *problemDataFile = malloc(strlen(problem) + 4);
	char *problemOutput = malloc(strlen(problem) + 10);
	char buffer[4] = {0};
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
	double start_time = 0;
	double end_time = 0;
	double visual_n = 1;

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
	strcat(problemOutput, "_output.");
	sprintf(buffer, "%i", myrank);
	strcat(problemOutput, buffer);

	/* Setup input file var for reading */
	strcpy(problemDataFile, problem);
	strcat(problemDataFile, ".dat");

	/* Read passed Reynolds number */
	if(argc > 1)
	{
		Re = atof(argv[1]);
		if(Re <= 0)
		{
			Programm_Stop("Reynolds number must be greater than zero!");
			return 0;
		}
	}

	/* Perform necessary initializations in the main process */
	if(myrank == 0)
	{
		/* Start timer */
		start_time = MPI_Wtime();

		/* Read parameters from DAT file, store locally, and check for potential error */
		readParamError = read_parameters(problemDataFile, &UI, &VI, &PI, &GX, &GY,
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
	broadcast_parameters(myrank, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt,
			&dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value, &wl, &wr,
			&wt, &wb, &iproc, &jproc);

	/* Initialize process dependent variables */
	init_parallel(iproc, jproc, imax, jmax, &myrank, &il, &ir, &jb, &jt, &rank_l,
            &rank_r, &rank_b, &rank_t, &omg_i, &omg_j, num_proc);

	/* Set dimensions of local process */
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
		boundaryvalues(x_dim, y_dim, U, V, wl, wr, wt, wb, flag, rank_l, rank_r, rank_b, rank_t);
		spec_boundary_val(imax, jmax, y_dim, U, V, UI, VI, omg_j, jb, jt, rank_l);

		calculate_fg(Re, GX, GY, alpha, dt, dx, dy, x_dim, y_dim, U, V, F, G, flag);
		calculate_rs(dt, dx, dy, x_dim, y_dim, F, G, RS);

		it = 0;
		do
		{
			sor(omg, dx, dy, P, RS, &res, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, imax, jmax, flag);
			it++;
		}
		while(it < itermax && res > eps);

		calculate_uv(dt, dx, dy, x_dim, y_dim, U, V, F, G, P, flag, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t);

		/* Visualize U, V, and P depending on dt_value */
		if((t / dt_value) >= visual_n)
		{
			output_uvp(U, V, P, flag, il, ir, jb, jt, omg_i, omg_j, problemOutput, visual_n);
			visual_n++;
		}

		// output sim stats to user by master rank
		if (myrank == 0)
		{
			printf("res=%f, it=%u, t=%f, dt=%f\n", res, it, t, dt);
		}

		calculate_dt(Re, tau, &dt, dx, dy, x_dim, y_dim, U, V);

		t += dt;
		n++;
	}

	/* Visualize last output of U, V, and P */
	output_uvp(U, V, P, flag, il, ir, jb, jt, omg_i, omg_j, problemOutput, visual_n);

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
