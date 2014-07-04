#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"
#include "shear_stress.h"
#include "string.h"

#include <stdio.h>

int main(int argc, char *argv[])
{
	/* Input file with user parameters */
	const char *problem = "flow_over_a_step";
	char *problemDataFile = malloc(strlen(problem) + 4);
	char *problemOutput = malloc(strlen(problem) + 10);
	char *simOutput = malloc(strlen("sim_") + 20);
	char buffer[10] = {0};
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
	double dt 		= 0;
	double tau 		= 0;
	double dt_value = 0;
	int n 			= 0;
//	double visual_n = 1;

	/* Pressure iteration data */
	int itermax 	= 0;
	int it 			= 0;
	double res 		= 1;
	double eps 		= 0;
	double omg 		= 0;
	double alpha 	= 0;

	/* Problem dependent quantities */
	double viscosity = 0;
	double rho = 1.225; /* standard air density in units of kg/m^3*/
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
	int flag_Re = 1;	/*	flag for Re or viscosity input, default value (1) is set for Re input*/

	/*	Parameters for flow reattachment and output */
	double x_loc = 0;
	int mc_id = 0;

	/* Read passed input parameters */
	if(read_args(argc, argv, &flag_Re, &Re, &viscosity, &mc_id, &imax, &jmax, &rho) == 0)
	{
		return 0;
	}

	/* Setup retrieval of data configuration file */
	strcpy(problemDataFile, problem);
	strcat(problemDataFile, ".dat");

	/* Append name for output vtk file */
	strcpy(problemOutput, problem);
	strcat(problemOutput, "_output");

	/* Setup simulation output file */
	strcpy(simOutput, "ns_sim_");
	sprintf(buffer, "%i", mc_id);
	strcat(simOutput, buffer);
	strcat(simOutput, ".mc");

	/* Read parameters from DAT file, store locally, and check for potential error */
	readParamError = read_parameters(problemDataFile, &UI, &VI, &PI, &GX, &GY,
			&xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau,
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

	/* Init flags indicating cell type and neighbours cell types */
	init_flag(problem, imax, jmax, &Flag);

	/* Begin the time iteration process */
	printf("Begin the simulation for id = %u\n", mc_id);
	while(res > eps)
	{
		calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V);
		boundaryvalues(imax, jmax, U, V, wl, wr, wt, wb, Flag);

		spec_boundary_val(imax, jmax, U, V, UI, VI);

		calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G, Flag);
		calculate_rs(dt, dx, dy, imax, jmax, F, G, RS);

		it = 0;
		do
		{
			sor(omg, dx, dy, imax, jmax, P, RS, &res, Flag);
			it++;
		}
		while(it < itermax && res > eps);
//		printf("n=%u, res=%f, it=%u, t=%f, dt=%f\n", n, res, it, t, dt);

		calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P, Flag);
//		if((t / dt_value) >= visual_n)
//		{
//			write_vtkFile(problemOutput, visual_n, xlength, ylength, imax, jmax, dx, dy, U, V, P, Flag);
//			visual_n++;
//		}

		n++;
		t += dt;
	}

	/* Visualize U, V, and P */
//	write_vtkFile(problemOutput, visual_n, xlength, ylength, imax, jmax, dx, dy, U, V, P, Flag);

	separation_point_U(&x_loc, dx, dy, imax, U, V);
//	printf("re-attachment point: %f \t time: %f\n", x_loc, t);

	/* Write simulation output values */

	/* Reynolds number */
	if(write_to_file((const char*)simOutput, Re) == 0)
	{
		return 0;
	}

	/* Re-attachment point location*/
	if(write_to_file((const char*)simOutput, x_loc) == 0)
	{
		return 0;
	}

	/* Re-attachment point time (steady state)*/
	if(write_to_file((const char*)simOutput, t) == 0)
	{
		return 0;
	}

	/* Deallocate heap memory */
	free_matrix(U, 0, imax+1, 0, jmax+1);
	free_matrix(V, 0, imax+1, 0, jmax+1);
	free_matrix(P, 0, imax+1, 0, jmax+1);
	free_matrix(RS, 0, imax+1, 0, jmax+1);
	free_matrix(F, 0, imax+1, 0, jmax+1);
	free_matrix(G, 0, imax+1, 0, jmax+1);
	free_imatrix(Flag, 0, imax+1, 0, jmax+1);

	return 1;
}
