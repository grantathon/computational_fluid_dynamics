#include "solver.hpp"

void rosie(int* samples_per_proc, int* dumb_variable, double* Re)
{
	double mysum = 0.0;

	for(int i = 1 ; i < *samples_per_proc ; i++)
	{
		mysum += *Re;
	}

	MPI_Reduce(&mysum, &dumb_variable, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}