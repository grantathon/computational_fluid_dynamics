//============================================================================
// Name        : Monte.cpp
// Author      : Ionut
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "Monte_Carlo.hpp"
#include "solver.hpp"
#include <mpi.h>
#include <cstring>

int main(int argc, char *argv[]) 
{
	/* MPI variables */
	int num_proc = 0;
	int myrank = 0;

	int nsamples = 4;
	int samples_per_proc = 0;
	char buffer[10];

	std::string eos = "End of simulation";

	/* for Gaussian distribution */
	double mean = 0.0;
	double stddev = 10.0;

	/* statistcs */
	double mean_nd = 0.0;
	double var_nd = 0.0;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);	


	MonteCarlo *m1 = new MonteCarlo(mean, stddev);


	std::cout << "*************************************************************" << std::endl;

	std::vector<double> ndsamples;
	std::vector<double> nd_output;

	ndsamples = m1->generate_nd_samples(mean, stddev, nsamples);

	m1->data_decomposition(myrank, &nsamples, &num_proc, &samples_per_proc);

	for(int i = 0 ; i < samples_per_proc ; i++)
	{
		char test[30] = "mpirun -np 4 ./sim ";
		snprintf(buffer, sizeof(buffer), "%g %d", ndsamples[i], myrank*samples_per_proc + i + 1);	
		strcat(test, buffer);
		std :: cout << "test = " << test << std::endl;
		system(test);
		nd_output.push_back(1.0);
	}

	if(myrank == 0)
	{	
		mean_nd = m1->compute_mean(nd_output);
		var_nd = m1->compute_variance(nd_output, mean_nd);

		std::cout << "Proc " << myrank << ": mean_nd = " << mean_nd << std::endl;
		std::cout << "Proc " << myrank << ": var_nd = " << var_nd << std::endl;
	}


	std::cout << "MC simulation done; start the NS solver" << std::endl;

	m1->MCSimulation_Stop(eos);

	delete m1;

	return 0;
}
