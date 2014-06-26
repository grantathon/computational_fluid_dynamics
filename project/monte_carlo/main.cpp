//============================================================================
// Name        : Monte.cpp
// Author      : Ionut
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "Monte_Carlo.hpp"

int main(int argc, char *argv[]) 
{
	/* MPI variables */
	int num_proc = 0;
	int myrank = 0;

	int nsamples = 10;
	int samples_per_proc = 0;

	/* for Gaussian distribution */
	double mean = 0.0;
	double stddev = 10.0;

	/* statistcs */
	double mean_nd = 0.0;
	double var_nd = 0.0;

	std::string MC_eos = "End of simulation";

	double start_time = 0.0, end_time = 0.0;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	if(myrank == 0)	
	{
		start_time = MPI_Wtime();
	}

	MonteCarlo *m = new MonteCarlo(mean, stddev);

	std::vector<double> nd_samples;
	std::vector<double> nd_qoi;

	nd_samples = m->generate_nd_samples(mean, stddev, nsamples);

	m->data_decomposition(&nsamples, &num_proc, &samples_per_proc);

	m->get_NS_solution(&samples_per_proc, nd_samples);
	m->get_QoI(&samples_per_proc, nd_qoi);

	if(myrank == 0)
	{	
		mean_nd = m->compute_mean(nd_qoi);
		var_nd = m->compute_variance(nd_qoi, mean_nd);

		std::cout << "The mean of the re-attachement point is: " << mean_nd << std::endl;
		std::cout << "The variance of the re-attachement point is:  " << var_nd << std::endl;
	}

	if (myrank == 0)
	{
		end_time = MPI_Wtime();
		printf("Elapsed time for MC simulation(%d samples) using %d processors is: %f seconds\n", nsamples, num_proc, end_time - start_time);
	}

	m->MCSimulation_Stop(MC_eos);

	delete m;

	return 0;
}
