//============================================================================
// Name        : Monte.cpp
// Author      : Ionut
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "Monte_Carlo.hpp"
#include <mpi.h>

double dumb_function(double Re)
{
	double a = 10.0, b = 5.0;

	return a*Re + b;
}

int main(int argc, char *argv[]) 
{
	/* MPI variables */
	int nproc = 0;
	int myrank = 0;
	int il = 0;
	int ir = 0;
	int rank_l = 0;
	int rank_r = 0;
	int omg_i = 0;
	int num_proc = 0;
	int nsamples = 0;

	std::string eos = "End of simulation";

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	double temp;

	/* for Gaussian distribution */
	double mean = 0.0;
	double stddev = 10.0;

	/* for uniform distribution */
	int x_min = 0;
	int x_max = 10;

	/* statistcs */
	double mean_nd = 0.0;
	double var_nd = 0.0;

	double mean_ud = 0.0;
	double var_ud = 0.0;

	nsamples = 100;
	nproc = 4;


	MonteCarlo *m1 = new MonteCarlo(mean, stddev);
	MonteCarlo *m2 = new MonteCarlo(x_min, x_max);

	m1->init_parallel(nsamples, nproc, &myrank, &il, &ir, &rank_l, &rank_r, &omg_i, num_proc);
	m2->init_parallel(nsamples, nproc, &myrank, &il, &ir, &rank_l, &rank_r, &omg_i, num_proc);


	std::cout << "*************************" << std::endl;

	std::vector<double> ndsamples, udsamples;
	std::vector<double> nd_output, ud_output;

	ndsamples = m1->generate_nd_samples(mean, stddev, nsamples);
	udsamples = m2->generate_ud_samples(x_min, x_max, nsamples);

	for(unsigned int i = 0 ; i < ndsamples.size() ; i++)
	{
		temp = dumb_function(ndsamples[i]);
		nd_output.push_back(temp);
	}

	for(unsigned int i = 0 ; i < udsamples.size() ; i++)
	{
		temp = dumb_function(udsamples[i]);
		ud_output.push_back(temp);
	}

	if(myrank == 0)
	{	
		mean_nd = m1->compute_mean(nd_output);
		var_nd = m1->compute_variance(nd_output, mean_nd);

		mean_ud = m2->compute_mean(ud_output);
		var_ud = m2->compute_variance(ud_output, mean_ud);

		std::cout << "Proc " << myrank << ": mean_nd = " << mean_nd << std::endl;
		std::cout << "Proc " << myrank << ": var_nd = " << var_nd << std::endl;

		std::cout<< "**********************************************" << std::endl;

		std::cout << "Proc " << myrank << ": mean_ud = " << mean_ud << std::endl;
		std::cout << "Proc " << myrank << ": var_ud = " << var_ud << std::endl;
	}

	delete m1;
	delete m2;

	m1->MCSimulation_Stop(eos);
	m2->MCSimulation_Stop(eos);

	return 0;
}
