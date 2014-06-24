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
	/*int il = 0;
	int ir = 0;
	int rank_l = 0;
	int rank_r = 0;
	int omg_i = 0; */
	int num_proc = 0;

	int myrank = 0;
	int nsamples = 0;

	std::string eos = "End of simulation";
	char test[30] = "./sim_ns ";
	//char buffer[10];

	//const char* launch_program = "mpirun -np 4 ./sim_ns 100";

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

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

	nsamples = 10;

	//int dumb_variable_nd;
	int dumb_variable_ud;
	int samples_per_proc;


	MonteCarlo *m1 = new MonteCarlo(mean, stddev);
	MonteCarlo *m2 = new MonteCarlo(x_min, x_max);

	/*m1->init_parallel(nsamples, nproc, &myrank, &il, &ir, &rank_l, &rank_r, &omg_i, num_proc);
	m2->init_parallel(nsamples, nproc, &myrank, &il, &ir, &rank_l, &rank_r, &omg_i, num_proc);*/


	std::cout << "*************************************************************" << std::endl;

	std::vector<double> ndsamples, udsamples;
	std::vector<double> nd_output, ud_output;

	ndsamples = m1->generate_nd_samples(mean, stddev, nsamples);
	udsamples = m2->generate_ud_samples(x_min, x_max, nsamples);

	m1->data_decomposition(myrank, &nsamples, &num_proc, &samples_per_proc);

	for(int i = 0 ; i < samples_per_proc ; i++)
	{
		//snprintf(buffer, sizeof(buffer), "%g", ndsamples[i]);
		//strcat(test, buffer);
		//rosie(&samples_per_proc, &dumb_variable_nd, &ndsamples[i]);
		system(test);
		nd_output.push_back(1.0);
	}

	for(int i = 0 ; i < nsamples ; i++)
	{
		rosie(&samples_per_proc, &dumb_variable_ud, &udsamples[i]);
		ud_output.push_back(dumb_variable_ud);
	}

	if(myrank == 0)
	{	
		//std::cout << "something nd = " << " " << dumb_variable_nd << std::endl;
		std::cout << "something ud = " << " " << dumb_variable_ud << std::endl;

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


	std::cout << "MC simulation done; start the NS solver :>" << std::endl;
	//system(test);

	m1->MCSimulation_Stop(eos);
	m2->MCSimulation_Stop(eos);

	delete m1;
	delete m2;

	return 0;
}
