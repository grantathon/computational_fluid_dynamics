/*
 * Monte_Carlo.cpp
 *
 *  Created on: Jun 19, 2014
 *      Author: ionut
 */

#include "Monte_Carlo.hpp"

MonteCarlo::MonteCarlo(double mean, double stddev)
{
	rng = boost::mt19937();

	normal_distr = boost::normal_distribution<>(mean, stddev);

	var_normal = new boost::variate_generator<boost::mt19937&,
			boost::normal_distribution<> >(rng, normal_distr);
}

MonteCarlo::MonteCarlo(int x_min, int x_max)
{
	rng = boost::mt19937();

	uniform_distr = boost::uniform_int<>(x_min, x_max);

	var_uniform = new boost::variate_generator<boost::mt19937&,
			boost::uniform_int<> >(rng, uniform_distr);
}

double MonteCarlo::get_normal() const {
	return (*var_normal)();
}

double MonteCarlo::get_uniform() const {
	return (*var_uniform)();
}

std::vector<double> MonteCarlo::generate_nd_samples(double mean,
		double stddev, int nsamples) {
	std::vector<double> ndsamples;
	double temp;

	for (int i = 0; i < nsamples; i++) {
		temp = get_normal();
		ndsamples.push_back(temp);
	}

	return ndsamples;
}

std::vector<double> MonteCarlo::generate_ud_samples(int x_min, int x_max,
		int nsamples) {
	std::vector<double> udsamples;
	double temp;

	for (int i = 0; i < nsamples; i++) {
		temp = get_uniform();
		udsamples.push_back(temp);
	}

	return udsamples;
}

double MonteCarlo::compute_mean(std::vector<double> v) {
	double sum = 0.0, mean = 0.0;

	sum = std::accumulate(v.begin(), v.end(), 0.0);
	mean = sum / v.size();

	return mean;
}

double MonteCarlo::compute_variance(std::vector<double> v, double mean) {
	double square_sum = 0.0, stddev = 0.0;

	std::vector<double> diff(v.size());
	std::transform(v.begin(), v.end(), diff.begin(),
			std::bind2nd(std::minus<double>(), mean));
	square_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(),
			0.0);
	stddev = square_sum / (v.size() - 1);

	return stddev;
}

void MonteCarlo::MCSimulation_Stop(std::string txt)
/* all processes will produce a text output, be synchronized and finished */
{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                           /* synchronize output */
   std::cout << "-STOP- P:" << myrank << " " << txt << std::endl;
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(1);
}

void MonteCarlo::init_parallel(int nsamples, int nproc, int *myrank, int *il, int *ir, int *rank_l, int *rank_r, int *omg_i, int num_proc)
{
	int i_per_iproc, i_rem;

	/* Set sub-domain coordinates */
	*omg_i = ((*myrank) % nproc) + 1;
	
	/* Compute il, ir for each sub-domain*/
	i_per_iproc = (nsamples / nproc);
	i_rem = (nsamples % nproc);

	/* for left and right*/
	if ((*omg_i) == 1)   /* to rank zero assign the remainder*/
	{
	  *il = (((*omg_i) -1) * i_per_iproc) + 1;
	}
	else
	{
	  *il = (((*omg_i) -1) * i_per_iproc) + i_rem + 1;
	}
	*ir = ((*omg_i) * i_per_iproc) + i_rem;

	/* Assign rank of neighbour to each sub-domain: rank_l, rank_r*/
	/* Left boundary*/
	if ((*il) == 1)
	{
	  *rank_l = MPI_PROC_NULL;
	}
	else
	{
	  *rank_l = (*myrank) - 1;
	}

	/* Right boundary*/
	if ((*ir) == nsamples)
	{
	  *rank_r = MPI_PROC_NULL;
	}
	else
	{
	  *rank_r = (*myrank) + 1;
	}

	std::cout <<"rank= "<< *myrank << std::endl;
	std::cout << "omg_i= " <<*omg_i << std::endl;
	std::cout << "il= " << *il << " " << "ir= " << *ir << " " << "rank_l= " <<*rank_l << " " << "rank_r= " << *rank_r << std::endl;
}

MonteCarlo::~MonteCarlo() {
	delete var_normal;
	delete var_uniform;
}

