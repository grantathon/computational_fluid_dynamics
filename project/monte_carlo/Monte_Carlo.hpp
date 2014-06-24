/*
 * MonteCarlo.hpp
 *
 *  Created on: Jun 19, 2014
 *      Author: ionut
 */

#ifndef MONTECARLO_HPP_
#define MONTECARLO_HPP_

#include <boost/random.hpp>
#include <iostream>
#include <numeric>
#include <mpi.h>

/* Monte Carlo method for UQ */
class MonteCarlo
{
private:

	/* mean and standard deviation for the Gaussian random variable */
	double mean, stddev;
	/* min and max values for the uniform distribute random variable */
	int x_min, x_max;

	/* mersenne twister random number generator */
	boost::mt19937 rng;
	/* normal(Gaussian) distribution */
	boost::normal_distribution<> normal_distr;
	/* uniform distribution */
	boost::uniform_int<> uniform_distr;

	boost::variate_generator<boost::mt19937&,
		                           boost::normal_distribution<> >* var_normal;

	boost::variate_generator<boost::mt19937&,
			                           boost::uniform_int<> >* var_uniform;

public:

	/* constructors */
	MonteCarlo(double mean, double stddev);
	MonteCarlo(int x_min, int x_max);

	/* used to generate nsamples of normal distributed random variables */
	std::vector<double> generate_nd_samples(double mean, double stddev, int nsamples);
	/* used to generate nsamples of uniformly distributed random variables */
	std::vector<double> generate_ud_samples(int x_min, int x_max, int nsamples);

	/* used to compute the first two statistical moments */
	double compute_mean(std::vector<double> v);
	double compute_variance(std::vector<double> v, double mean);

	/* used to get a normal and uniform distributed RV */
	double get_normal() const;
	double get_uniform() const;

	/* helper function for ending the mpi parallelization */
	void MCSimulation_Stop(std::string txt);
	/* domain decomposition for the second layer of parallelism */
	void init_parallel(int nsamples, int nproc, int *myrank, int *il, int *ir, int *rank_l, int *rank_r, int *omg_i, int num_proc);
	/* data decomposition - peasent style */
	void data_decomposition(int rank, int* nsampels, int* nprocs, int* samples_per_proc);

	~MonteCarlo();
};

#endif /* MONTECARLO_HPP_ */
