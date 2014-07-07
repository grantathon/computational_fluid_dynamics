#ifndef MONTECARLO_HPP_
#define MONTECARLO_HPP_

#include <mpi.h>
#include <boost/random.hpp>
#include <iostream>
#include <numeric>
#include <cmath>
#include <cstring>
#include <fstream>

#include "helper.hpp"

/* Monte Carlo method for UQ */
 class MonteCarlo
 {
 private:
	/* the Raynolds number will be of the form mean + stddev*sigma, where
	sigma is ~N(0,1) or ~U(0,1) */
	double par1, par2; 

	/* mersenne twister random number generator */
	boost::mt19937 rng;
	/* normal(Gaussian) distribution */
	boost::normal_distribution<> normal_distr;
	/* uniform distribution */
	boost::uniform_real<> uniform_distr;

	boost::variate_generator<boost::mt19937&,
	boost::normal_distribution<> >* var_normal;

	boost::variate_generator<boost::mt19937&,
	boost::uniform_real<> >* var_uniform;

public:

	/* constructor */
	MonteCarlo(double par1, double par2, int distr_flag);

	/** Uncertainty (i.e. Random variables) related methods **/

	/* generate nsamples samples of normal distributed random variables */
	std::vector<double> generate_nd_samples(double mean_nd, double sttdev_nd, int* nsamples);
	/* generate nsamples samples of uniformly distributed random variables */
	std::vector<double> generate_ud_samples(double mean_ud, double sttdev_ud, int* nsamples);

	/* get a normal and uniform distributed RV */
	double get_normal() const;
	double get_uniform() const;

	/***********************************************************/

	/** Parallelization related methods **/

	/* data decomposition among processes*/
	void data_decomposition(int samples_per_proc, int* nsampels, int* nprocs, int *myrank, int *il, int *ir);

	/* call the NS solver for each generated sample */
	void monte_carlo_simulation(int *myrank, int* nsamples, int* samples_per_proc, int *il, int *ir
		, std::vector<double> &rv, int rv_flag, int imax, int jmax, double *mean, double* variance, int* flag_prog);


	/***********************************************************/

	/* destructor */
	~MonteCarlo();
};

#endif /* MONTECARLO_HPP_ */
