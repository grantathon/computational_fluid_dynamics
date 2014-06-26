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
#include <cmath>
#include <cstring>
#include <fstream>
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

	/** Uncertainty (i.e. Random variables) related methods **/

	/* generate nsamples samples of normal distributed random variables */
	std::vector<double> generate_nd_samples(double mean, double stddev, int nsamples);
	/* generate nsamples samples of uniformly distributed random variables */
	std::vector<double> generate_ud_samples(int x_min, int x_max, int nsamples);

	/* compute the first two statistical moments */
	double compute_mean(const std::vector<double> &v) const;
	double compute_variance(const std::vector<double> &v, double mean) const;

	/* get a normal and uniform distributed RV */
	double get_normal() const;
	double get_uniform() const;

	/***********************************************************/

	/** Parallelization related methods **/

	/* data decomposition among processes*/
	void data_decomposition(int* nsampels, int* nprocs, int* samples_per_proc);

	/* call the NS solver for each generated sample */
	void get_NS_solution(int* samples_per_proc, const std::vector<double> &Re);

	/* get the QoI (Quantities of interest - the desired output parameters, from a UQ point of view) */
	void get_QoI(int* samples_per_proc, std::vector<double> &qoi);

	/* helper function for ending the mpi parallelization */
	void MCSimulation_Stop(const std::string &txt);

	/***********************************************************/

	/* destructor */
	~MonteCarlo();
};

#endif /* MONTECARLO_HPP_ */
