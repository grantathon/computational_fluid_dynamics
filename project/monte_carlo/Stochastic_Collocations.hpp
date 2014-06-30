/*
 * StochasticCollocations.hpp
 *
 *  Created on: Jun 27, 2014
 *      Author: ionut
 */

#ifndef STOCHASTICCOLLOCATIONS_HPP_
#define STOCHASTICCOLLOCATIONS_HPP_

#include <boost/random.hpp>
#include <iostream>
#include <numeric>
#include <cmath>
#include <cstring>
#include <fstream>
#include <mpi.h>


 class StochasticCollocations
 {
 private:
    /* the Raynolds number will be of the form mean + stddev*sigma, where
    sigma is ~N(0,1) or ~U(0,1) */
    double u_gauss, s_gauss; // u_gauss = 0, s_gauss = 1;
    int u_uniform , s_uniform; // u_uniform = 0; s_uniform = 1;

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
    StochasticCollocations(double u_gauss, double s_gauss);
    StochasticCollocations(int u_uniform, int s_uniform);

    /** Uncertainty (i.e. Random variables) related methods **/

    /* generate nsamples samples of normal distributed random variables */
    std::vector<double> generate_nd_samples(double mean_nd, double sttdev_nd, int nsamples);
    /* generate nsamples samples of uniformly distributed random variables */
    std::vector<double> generate_ud_samples(double mean_ud, double sttdev_ud, int nsamples);

    /* compute the first two statistical moments */
    double compute_mean(const std::vector<double> &v) const;
    double compute_variance(const std::vector<double> &v, double mean) const;

    /* get a normal and uniform distributed RV */
    double get_normal() const;
    double get_uniform() const;

    double hermite_poly(int degree, double &var);
    void gauss_hermite_quad(int quad_degree, std::vector<double> &nodes, std::vector<double> &weights);
    std::vector<double> &get_coefficiants(int quad_degree, int no_coeff, std::vector<double> &nodes, 
        std::vector<double> &weights, int rv_choice, int &imax, int &jmax);

    /***********************************************************/

    /** Parallelization related methods **/

    /* data decomposition among processes*/
    void data_decomposition(int* nsampels, int* nprocs, int* samples_per_proc);

    /* call the NS solver for each generated sample */
    void get_NS_solution(int* samples_per_proc, const std::vector<double> &Re);

    /* get the QoI (Quantities of interest - the desired output parameters, from a UQ point of view) */
    void get_QoI(int* samples_per_proc, std::vector<double> &qoi);

    /* helper function for ending the mpi parallelization */
    void SCSimulation_Stop(const std::string &txt);

    /***********************************************************/

    /* destructor */
    ~StochasticCollocations();
};

#endif /* MONTECARLO_HPP_ */