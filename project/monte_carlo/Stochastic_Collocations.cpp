/*
 * Stochastic_Collocations.cpp
 *
 *  Created on: Jun 27, 2014
 *      Author: ionut
 */

#include "Stochastic_Collocations.hpp"


 StochasticCollocations::StochasticCollocations(double u_gauss, double s_gauss)
 {
 	rng = boost::mt19937();

 	normal_distr = boost::normal_distribution<>(u_gauss, s_gauss);

 	var_normal = new boost::variate_generator<boost::mt19937&,
 	boost::normal_distribution<> >(rng, normal_distr);
 }

 StochasticCollocations::StochasticCollocations(int u_uniform, int s_uniform)
 {
 	rng = boost::mt19937();

 	uniform_distr = boost::uniform_int<>(u_uniform, s_uniform);

 	var_uniform = new boost::variate_generator<boost::mt19937&,
 	boost::uniform_int<> >(rng, uniform_distr);
 }

 double StochasticCollocations::get_normal() const {
 	return (*var_normal)();
 }

 double StochasticCollocations::get_uniform() const {
 	return (*var_uniform)();
 }

 double StochasticCollocations::hermite_poly(int degree, double &var)
 {
 	double poly_eval;
 	switch(degree)
 	{
 		case 0:
 		poly_eval = 1.0;
 		break;
 		case 1:
 		poly_eval = var;
 		break;
 		case 2:
 		poly_eval = pow(var, 2) - 1;
 		break;
 		case 3:
 		poly_eval = pow(var, 3) - 3*var;
 		break;
 		case 4:
 		poly_eval = pow(var, 4) - 6*pow(var, 2) + 3;
 		break;
 		case 5:
 		poly_eval = pow(var, 5) - 10*pow(var, 3) + 15*var;
 		break;
 		case 6:
 		poly_eval = pow(var, 6) - 15*pow(var, 4) + 45*pow(var, 2) - 15;
 		break;
 		case 7:
 		poly_eval = pow(var, 7) - 21*pow(var, 5) + 105*pow(var, 3) - 105*var;
 		break;
 		case 8:
 		poly_eval = pow(var, 8) - 28*pow(var, 6) + 210*pow(var, 4) - 420*pow(var, 2) + 105;
 		break;
 		case 9:
 		poly_eval = pow(var, 9) - 36*pow(var, 7) + 378*pow(var, 5) - 1260*pow(var, 3) + 945*var;
 		break;
 		case 10:
 		poly_eval = pow(var, 10) - 45*pow(var, 8) + 630*pow(var, 6) - 3150*pow(var, 4) + 4725*pow(var, 2) - 945;
 		break;
 		default:
 		std::cout << "Please input a degree less or equal to 10!" << std::endl;
 	}

 	return poly_eval;
 }

 void StochasticCollocations::gauss_hermite_quad(int quad_degree, std::vector<double> &nodes, std::vector<double> &weights)
 {
 	switch(quad_degree)
 	{
 		case 1:
 		nodes.push_back(0.0);
 		weights.push_back(1.7725);
 		break;

 		case 2:
 		nodes.push_back(-0.7071);
 		nodes.push_back(0.7071);
 		weights.push_back(0.8862);
 		weights.push_back(0.8862);
 		break;

 		case 3:
 		nodes.push_back(-1.2247);
 		nodes.push_back(0.0000);
 		nodes.push_back(1.2247);
 		weights.push_back(0.2954);
 		weights.push_back(1.1816);
 		weights.push_back(0.2954);
 		break;

 		case 4:
 		nodes.push_back(1.6507);
 		nodes.push_back(-0.5246);
 		nodes.push_back(0.5246);
 		nodes.push_back(1.6507);
 		weights.push_back(0.0813);
 		weights.push_back(0.8049);
 		weights.push_back(0.8049);
 		weights.push_back(0.0813);
 		break;

 		case 5:
 		nodes.push_back(-2.0202);
 		nodes.push_back(-0.9586);
 		nodes.push_back(0.0);
 		nodes.push_back(0.9586);
 		nodes.push_back(2.0202);
 		weights.push_back(0.0200);
 		weights.push_back(0.3936);
 		weights.push_back(0.9453);
 		weights.push_back(0.3936);
 		weights.push_back(0.0200);
 		break;

 		case 6:
 		nodes.push_back(-2.3506);
 		nodes.push_back(-1.3358);
 		nodes.push_back(-0.4361);
 		nodes.push_back(0.4361);
 		nodes.push_back(1.3358);
 		nodes.push_back(2.3506);
 		weights.push_back(0.0045);
 		weights.push_back(0.1571);
 		weights.push_back(0.7246);
 		weights.push_back(0.7246);
 		weights.push_back(0.1571);
 		weights.push_back(0.0045);
 		break;

 		case 7:
 		nodes.push_back(-2.6520);
 		nodes.push_back(-1.6736);
 		nodes.push_back(-0.8163);
 		nodes.push_back(0.0000);
 		nodes.push_back(0.8163);
 		nodes.push_back(1.6736);
 		nodes.push_back(2.6520);
 		weights.push_back(0.0010);
 		weights.push_back(0.0545);
 		weights.push_back(0.4256);
 		weights.push_back(0.8103);
 		weights.push_back(0.4256);
 		weights.push_back(0.0545);
 		weights.push_back(0.0010);
 		break;

 		case 8:
 		nodes.push_back(-2.9306);
 		nodes.push_back(-1.9817);
 		nodes.push_back(-1.1572);
 		nodes.push_back(-0.3812);
 		nodes.push_back(0.3812);
 		nodes.push_back(1.1572);
 		nodes.push_back(1.9817);
 		nodes.push_back(2.9306);
 		weights.push_back(0.0002);
 		weights.push_back(0.0171);
 		weights.push_back(0.2078);
 		weights.push_back(0.6611);
 		weights.push_back(0.6611);
 		weights.push_back(0.2078);
 		weights.push_back(0.0171);
 		weights.push_back(0.0002);
 		break;

 		case 9:
 		nodes.push_back(-3.1910);
 		nodes.push_back(-2.2666);
 		nodes.push_back(-1.4686);
 		nodes.push_back(-0.7236);
 		nodes.push_back(0.0);
 		nodes.push_back(0.7236);
 		nodes.push_back(1.4686);
 		nodes.push_back(2.2666);
 		nodes.push_back(3.1910);
 		weights.push_back(0.0000);
 		weights.push_back(0.0049);
 		weights.push_back(0.0885);
 		weights.push_back(0.4327);
 		weights.push_back(0.7202);
 		weights.push_back(0.4327);
 		weights.push_back(0.0885);
 		weights.push_back(0.0049);
 		weights.push_back(0.0000);
 		break;

 		case 10:
 		nodes.push_back(-3.4362);
 		nodes.push_back(-2.5327);
 		nodes.push_back(-1.7567);
 		nodes.push_back(-1.0366);
 		nodes.push_back(-0.3429);
 		nodes.push_back(0.3429);
 		nodes.push_back(1.0366);
 		nodes.push_back(1.7567);
 		nodes.push_back(2.5327);
 		nodes.push_back(3.4362);
 		weights.push_back(0.0000);
 		weights.push_back(0.0013);
 		weights.push_back(0.0339);
 		weights.push_back(0.2401);
 		weights.push_back(0.6109);
 		weights.push_back(0.6109);
 		weights.push_back(0.2401);
 		weights.push_back(0.0339);
 		weights.push_back(0.0013);
 		weights.push_back(0.0000);
 		break;

 		default:
 		std::cout << "Please input a degree less or equal to 10" << std::endl;
 	}
 }
 std::vector<double> StochasticCollocations::get_coefficiants(int quad_degree, int no_coeff, double mean, double stddev, std::vector<double> &nodes, std::vector<double> &weights)
 {
 	std::vector<double> coeff;
 	double temp = 0.0;
 	double var = 0.0;

 	char buffer[15];

 	for(int j = 0 ; j < no_coeff ; j++)
 	{
 		double nodes_quad, x_global, t_reattach;
 		char datafile_name[30] = "ns_sim_";

 		snprintf(buffer, sizeof(buffer), "%d%s", j + 1, ".mc");	
 		strcat(datafile_name, buffer);

 		std::ifstream MC_data(datafile_name);
 		MC_data >> nodes_quad >> x_global >> t_reattach;

 		for(int i = 0 ; i < quad_degree ; i++)
 		{
 			var = sqrt(2)*stddev*nodes[i] + mean;
 			temp = temp + t_reattach * hermite_poly(j, var) * weights[i];
 		}
 		temp = temp/M_PI;
 		coeff.push_back(temp);
 	}
 	return coeff;
 }

 void StochasticCollocations::data_decomposition(int* ncoeff, int* nprocs, int* coeff_per_proc)
 {
 	int rank;

 	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

 	if (rank == 0) 
 	{
 		*coeff_per_proc = *ncoeff - (*nprocs - 1)*((*ncoeff)/(*nprocs));
 	} 
 	else 
 	{
 		*coeff_per_proc = (*ncoeff)/(*nprocs);
 	}
 }

void StochasticCollocations::get_NS_solution(int* coeff_per_proc, const std::vector<double> &nodes, int rv_flag, int imax, int jmax)
 {
 	int rank;
 	char buffer[20];

 	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

 	for(int i = 0 ; i < *coeff_per_proc ; i++)
 	{
 		char call_NS_solver[30] = "mpirun -np 4 ./sim ";
 		snprintf(buffer, sizeof(buffer), "%d %g %d %d %d", rv_flag, nodes[i], rank*(*coeff_per_proc) + i + 1, imax, jmax);	
 		strcat(call_NS_solver, buffer);

 		system(call_NS_solver);
 	}
 }

  StochasticCollocations::~StochasticCollocations() {
 	delete var_normal;
 	delete var_uniform;
 }	
 