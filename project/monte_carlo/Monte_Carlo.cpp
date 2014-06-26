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
 		temp = fabs(get_normal());
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

 double MonteCarlo::compute_mean(const std::vector<double> &v) const {
 	double sum = 0.0, mean = 0.0;

 	sum = std::accumulate(v.begin(), v.end(), 0.0);
 	mean = sum / v.size();

 	return mean;
 }

 double MonteCarlo::compute_variance(const std::vector<double> &v, double mean) const {
 	double square_sum = 0.0, stddev = 0.0;

 	std::vector<double> diff(v.size());
 	std::transform(v.begin(), v.end(), diff.begin(),
 		std::bind2nd(std::minus<double>(), mean));
 	square_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(),
 		0.0);
 	stddev = square_sum / (v.size() - 1);

 	return stddev;
 }


 void MonteCarlo::data_decomposition(int* nsampels, int* nprocs, int *samples_per_proc)
 {
 	int rank;

 	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

 	if (rank == 0) 
 	{
 		*samples_per_proc = *nsampels - (*nprocs - 1)*((*nsampels)/(*nprocs));
 	} 
 	else 
 	{
 		*samples_per_proc = (*nsampels)/(*nprocs);
 	}
 }

 void MonteCarlo::get_NS_solution(int* samples_per_proc, const std::vector<double> &Re)
 {
 	int rank;
 	char buffer[10];

 	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

 	for(int i = 0 ; i < *samples_per_proc ; i++)
	{
		char call_NS_solver[30] = "mpirun -np 4 ./sim ";
		snprintf(buffer, sizeof(buffer), "%g %d", Re[i], rank*(*samples_per_proc) + i + 1);	
		strcat(call_NS_solver, buffer);

		system(call_NS_solver);

		//MPI_Barrier( MPI_COMM_WORLD);
	}
 }

 void MonteCarlo::get_QoI(int* samples_per_proc, std::vector<double> &qoi)
 {
 	int rank;
 	char buffer[10];

 	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

 	for(int i = 0 ; i < *samples_per_proc ; i++)
	{
		double Re, x_global, t_reattach;
		char datafile_name[30] = "ns_sim_";

		snprintf(buffer, sizeof(buffer), "%d%s", rank*(*samples_per_proc) + i + 1, ".mc");	
		strcat(datafile_name, buffer);

		std::ifstream MC_data(datafile_name);
		MC_data >> Re >> x_global >> t_reattach;

		qoi.push_back(t_reattach);

		//MPI_Barrier( MPI_COMM_WORLD);
	}
 }

 void MonteCarlo::MCSimulation_Stop(const std::string &txt)
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

 MonteCarlo::~MonteCarlo() {
 	delete var_normal;
 	delete var_uniform;
 }