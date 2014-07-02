#include "Monte_Carlo.hpp"

 MonteCarlo::MonteCarlo(double u_gauss, double s_gauss)
 {
 	rng = boost::mt19937();

 	normal_distr = boost::normal_distribution<>(u_gauss, s_gauss);

 	var_normal = new boost::variate_generator<boost::mt19937&,
 	boost::normal_distribution<> >(rng, normal_distr);
 }

 MonteCarlo::MonteCarlo(int u_uniform, int s_uniform)
 {
 	rng = boost::mt19937();

 	uniform_distr = boost::uniform_int<>(u_uniform, s_uniform);

 	var_uniform = new boost::variate_generator<boost::mt19937&,
 	boost::uniform_int<> >(rng, uniform_distr);
 }

 double MonteCarlo::get_normal() const {
 	return (*var_normal)();
 }

 double MonteCarlo::get_uniform() const {
 	return (*var_uniform)();
 }

 std::vector<double> MonteCarlo::generate_nd_samples(double mean_nd,
 	double stddev_nd, int nsamples) {
 	std::vector<double> ndsamples;
 	double temp;

 	for (int i = 0; i < nsamples; i++) {
 		temp = mean_nd + stddev_nd*get_normal();
 		ndsamples.push_back(temp);
 	}

 	return ndsamples;
 }

 std::vector<double> MonteCarlo::generate_ud_samples(double mean_ud, double stddev_ud,
 	int nsamples) {
 	std::vector<double> udsamples;
 	double temp;

 	for (int i = 0; i < nsamples; i++) {
 		temp = mean_ud + stddev_ud*get_uniform();
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

 void MonteCarlo::get_NS_solution(int* samples_per_proc, const std::vector<double> &Re, int rv_flag, int imax, int jmax)
 {
 	int rank;
 	char buffer[20];

 	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

 	for(int i = 0 ; i < *samples_per_proc ; i++)
 	{
 		char call_NS_solver[30] = "./sim ";
 		snprintf(buffer, sizeof(buffer), "%d %g %d %d %d", rv_flag, Re[i], rank*(*samples_per_proc) + i + 1, imax, jmax);	
 		strcat(call_NS_solver, buffer);

 		system(call_NS_solver);
 	}
 }

 void MonteCarlo::get_QoI(int* samples_per_proc, std::vector<double> &qoi)
 {
 	int rank;
 	char buffer[15];

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
 	}
 }

 MonteCarlo::~MonteCarlo() {
 	delete var_normal;
 	delete var_uniform;
 }
