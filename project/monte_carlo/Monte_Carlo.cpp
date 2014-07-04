#include "Monte_Carlo.hpp"
#include <stdio.h>

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
	double stddev_nd, int* nsamples)
{
	std::vector<double> ndsamples;
	double temp;

	for (int i = 0  ; i < *nsamples ; i++) 
	{
		temp = fabs(mean_nd + stddev_nd*get_normal());
		ndsamples.push_back(temp);
		/*printf("rank: %i \t sample: %f \n", rank, temp);*/
	}

	return ndsamples;
}


std::vector<double> MonteCarlo::generate_ud_samples(double mean_ud, double stddev_ud,
	int *nsamples) {
	std::vector<double> udsamples;
	double temp;

	for (int i = 0; i < *nsamples; i++) {
		temp = fabs(mean_ud + stddev_ud*get_uniform());
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


void MonteCarlo::data_decomposition(int samples_per_proc, int* nsampels, int* nprocs, int *myrank, int *il, int *ir)
{
	int i_rem;

	/* Compute il, ir for each sub-domain*/
	i_rem = (*nsampels) % (*nprocs);

	*il = ((*myrank) * samples_per_proc);
	if ( *myrank == ((*nprocs) - 1))   /* make last rank the monḱey and put a banana in its mouth*/
	{
		*ir = ((*myrank + 1) * samples_per_proc) + i_rem - 1;
	}
	else
	{
		*ir = ((*myrank + 1) * samples_per_proc) - 1;
	}
	/*printf("rank: %i \t il: %i \t ir: %i \n", *myrank, *il, *ir);*/
}


void MonteCarlo::get_NS_solution(int *myrank, int* samples_per_proc, std::vector<double> &Re, int *il, int *ir,
									std::ofstream &outputFile, std::vector<double> &qoi, int rv_flag, int imax, int jmax)
{
	int sample_size, global_id;
	char buffer_solver[20];
	char buffer_datafile[20];

	sample_size = (*ir) - (*il) + 1;

	for(int i = 0; i < sample_size ; i++)
	{
		global_id = i + (*myrank) * (*samples_per_proc);
		
		char call_NS_solver[30] = "./sim ";
		/*char call_NS_solver[30] = "./t_rex_balls ";*/
		
		printf("rank: %i \t NS_inputs: \t %d %g %d %d %d \n", *myrank, rv_flag, Re[global_id], global_id + 1, imax, jmax);

		snprintf(buffer_solver, sizeof(buffer_solver), "%u %g %u %u %u", rv_flag, Re[global_id], global_id + 1, imax, jmax);
		strcat(call_NS_solver, buffer_solver);

		system(call_NS_solver);

		double Re, x_global, t_reattach;
		char datafile_name[30] = "ns_sim_";

		snprintf(buffer_datafile, sizeof(buffer_datafile), "%d%s", global_id + 1, ".mc");	
		strcat(datafile_name, buffer_datafile);

		std::ifstream MC_data(datafile_name);
		MC_data >> Re >> x_global >> t_reattach;

		/*printf("rank: %i \t %s \t \n", rank, buffer_datafile);*/

		write_to_file(outputFile, x_global, t_reattach);
		qoi.push_back(t_reattach);
		
	}
}

void MonteCarlo::get_QoI(int *myrank, int *samples_per_proc, int *il, int *ir, std::vector<double> &qoi, std::ofstream &outputFile)
{
	int sample_size, global_id;

	sample_size = (*ir) - (*il) + 1;

	char buffer[15];
	
/*	for(int i = rank*(*samples_per_proc) ; i < (rank + 1)*(*samples_per_proc) ; i++) */
	for(int i = 0; i < sample_size ; i++)
	{
		global_id = i + (*myrank) * (*samples_per_proc);

		double Re, x_global, t_reattach;
		char datafile_name[30] = "ns_sim_";

		snprintf(buffer, sizeof(buffer), "%d%s", global_id + 1, ".mc");	
		strcat(datafile_name, buffer);

		std::ifstream MC_data(datafile_name);
		MC_data >> Re >> x_global >> t_reattach;
		qoi.push_back(t_reattach);

	}
	MPI_Barrier(MPI_COMM_WORLD);
}

MonteCarlo::~MonteCarlo() {
	delete var_normal;
	delete var_uniform;
}
