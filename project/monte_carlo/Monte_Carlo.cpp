#include "Monte_Carlo.hpp"
#include <stdio.h>

MonteCarlo::MonteCarlo(double par1, double par2, int distr_flag)
{
	rng = boost::mt19937();

	if(distr_flag == 1)
	{
		normal_distr = boost::normal_distribution<>(par1, par2);

		var_normal = new boost::variate_generator<boost::mt19937&,
		boost::normal_distribution<> >(rng, normal_distr);
	}
	else
	{
		uniform_distr = boost::uniform_real<>(par1, par2);

		var_uniform = new boost::variate_generator<boost::mt19937&,
		boost::uniform_real<> >(rng, uniform_distr);
	}
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


void MonteCarlo::data_decomposition(int samples_per_proc, int* nsampels, int* nprocs, int *myrank, int *il, int *ir)
{
	int i_rem;

	/* Compute il, ir for each sub-domain*/
	i_rem = (*nsampels) % (*nprocs);

	/*	sprinkle the remainder over all procs */
	if ( *myrank < i_rem)   
	{
		*il = ((*myrank) * (samples_per_proc + 1));
		*ir = (*il) + samples_per_proc;
	}
	else
	{
		*il = ((*myrank) * samples_per_proc) + i_rem;
		*ir = (*il) + samples_per_proc - 1;
	}

	printf("rank: %i \t il: %i \t ir: %i \n", *myrank, *il, *ir);
}


void MonteCarlo::monte_carlo_simulation(int *myrank, int* nsamples, int* samples_per_proc, int *il, int *ir
	, std::vector<double> &rv, int rv_flag, int imax, int jmax, double *mean, double* variance)
{
	int sample_size, global_id;
	char buffer_solver[30];
	char buffer_datafile[20];

	sample_size = (*ir) - (*il) + 1;
	double partial_sum = 0.0;
	double sum = 0.0;

	double partial_sum2 = 0.0;
	double sum2 = 0.0;

	double Re, x_global, t_reattach;

	for(int i = 0; i < sample_size ; i++)
	{
		/*********************************************/
		/* call the NS solver */
		global_id = i + (*myrank) * (*samples_per_proc);
		
		char call_NS_solver[30] = "./sim ";
		snprintf(buffer_solver, sizeof(buffer_solver), "%d %g %d %d %d", rv_flag, rv[global_id], global_id + 1, imax, jmax);	
		strcat(call_NS_solver, buffer_solver);

		system(call_NS_solver);
		/*********************************************/

		/*********************************************/
		/* obtain the QoI from the data files */
		char datafile_name[30] = "ns_sim_";
		snprintf(buffer_datafile, sizeof(buffer_datafile), "%d%s", global_id + 1, ".mc");	
		strcat(datafile_name, buffer_datafile);

		std::ifstream MC_data(datafile_name);
		MC_data >> Re >> x_global >> t_reattach;
		/*********************************************/

		/*********************************************/
		/* compute the partial sum for each process - used to compute the mean */
		partial_sum += t_reattach;
		/*********************************************/

		/*********************************************/
		/* compute the partial sum squared for each process - used to compute the variance */
		partial_sum2 += t_reattach*t_reattach;
		/*********************************************/
	}

	/* compute global mean */
	MPI_Allreduce(&partial_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	*mean = sum/(*nsamples);

	/* compute global variance */
	MPI_Allreduce(&partial_sum2, &sum2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	*variance = (sum2 - (sum*sum)/(*nsamples))/(*nsamples - 1.0);
}

MonteCarlo::~MonteCarlo() {
	delete var_normal;
	delete var_uniform;
}
