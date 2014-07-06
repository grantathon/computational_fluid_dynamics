#include "Monte_Carlo.hpp"
#include "Stochastic_Collocations.hpp"
#include "helper.hpp"
#include "PBM_File.hpp"

int main(int argc, char *argv[]) 
{
    /* MPI variables */
    int num_proc = 0;
    int myrank = 0;
    int il = 0;
    int ir = 0;

    /* Runtime measurement */
    double start_time = 0.0;
    double end_time = 0.0;

    /* Command line arguments */
    int flag_uq = 0;
    int flag_distr = 0;
    int flag_rv = 0;
    int nsamples = 0;
    double mean = 0.0;
    double stddev = 0.0; 
    int imax = 0;
    int jmax = 0;

    /* Data decomposition */
    int samples_per_proc = 0;

    /* for Gaussian distribution*/
    const double u_gauss = 0.0;
    const double s_gauss = 1.0;

    /* for Uniform distribution */
    const double x_uniform = 0.0;
    const double y_uniform = 1.0;

    /* statistics*/
    double expectation_mc = 0.0;
    double var_mc = 0.0;

    /* for finalizing MPI */
    const std::string eos = "End of simulation";
    /* used for printing purposes */
    /*const char* data_file = "UQ_simulation.cvs";*/

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    /* read the command line arguments */
    if(read_parameters(argc, argv, &flag_uq, &flag_distr, &flag_rv, &nsamples, &mean, &stddev, &imax, &jmax) == 0)
    {
        return 0;
    } 

    /* samples per processor*/
    samples_per_proc = (nsamples / num_proc);

    if(myrank == 0) 
    {
        start_time = MPI_Wtime();

        PBMFile *pbm = new PBMFile(imax, jmax, 0, 1, "flow_over_a_step.pbm");
        pbm->OutputStep(0.5, 0.1);
        delete pbm;
    }

    /* Declare a MC object */
    MonteCarlo *m;
    std::vector<double> samples;

    /* MC + Gaussian distribution + Re */
    if(flag_uq == 1 && flag_distr == 1 && flag_rv == 1)
    {
        m = new MonteCarlo(u_gauss, s_gauss, flag_distr);
        m->data_decomposition(samples_per_proc, &nsamples, &num_proc, &myrank, &il, &ir);
        samples = m->generate_nd_samples(mean, stddev, &nsamples);
        m->monte_carlo_simulation(&myrank, &nsamples, &samples_per_proc, &il, &ir, samples, flag_rv, imax, jmax, &expectation_mc, &var_mc);

        if(myrank == 0)
        {
            std::cout << "The mean of the re-attachement(when Re is normal distributed) point is: " << expectation_mc << std::endl;
            std::cout << "The variance of the re-attachement(when Re is normal distributed) point is:  " << var_mc << std::endl;
        }
    }

    /* MC + Gaussian distribution + viscosity */
    if(flag_uq == 1 && flag_distr == 1 && flag_rv != 1)
    {
        m = new MonteCarlo(u_gauss, s_gauss, flag_distr);
        m->data_decomposition(samples_per_proc, &nsamples, &num_proc, &myrank, &il, &ir);
        samples = m->generate_nd_samples(mean, stddev, &nsamples);
        m->monte_carlo_simulation(&myrank, &nsamples, &samples_per_proc, &il, &ir, samples, flag_rv, imax, jmax, &expectation_mc, &var_mc);

        if(myrank == 0)
        {
            std::cout << "The mean of the re-attachement(when viscosity is normal distributed) point is: " << expectation_mc << std::endl;
            std::cout << "The variance of the re-attachement point(when viscosity is normal distributed) is:  " << var_mc << std::endl;
        }     
    }

    /* MC + Uniform distribution + Re */
    else if(flag_uq == 1 && flag_distr != 1 && flag_rv == 1)
    {
        m = new MonteCarlo(x_uniform, y_uniform, flag_distr);
        m->data_decomposition(samples_per_proc, &nsamples, &num_proc, &myrank, &il, &ir);
        samples = m->generate_ud_samples(mean, stddev, &nsamples);
        m->monte_carlo_simulation(&myrank, &nsamples, &samples_per_proc, &il, &ir, samples, flag_rv, imax, jmax, &expectation_mc, &var_mc);

        if(myrank == 0)
        {
            std::cout << "The mean of the re-attachement(when Re is uniform distributed) point is: " << expectation_mc << std::endl;
            std::cout << "The variance of the re-attachement point(when Re is uniform distributed) is:  " << var_mc << std::endl;
        }      
    }

    /* MC + Uniform distribution + viscosity */
    else if(flag_uq == 1 && flag_distr != 1 && flag_rv != 1)
    {
         m = new MonteCarlo(x_uniform, y_uniform, flag_distr);
         m->data_decomposition(samples_per_proc, &nsamples, &num_proc, &myrank, &il, &ir);
         samples = m->generate_ud_samples(mean, stddev, &nsamples);
         m->monte_carlo_simulation(&myrank, &nsamples, &samples_per_proc, &il, &ir, samples, flag_rv, imax, jmax, &expectation_mc, &var_mc);

         if(myrank == 0)
         {
            std::cout << "The mean of the re-attachement(when viscosity is uniform distributed) point is: " << expectation_mc << std::endl;
            std::cout << "The variance of the re-attachement point(when viscosity is uniform distributed) is:  " << var_mc << std::endl;
          }        
    }


    if (myrank == 0)
    {
        /*write_data_file(data_file, &nsamples);*/
        end_time = MPI_Wtime();
        printf("Elapsed time for MC simulation (%d samples) using %d processors is: %f seconds\n", nsamples, num_proc, end_time - start_time);
    }

    Simulation_Stop(eos);
    delete m;

    return 0;
}
