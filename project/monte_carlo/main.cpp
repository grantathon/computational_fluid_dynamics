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
    const int x_uniform = 0.0;
    const int y_uniform = 1.0;

    /* statistcs for GRV*/
    double expectation_mc = 0.0;
    double var_mc = 0.0;

    const char* file_name = "UQ_simulation.dat";
    const std::string eos = "End of simulation";

    /* write data for ploting purposes */
    std::ofstream output_file;
    output_file.open(file_name);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    
    output_file << "x_reattach" << " " << "t_reattach" << std::endl;
    
    
    if(read_parameters(argc, argv, &flag_uq, &flag_distr, &flag_rv, &nsamples, &mean, &stddev, &imax, &jmax) == 0)
    {
        return 0;
    } 

    /* samples per processor*/
    samples_per_proc = (nsamples / num_proc);

    /*printf("UQ: %i \t dist: %i \t RV: %i \t N: %i \t mean: %f \t var: %f \t imax: %i \t jmax: %i \n", flag_uq, flag_distr, flag_rv, nsamples, mean, stddev, imax, jmax);*/

    if(myrank == 0) 
    {
        start_time = MPI_Wtime();

        // Construct PBM file based on user input
        PBMFile *pbm = new PBMFile(imax, jmax, 0, 1, "flow_over_a_step.pbm");
        pbm->OutputStep(0.5, 0.1);
        delete pbm;
    }

    MonteCarlo *m;

    std::vector<double> samples;
    std::vector<double> qoi;

    /* MC + Gaussian distribution + Re */
    if(flag_uq == 1 && flag_distr == 1 && flag_rv == 1)
    {
        m = new MonteCarlo(u_gauss, s_gauss);
        m->data_decomposition(samples_per_proc, &nsamples, &num_proc, &myrank, &il, &ir);
        printf("rank: %i \t il: %i \t ir: %i \n", myrank, il, ir);

        samples = m->generate_nd_samples(mean, stddev, &nsamples);

        m->get_NS_solution(&myrank, &samples_per_proc, samples, &il, &ir, output_file, qoi, flag_rv, imax, jmax);

        /*if(myrank == 0)
        {
            for(int i = 0 ; i < nsamples ; i++)
            {
                std::cout << "qoi[" << i << "]=" << qoi[i] << std::endl;
            }
        }*/
        
        m->get_QoI(&myrank, &samples_per_proc, &il, &ir, qoi, output_file);
    }

    
    if(flag_uq == 1 && flag_distr == 1 && flag_rv != 1)
    {
        m = new MonteCarlo(u_gauss, s_gauss);
        m->data_decomposition(samples_per_proc, &nsamples, &num_proc, &myrank, &il, &ir);
        printf("rank: %i \t il: %i \t ir: %i \n", myrank, il, ir);

        samples = m->generate_nd_samples(mean, stddev, &nsamples);

        m->get_NS_solution(&myrank, &samples_per_proc, samples, &il, &ir, output_file, qoi, flag_rv, imax, jmax);

        m->get_QoI(&myrank, &samples_per_proc, &il, &ir, qoi, output_file);     
    }


    else if(flag_uq == 1 && flag_distr != 1 && flag_rv == 1)
    {
        m = new MonteCarlo(x_uniform, y_uniform);
        m->data_decomposition(samples_per_proc, &nsamples, &num_proc, &myrank, &il, &ir);
        printf("rank: %i \t il: %i \t ir: %i \n", myrank, il, ir);

        samples = m->generate_ud_samples(mean, stddev, &nsamples);

        m->get_NS_solution(&myrank, &samples_per_proc, samples, &il, &ir, output_file, qoi, flag_rv, imax, jmax);

        m->get_QoI(&myrank, &samples_per_proc, &il, &ir, qoi, output_file); 
    }

    
    if(flag_uq == 1 && flag_distr != 1 && flag_rv != 1)
    {
        m = new MonteCarlo(x_uniform, y_uniform);
        m->data_decomposition(samples_per_proc, &nsamples, &num_proc, &myrank, &il, &ir);
        printf("rank: %i \t il: %i \t ir: %i \n", myrank, il, ir);

        samples = m->generate_ud_samples(mean, stddev, &nsamples);

        m->get_NS_solution(&myrank, &samples_per_proc, samples, &il, &ir, output_file, qoi, flag_rv, imax, jmax);

        m->get_QoI(&myrank, &samples_per_proc, &il, &ir, qoi, output_file);    
    }


    if(myrank == 0)
    { 
        expectation_mc = m->compute_mean(qoi);
        var_mc = m->compute_variance(qoi, expectation_mc);
        std::cout << "The mean of the re-attachement point is: " << expectation_mc << std::endl;
        std::cout << "The variance of the re-attachement point is:  " << var_mc << std::endl;
    }

    if (myrank == 0)
    {
        end_time = MPI_Wtime();
        printf("Elapsed time for MC simulation (%d samples) using %d processors is: %f seconds\n", nsamples, num_proc, end_time - start_time);
    }

    Simulation_Stop(eos);
    output_file.close();

    delete m;
    return 0;
}
