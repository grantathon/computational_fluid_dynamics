#include "Monte_Carlo.hpp"
#include "Stochastic_Collocations.hpp"
#include "helper.hpp"
#include "PBM_File.hpp"

int main(int argc, char *argv[]) 
{
    /* MPI variables */
    int num_proc = 0;
    int myrank = 0;

    double start_time = 0.0, end_time = 0.0;

    int samples_per_proc = 0;

    int nsamples = 0;
    double mean_nd = 0.0, stddev_nd = 0.0;
    int uq_method = 0 , distribution = 0;
    int imax = 0, jmax = 0;
    int rv_flag = 1;

    /* for Gaussian distribution - consider only ~N(0,1)*/
    const double u_gauss = 0.0;
    const double s_gauss = 1.0;

    /* statistcs */
    double expectation_mc = 0.0;
    double var_mc = 0.0;

    std::string eos = "End of simulation";

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    if(argc != 7)
    {
        uq_method = atoi(argv[1]);
        distribution = atoi(argv[2]);
        nsamples = atoi(argv[3]);
        mean_nd = atof(argv[4]);
        stddev_nd = atof(argv[5]);
        imax = atoi(argv[6]);
        jmax = atoi(argv[7]);
    }
    else
    {
        Simulation_Stop("Run the program as mpirun -np x ./sim_UQ arg1 arg2 arg3 arg4 arg5, where arg1 = mc or sc, arg 2 = g or u arg3 is the number of samples(mc) or coeff(sc) and arg4 and arg5 are the mean and stddev of Re ");
        return 0;
    }
    

    std::cout << uq_method << distribution << std::endl;

    if(myrank == 0) 
    {
        start_time = MPI_Wtime();

        // Construct PBM file based on user input
        PBMFile *pbm = new PBMFile(imax, jmax, 0, 1, "flow_over_a_step.pbm");
        pbm->OutputStep(0.5, 0.1);
        delete pbm;
    }

    MonteCarlo *m = new MonteCarlo(u_gauss, s_gauss);

    std::vector<double> nd_samples;
    std::vector<double> nd_qoi;

    nd_samples = m->generate_nd_samples(mean_nd, stddev_nd, nsamples);
    m->data_decomposition(&nsamples, &num_proc, &samples_per_proc);
    std:: cout << samples_per_proc << std::endl;
    m->get_NS_solution(&samples_per_proc, nd_samples, rv_flag, imax, jmax);
    m->get_QoI(&samples_per_proc, nd_qoi);

    if(myrank == 0)
    {   
        expectation_mc = m->compute_mean(nd_qoi);
        var_mc = m->compute_variance(nd_qoi, expectation_mc);

        std::cout << "The mean of the re-attachement point is: " << expectation_mc << std::endl;
        std::cout << "The variance of the re-attachement point is:  " << var_mc << std::endl;
    }

    if (myrank == 0)
    {
        end_time = MPI_Wtime();
        printf("Elapsed time for MC simulation (%d samples) using %d processors is: %f seconds\n", nsamples, num_proc, end_time - start_time);
    }

    Simulation_Stop(eos);

    delete m;

    return 0;
}
