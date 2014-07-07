#include "helper.hpp"

void write_data_file(const char* file_name, int* nsampels)
{
   char buffer_datafile[15];
   double Re, x_global, t_reattach;

   std::ofstream data_file;

   data_file.open(file_name);

   data_file << "x_reattach" << " " << "t_reattach" << std::endl;

   for(int i = 0 ; i < *nsampels ; i++)
   {
      char datafile_name[30] = "ns_sim_";
      snprintf(buffer_datafile, sizeof(buffer_datafile), "%d%s", i + 1, ".mc");   
      strcat(datafile_name, buffer_datafile);

      std::ifstream MC_data(datafile_name);
      MC_data >> Re >> x_global >> t_reattach;

      data_file << x_global << ", " << t_reattach << std::endl;
   }

   data_file.close();
}

int read_parameters(int argc, char** argv, int *flag_UQ, int* flag_distr, int* flag_RV, int* flag_prog, int *npoints, double* mean, double* stddev, int* imax, int *jmax)
{
   if(argc == 10)
   {
      /*flag_UQ flag_distr flag_RV npoints mean stddev imax jmax*/
      /* 1 1 1 5 0 1 10 2*/
   /* flag_UQ = 1 for MC, everything else (>=0) for SC*/
      *flag_UQ = atoi(argv[1]);
   /* flag_distr = 1 for Gaussian distribution, everything else (>=0) for normal distribution */
      *flag_distr = atoi(argv[2]);
   /* flag_RV = 1 for Reynolds number, everything else (>=0) for viscosity */
      *flag_RV = atoi(argv[3]);
   /* flag_prog = 1 for serial, everything else for parallel */
      *flag_prog = atoi(argv[4]);
   /* number of samples (if MC) or no of coefficiants (SC) */
      *npoints = atoi(argv[5]);
   /* mean for the RV ; both are of the form mean + zeta*stddev, where zeta is ~N(0,1) or ~U(0,1) */
      *mean = atof(argv[6]);
   /* sttdev for RV */
      *stddev = atof(argv[7]);
   /* imax for the pbm file */
      *imax = atoi(argv[8]);
   /* jmax for the pbm file */
      *jmax = atoi(argv[9]);

      if(*flag_UQ < 0)
      {
         Simulation_Stop("Flag for UQ must be positive.");
         return 0;
      }

      if(*flag_distr < 0)
      {
         Simulation_Stop("Flag for pdf must be positive.");
         return 0;
      }

      if(*flag_RV < 0)
      {
         Simulation_Stop("Flag for RV must be positive.");
         return 0;
      }

      if(*flag_prog < 0)
      {
         Simulation_Stop("Flag for program must be positive.");
         return 0;
      }

      if(*npoints < 0)
      {
         Simulation_Stop("The number of points must be positive.");
         return 0;
      }
      if(*imax < 0)
      {
         Simulation_Stop("imax must be a positive integer.");
         return 0;
      }
      if(*jmax < 0)
      {
         Simulation_Stop("jmax must be a positive integer.");
         return 0;
      }
   }
   else
   {
      Simulation_Stop("Please pass the correct number and type of arguments (flag_UQ flag_distr flag_RV npoints mean stddev imax jmax");
         return 0;
      }
      return 1;
   }

   void Simulation_Message(const std::string &txt)
/* produces a stderr text output  */

   {
      int myrank;

      MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
      std::cout << "-MESSAGE- P:" << myrank << " " << txt << std::endl;
      fflush(stdout);
      fflush(stderr);
   }


   void Simulation_Sync(const std::string &txt)
/* produces a stderr textoutput and synchronize all processes */

   {
      int myrank;

      MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
MPI_Barrier(MPI_COMM_WORLD);                             /* synchronize output */  
      std::cout << "-MESSAGE- P:" << myrank << " " << txt << std::endl;
      fflush(stdout);
      fflush(stderr);
      MPI_Barrier(MPI_COMM_WORLD);
   }


   void Simulation_Stop(const std::string &txt)
   {
    int myrank;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
 MPI_Barrier(MPI_COMM_WORLD); /* synchronize output */
    std::cout << "-STOP- P:" << myrank << " " << txt << std::endl;
    fflush(stdout);
    fflush(stderr);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
 }
