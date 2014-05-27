#include "parallel.h"


void Program_Message(char *txt)
/* produces a stderr text output  */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
}


void Programm_Sync(char *txt)
/* produces a stderr textoutput and synchronize all processes */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                             /* synchronize output */  
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
}


void Programm_Stop(char *txt)
/* all processes will produce a text output, be synchronized and finished */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                           /* synchronize output */
   fprintf(stderr,"-STOP- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(1);
}

void init_parallel(int iproc,
               int jproc,
               int imax,
               int jmax,
               int *myrank,
               int *il,
               int *ir,
               int *jb,
               int *jt,
               int *rank_l,
               int *rank_r,
               int *rank_b,
               int *rank_t,
               int *omg_i,
               int *omg_j,
               int num_proc)
{
	/* Set domain coordinates */
	*omg_i = ((*myrank) % iproc) + 1;
	*omg_j = (((*myrank) - (*omg_i) + 1) / iproc) + 1;

	printf("omg_i = %u, omg_j = %u\n", *omg_i, *omg_j);
}

void pressure_comm(double **P,
            int il,
            int ir,
            int jb,
            int jt,
            int rank_l,
            int rank_r,
            int rank_b,
            int rank_t,
            double *bufSend,
            double *bufRecv, 
            MPI_Status *status, 
            int chunk)
{
   /* TODO: */
}

void uv_comm(double **U,
         double **V,
         int il,
         int ir,
         int jb,
         int jt,
         int rank_l,
         int rank_r,
         int rank_b,
         int rank_t,
         double *bufSend,
         double *bufRecv,
         MPI_Status *status,
         int chunk)
{
   /* TODO: */
}
