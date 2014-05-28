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
   int i_per_iproc, j_per_jproc, i_rem, j_rem;

	/* Set domain coordinates */
	*omg_i = ((*myrank) % iproc) + 1;
	*omg_j = (((*myrank) - (*omg_i) + 1) / iproc) + 1;

	/*printf("rank: %u, omg_i: %u, omg_j: %u, iproc: %u, jproc: %u\n", *myrank, *omg_i, *omg_j, iproc, jproc);*/

   /* Compute il, ir for each sub-domain*/
   i_per_iproc = (imax / iproc);
   j_per_jproc = (jmax / jproc);
   i_rem = (imax % iproc);
   j_rem = (jmax % jproc);

   /* for il ir*/
   if ((*omg_i) == 1)   /* to rank zero assign the remainder*/
   {
      *il = (((*omg_i) -1) * i_per_iproc) + 1;
      *ir = ((*omg_i) * i_per_iproc) + i_rem;
   }
   else
   {
      *il = (((*omg_i) -1) * i_per_iproc) + i_rem + 1;
      *ir = ((*omg_i) * i_per_iproc) + i_rem;
   }

   /* for jb jt*/
   if ((*omg_j) == 1)
   {
      *jb = (((*omg_j) -1) * j_per_jproc) + 1;
      *jt = ((*omg_j) * j_per_jproc) + j_rem;
   }
   else
   {
      *jb = (((*omg_j) -1) * j_per_jproc) + j_rem + 1;
      *jt = ((*omg_j) * j_per_jproc) + j_rem;
   }


   /* Assign rank of neighbour to each sub-domain: rank_l, rank_r, rank_b, rank_t*/
   /* Left boundary*/
   if ((*il) == 1)
   {
      *rank_l = MPI_PROC_NULL;
   }
   else
   {
      *rank_l = (*myrank) - 1;
   }

   /* Right boundary*/
   if ((*ir) == imax)
   {
      *rank_r = MPI_PROC_NULL;
   }
   else
   {
      *rank_r = (*myrank) + 1;
   }

   /* Bottom boundary*/
   if ((*jb) == 1)
   {
      *rank_b = MPI_PROC_NULL;
   }
   else
   {
      *rank_b = (*myrank) - iproc;
   }

   /* Top boundary*/
   if ((*jt) == jmax)
   {
      *rank_t = MPI_PROC_NULL;
   }
   else
   {
      *rank_t = (*myrank) + iproc;  
   }

   printf("rank: %u, omg_i: %u, omg_j: %u, il: %u\t, ir: %u, jb: %u, jt: %u\t, rank_l: %u, rank_r: %u\t, rank_b: %u, rank_t: %u\n",
            *myrank, *omg_i, *omg_j, *il, *ir, *jb, *jt, *rank_l, *rank_r, *rank_b, *rank_t);

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
