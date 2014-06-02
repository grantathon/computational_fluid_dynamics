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

	/* Set sub-domain coordinates */
	*omg_i = ((*myrank) % iproc) + 1;
	*omg_j = (((*myrank) - (*omg_i) + 1) / iproc) + 1;

	/* Compute il, ir for each sub-domain*/
	i_per_iproc = (imax / iproc);
	j_per_jproc = (jmax / jproc);
	i_rem = (imax % iproc);
	j_rem = (jmax % jproc);

	/* for il ir*/
	if ((*omg_i) == 1)   /* to rank zero assign the remainder*/
	{
	  *il = (((*omg_i) -1) * i_per_iproc) + 1;
	}
	else
	{
	  *il = (((*omg_i) -1) * i_per_iproc) + i_rem + 1;
	}
	*ir = ((*omg_i) * i_per_iproc) + i_rem;

	/* for jb jt*/
	if ((*omg_j) == 1)
	{
	  *jb = (((*omg_j) -1) * j_per_jproc) + 1;
	}
	else
	{
	  *jb = (((*omg_j) -1) * j_per_jproc) + j_rem + 1;
	}
	*jt = ((*omg_j) * j_per_jproc) + j_rem;


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
	printf("Entered pressure_comm() \n");
	int i, j, size;

	/* Send to the left, receive from the right */
	/* Send to the right, receive from the left */
	if(rank_l != MPI_PROC_NULL || rank_r != MPI_PROC_NULL)
	{
		size = jt - jb + 1;
		printf("size = %d\n", size);

		if(rank_l != MPI_PROC_NULL && rank_r != MPI_PROC_NULL)  /* Perform both left-right transfers */
		{
			/* Need two buffers for data transfer */
			bufSend = malloc(size*sizeof(double));
			bufRecv = malloc(size*sizeof(double));

			/* Copy left values to send */
			for(j = 1; j <= size; j++)
			{
				bufSend[j - 1] = P[1][j];
			}

			/* Send left values, receive right values */
			printf("Before MPI_Sendrecv for rank_l=%u & rank_r=%u\n", rank_l, rank_r);
			MPI_Sendrecv(bufSend, size, MPI_DOUBLE, rank_l, 1, bufRecv, size, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, status);
			printf("After MPI_Sendrecv for rank_l=%u & rank_r=%u\n", rank_l, rank_r);

			/* Copy received right values */
			for(j = 1; j <= size; j++)
			{
				P[size + 1][j] = bufRecv[j - 1];
			}

			/* Copy right values to send */
			for(j = 1; j <= size; j++)
			{
				bufSend[j - 1] = P[size][j];
			}

			/* Send right values, receive left values */
			printf("Before MPI_Sendrecv for rank_l=%u & rank_r=%u\n", rank_l, rank_r);
			MPI_Sendrecv(bufSend, size, MPI_DOUBLE, rank_r, 1, bufRecv, size, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD, status);
			printf("After MPI_Sendrecv for rank_l=%u & rank_r=%u\n", rank_l, rank_r);

			/* Copy received left values */
			for(j = 1; j <= size; j++)
			{
				P[0][j] = bufRecv[j - 1];
			}

			free(bufSend);
			free(bufRecv);
		}
		else if(rank_l == MPI_PROC_NULL && rank_r != MPI_PROC_NULL)  /* Only receive/send data from/to right */
		{
			/* Need only one buffer for data transfer */
			bufSend = malloc(size*sizeof(double));

			/* Receive right values */
			printf("Before MPI_Recv from rank_r=%u\n", rank_r);
			MPI_Recv(bufSend, size, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, status);
			printf("After MPI_Recv from rank_r=%u\n", rank_r);

			/* Copy received right values */
			for(j = 1; j <= size; j++)
			{
				P[size + 1][j] = bufSend[j - 1];
			}

			/* Copy right values to send */
			for(j = 1; j <= size; j++)
			{
				bufSend[j - 1] = P[size][j];
			}

			/* Send right values */
			printf("Before MPI_Send to rank_r=%u\n", rank_r);
			MPI_Send(bufSend, size, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD);
			printf("After MPI_Send to rank_r=%u\n", rank_r);

			free(bufSend);
		}
		else if(rank_l != MPI_PROC_NULL && rank_r == MPI_PROC_NULL)  /* Only send/receive data to/from left */
		{
			/* Need only one buffer for data transfer */
			bufSend = malloc(size*sizeof(double));

			/* Copy left values to send */
			for(j = 1; j <= size; j++)
			{
				bufSend[j - 1] = P[1][j];
			}

			/* Send left values */
			printf("Before MPI_Send to rank_l=%u\n", rank_l);
			MPI_Send(bufSend, size, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD);
			printf("After MPI_Send to rank_l=%u\n", rank_l);

			/* Receive left values */
			printf("Before MPI_Recv from rank_l=%u\n", rank_l);
			MPI_Recv(bufSend, size, MPI_DOUBLE, rank_l, 1, MPI_COMM_WORLD, status);
			printf("After MPI_Recv from rank_l=%u\n", rank_l);

			/* Copy received left values */
			for(j = 1; j <= size; j++)
			{
				P[0][j] = bufSend[j - 1];
			}

			free(bufSend);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);  /* Wait for all processes */

	/* Send to the top, receive from the bottom */
	/* Send to the bottom, receive from the top */
	if(rank_t != MPI_PROC_NULL || rank_b != MPI_PROC_NULL)
	{
		size = ir - il + 1;
		printf("size = %d\n", size);

		if(rank_t != MPI_PROC_NULL && rank_b != MPI_PROC_NULL)  /* Perform both top-bottom transfers */
		{
			/* Need two buffers for data transfer */
			bufSend = malloc(size*sizeof(double));
			bufRecv = malloc(size*sizeof(double));

			/* Copy top values to send */
			for(i = 1; i <= size; i++)
			{
				bufSend[i - 1] = P[i][size];
			}

			/* Send top values, receive bottom values */
			printf("Before MPI_Sendrecv for rank_t=%u & rank_b=%u\n", rank_t, rank_b);
			MPI_Sendrecv(bufSend, size, MPI_DOUBLE, rank_t, 1, bufRecv, size, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD, status);
			printf("After MPI_Sendrecv for rank_t=%u & rank_b=%u\n", rank_t, rank_b);

			/* Copy received bottom values */
			for(i = 1; i <= size; i++)
			{
				P[i][0] = bufRecv[i - 1];
			}

			/* Copy bottom values to send */
			for(i = 1; i <= size; i++)
			{
				bufSend[i - 1] = P[i][1];
			}

			/* Send bottom values, receive top values */
			printf("Before MPI_Sendrecv for rank_t=%u & rank_b=%u\n", rank_t, rank_b);
			MPI_Sendrecv(bufSend, size, MPI_DOUBLE, rank_b, 1, bufRecv, size, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, status);
			printf("After MPI_Sendrecv for rank_t=%u & rank_b=%u\n", rank_t, rank_b);

			/* Copy received top values */
			for(i = 1; i <= size; i++)
			{
				P[i][size + 1] = bufRecv[i - 1];
			}

			free(bufSend);
			free(bufRecv);
		}
		else if(rank_t == MPI_PROC_NULL && rank_b != MPI_PROC_NULL)  /* Only receive/send data from/to bottom */
		{
			/* Need only one buffer for data transfer */
			bufSend = malloc(size*sizeof(double));

			/* Receive bottom values */
			printf("Before MPI_Recv from rank_b=%u\n", rank_b);
			MPI_Recv(bufSend, size, MPI_DOUBLE, rank_r, 1, MPI_COMM_WORLD, status);
			printf("After MPI_Recv from rank_b=%u\n", rank_b);

			/* Copy received bottom values */
			for(i = 1; i <= size; i++)
			{
				P[i][0] = bufSend[i - 1];
			}

			/* Copy bottom values to send */
			for(i = 1; i <= size; i++)
			{
				bufSend[i - 1] = P[i][1];
			}

			/* Send bottom values */
			printf("Before MPI_Send to rank_b=%u\n", rank_b);
			MPI_Send(bufSend, size, MPI_DOUBLE, rank_b, 1, MPI_COMM_WORLD);
			printf("After MPI_Send to rank_b=%u\n", rank_b);

			free(bufSend);
		}
		else if(rank_t != MPI_PROC_NULL && rank_b == MPI_PROC_NULL)  /* Only send/receive data to/from top */
		{
			/* Need only one buffer for data transfer */
			bufSend = malloc(size*sizeof(double));

			/* Copy top values to send */
			for(i = 1; i <= size; i++)
			{
				bufSend[i - 1] = P[i][size];
			}

			/* Send top values */
			printf("Before MPI_Send to rank_t=%u\n", rank_t);
			MPI_Send(bufSend, size, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD);
			printf("After MPI_Send to rank_t\n");

			/* Receive top values */
			printf("Before MPI_Recv from rank_t=%u\n", rank_t);
			MPI_Recv(bufSend, size, MPI_DOUBLE, rank_t, 1, MPI_COMM_WORLD, status);
			printf("After MPI_Recv from rank_t=%u\n", rank_t);

			/* Copy received top values */
			for(i = 1; i <= size; i++)
			{
				P[i][size + 1] = bufSend[i - 1];
			}

			free(bufSend);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);  /* Wait for all processes to finish */
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
	int i, j, size;

	/* Send to the left, receive from the right */
	if(rank_l != MPI_PROC_NULL)
	{
		/* Allocate memory for buffers */
		size = jt - jb + 1;
		bufSend = malloc(size*sizeof(double));
		bufRecv = malloc(size*sizeof(double));

		/* Copy left values to send */
		for(j = 1; j <= size; j++)
		{
			bufSend[j - 1] = P[1][j];
		}

		MPI_Sendrecv(bufSend, size, MPI_DOUBLE, rank_l, 0, bufRecv, size, MPI_DOUBLE, (rank_l + 1), 0, MPI_COMM_WORLD, status);

		/* Copy received right values */
		for(j = 1; j <= size; j++)
		{
			P[size + 1][j] = bufRecv[j - 1];
		}

		free(bufSend);
		free(bufRecv);
	}

}
