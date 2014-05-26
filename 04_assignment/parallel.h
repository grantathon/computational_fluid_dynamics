#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>

void Program_Message(char *txt);
/* produces a stderr text output  */



void Programm_Sync(char *txt);
/* produces a stderr textoutput and synchronize all processes */



void Programm_Stop(char *txt);
/* all processes will produce a text output, be synchronized and finished */

void init parallel (int iproc,
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
					int num proc);


void pressure comm(double **P,
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
				MPI Status *status, 
				int chunk);


void uv comm(double **U,
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
			MPI Status *status,
			int chunk);