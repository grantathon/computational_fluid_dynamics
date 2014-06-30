/*
 * helper.cpp
 *
 *  Created on: Jun 29, 2014
 *      Author: ionut
 */

#include "helper.hpp"

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