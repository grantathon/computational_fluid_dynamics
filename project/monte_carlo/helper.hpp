#ifndef HELPER_HPP_
#define HELPER_HPP_

#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cstring>

void write_data_file(const char* file_name, int* nsampels);

int read_parameters(int argc, char** argv, int *flag_UQ, int* flag_distr, int* flag_RV, int* flag_prog, int *npoints, double* mean, double* stddev, int* imax, int *jmax);

void Simulation_Message(const std::string &txt);
/* produces a stderr text output  */

void Simulation_Sync(const std::string &txt);
/* produces a stderr textoutput and synchronize all processes */

void Simulation_Stop(const std::string &txt);
/* all processes will produce a text output, be synchronized and finished */

#endif /* HELPER_HPP_ */
