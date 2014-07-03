/*
 * helper.hpp
 *
 *  Created on: Jun 29, 2014
 *      Author: ionut
 */

#ifndef HELPER_HPP_
#define HELPER_HPP_

#include <iostream>
#include <mpi.h>
#include <cstdlib>
#include <fstream>

void write_to_file(std::ofstream &outputFile, double reaatach_loc, double reaatch_time);

int read_parameters(int argc, char** argv, int *flag_UQ, int* flag_distr, int* flag_RV, int *npoints, double* mean, double* stddev, int* imax, int *jmax);

void Simulation_Message(const std::string &txt);
/* produces a stderr text output  */

void Simulation_Sync(const std::string &txt);
/* produces a stderr textoutput and synchronize all processes */

void Simulation_Stop(const std::string &txt);
/* all processes will produce a text output, be synchronized and finished */

#endif /* HELPER_HPP_ */