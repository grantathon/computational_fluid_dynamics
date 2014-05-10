#ifndef _MAIN_C_
#define _MAIN_C_

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include <stdio.h>

int main(int argc, char *argv[]){
	double *collideField = NULL;
	double *streamField = NULL;
	int *flagField = NULL;
	int xlength = 0;
	double tau = 0.0;
	double velocityWall[3] = {0};
	int timesteps = 0;
	int timestepsPerPlotting = 0;
	int t = 0;
	int readParamError = 0;
	const char *szProblem = "data_CFD_Assignment_02";

	readParamError = readParameters(&xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, argc, argv);
	if(readParamError != 0)
	{
		printf("Please provide only one argument to the program; the path to the input file.");
		return -1;
	}

	/* Initialize pointers according to the D3Q19 discretization scheme */
	collideField 	= malloc(19*(xlength+2)*(xlength+2)*(xlength+2)*sizeof(*collideField));
	streamField 	= malloc(19*(xlength+2)*(xlength+2)*(xlength+2)*sizeof(*streamField));
	flagField 		= malloc((xlength+2)*(xlength+2)*(xlength+2)*sizeof(*flagField));

	initialiseFields(collideField, streamField, flagField, xlength);

	for(t = 0; t < timesteps; t++)
	{
		double *swap = NULL;

		doStreaming(collideField, streamField, flagField, xlength);

		swap = collideField;
		collideField = streamField;
		streamField = swap;

		doCollision(collideField, flagField, &tau, xlength);

		treatBoundary(collideField, flagField, velocityWall, xlength);

		if (t % timestepsPerPlotting == 0)
		{
			writeVtkOutput(collideField, flagField, szProblem, t, xlength);
		}
	}

	/* Free heap memory */
	free(collideField);
	free(streamField);
	free(flagField);

	return 0;
}

#endif
