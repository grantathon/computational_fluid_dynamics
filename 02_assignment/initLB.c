#include "initLB.h"
#include "helper.h"
#include "LBDefinitions.h"

int readParameters(int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[])
{
	const char * filename = argv[1];
	double velocityWallX;
	double velocityWallY;
	double velocityWallZ;

	/* Make sure arg input count, not including the program call, is 1 */
	if(argc != 2)
	{
		return -1;
	}

	/* Read all parameters */
	READ_INT(filename, *xlength);
	READ_DOUBLE(filename, *tau);
	READ_DOUBLE(filename, velocityWallX);
	READ_DOUBLE(filename, velocityWallY);
	READ_DOUBLE(filename, velocityWallZ);
	READ_INT(filename, *timesteps);
	READ_INT(filename, *timestepsPerPlotting);

	/* Set wall velocity */
	velocityWall[0] = velocityWallX;
	velocityWall[1] = velocityWallY;
	velocityWall[2] = velocityWallZ;

	return 0;
}

void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength)
{
	unsigned long x, y, z;

	/* Loop through fluid cells, not boundary nodes */
	for(z = 0; z <= (xlength+1); z++)
	{
		for(y = 0; y <= (xlength+1); y++)
		{
			for(x = 0; x <= (xlength+1); x++)
			{
				const unsigned long idx = z*(xlength+2)*(xlength+2) + y*(xlength+2) + x;

				collideField[19*idx] 	= LATTICEWEIGHTS[0]; 	/* (0, -1, -1) */
				collideField[19*idx+1] 	= LATTICEWEIGHTS[1]; 	/* (-1, 0, -1) */
				collideField[19*idx+2] 	= LATTICEWEIGHTS[2]; 	/* (0, 0, -1) */
				collideField[19*idx+3] 	= LATTICEWEIGHTS[3]; 	/* (1, 0, -1) */
				collideField[19*idx+4] 	= LATTICEWEIGHTS[4]; 	/* (0, 1, -1) */
				collideField[19*idx+5] 	= LATTICEWEIGHTS[5]; 	/* (-1, -1, 0) */
				collideField[19*idx+6] 	= LATTICEWEIGHTS[6]; 	/* (0, -1, 0) */
				collideField[19*idx+7] 	= LATTICEWEIGHTS[7]; 	/* (1, -1, 0) */
				collideField[19*idx+8] 	= LATTICEWEIGHTS[8]; 	/* (-1, 0, 0) */
				collideField[19*idx+9] 	= LATTICEWEIGHTS[9]; 	/* (0, 0, 0) */
				collideField[19*idx+10] = LATTICEWEIGHTS[10]; 	/* (1, 0, 0) */
				collideField[19*idx+11] = LATTICEWEIGHTS[11]; 	/* (-1, 1, 0) */
				collideField[19*idx+12] = LATTICEWEIGHTS[12]; 	/* (0, 1, 0) */
				collideField[19*idx+13] = LATTICEWEIGHTS[13]; 	/* (1, 1, 0) */
				collideField[19*idx+14] = LATTICEWEIGHTS[14]; 	/* (0, -1, 1) */
				collideField[19*idx+15] = LATTICEWEIGHTS[15]; 	/* (-1, 0, 1) */
				collideField[19*idx+16] = LATTICEWEIGHTS[16]; 	/* (0, 0, 1) */
				collideField[19*idx+17] = LATTICEWEIGHTS[17]; 	/* (1, 0, 1) */
				collideField[19*idx+18] = LATTICEWEIGHTS[18]; 	/* (0, 1, 1) */

				streamField[19*idx] 	= LATTICEWEIGHTS[0]; 	/* (0, -1, -1) */
				streamField[19*idx+1] 	= LATTICEWEIGHTS[1]; 	/* (-1, 0, -1) */
				streamField[19*idx+2] 	= LATTICEWEIGHTS[2]; 	/* (0, 0, -1) */
				streamField[19*idx+3] 	= LATTICEWEIGHTS[3]; 	/* (1, 0, -1) */
				streamField[19*idx+4] 	= LATTICEWEIGHTS[4]; 	/* (0, 1, -1) */
				streamField[19*idx+5] 	= LATTICEWEIGHTS[5]; 	/* (-1, -1, 0) */
				streamField[19*idx+6] 	= LATTICEWEIGHTS[6]; 	/* (0, -1, 0) */
				streamField[19*idx+7] 	= LATTICEWEIGHTS[7]; 	/* (1, -1, 0) */
				streamField[19*idx+8] 	= LATTICEWEIGHTS[8]; 	/* (-1, 0, 0) */
				streamField[19*idx+9] 	= LATTICEWEIGHTS[9]; 	/* (0, 0, 0) */
				streamField[19*idx+10] 	= LATTICEWEIGHTS[10]; 	/* (1, 0, 0) */
				streamField[19*idx+11]	= LATTICEWEIGHTS[11]; 	/* (-1, 1, 0) */
				streamField[19*idx+12] 	= LATTICEWEIGHTS[12]; 	/* (0, 1, 0) */
				streamField[19*idx+13] 	= LATTICEWEIGHTS[13]; 	/* (1, 1, 0) */
				streamField[19*idx+14] 	= LATTICEWEIGHTS[14]; 	/* (0, -1, 1) */
				streamField[19*idx+15] 	= LATTICEWEIGHTS[15]; 	/* (-1, 0, 1) */
				streamField[19*idx+16] 	= LATTICEWEIGHTS[16]; 	/* (0, 0, 1) */
				streamField[19*idx+17] 	= LATTICEWEIGHTS[17]; 	/* (1, 0, 1) */
				streamField[19*idx+18] 	= LATTICEWEIGHTS[18]; 	/* (0, 1, 1) */

				if((x < xlength+1) && (x > 0) && (y < xlength+1) && (y > 0) && (z < xlength+1) && (z > 0))
				{
					flagField[idx] = FLUID;
				}
				else if((x == 0) || (y == 0) || (z == 0) || (x == xlength+1) || (y == xlength+1))
				{
					flagField[idx] = NO_SLIP;
				}
				else /* z == xlength+1 */
				{
					flagField[idx] = MOVING_WALL;
				}
			}
		}
	}
}

