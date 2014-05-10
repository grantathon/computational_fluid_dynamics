#include "computeCellValues.h"
#include "LBDefinitions.h"
#include "math.h"
#include <stdio.h>

void computeDensity(const double *const currentCell, double *density)
{
	*density = 	currentCell[0] + currentCell[1] + currentCell[2] + currentCell[3] + currentCell[4] +
				currentCell[5] + currentCell[6] + currentCell[7] + currentCell[8] + currentCell[9] +
				currentCell[10] + currentCell[11] + currentCell[12] + currentCell[13] + currentCell[14] +
				currentCell[15] + currentCell[16] + currentCell[17] + currentCell[18];
}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity)
{
	int i;

	if((*density) != 0)
	{
		for(i = 0; i < 3; i++)
		{
			velocity[i] = 	(	currentCell[0]*LATTICEVELOCITIES[0][i] + currentCell[1]*LATTICEVELOCITIES[1][i] +
								currentCell[2]*LATTICEVELOCITIES[2][i] + currentCell[3]*LATTICEVELOCITIES[3][i] +
								currentCell[4]*LATTICEVELOCITIES[4][i] + currentCell[5]*LATTICEVELOCITIES[5][i] +
								currentCell[6]*LATTICEVELOCITIES[6][i] + currentCell[7]*LATTICEVELOCITIES[7][i] +
								currentCell[8]*LATTICEVELOCITIES[8][i] + currentCell[9]*LATTICEVELOCITIES[9][i] +
								currentCell[10]*LATTICEVELOCITIES[10][i] + currentCell[11]*LATTICEVELOCITIES[11][i] +
								currentCell[12]*LATTICEVELOCITIES[12][i] + currentCell[13]*LATTICEVELOCITIES[13][i] +
								currentCell[14]*LATTICEVELOCITIES[14][i] + currentCell[15]*LATTICEVELOCITIES[15][i] +
								currentCell[16]*LATTICEVELOCITIES[16][i] + currentCell[17]*LATTICEVELOCITIES[17][i] +
								currentCell[18]*LATTICEVELOCITIES[18][i]	) / (*density);
		}
	}
}

void computeFeq(const double * const density, const double * const velocity, double *feq)
{
	int i;
	const double cs = 1.0/sqrt(3.0);
	const double uuDotProd = velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2];

	for(i = 0; i < 19; i++)
	{
		const double cuDotProd = LATTICEVELOCITIES[i][0]*velocity[0] + LATTICEVELOCITIES[i][1]*velocity[1] +
				LATTICEVELOCITIES[i][2]*velocity[2];

		feq[i] = LATTICEWEIGHTS[i]*(*density)*(1 + (cuDotProd/(cs*cs)) +
				((cuDotProd*cuDotProd)/(2*cs*cs*cs*cs)) -
				((uuDotProd)/(2*cs*cs)));
	}
}

