#include "collision.h"
#include <stdlib.h>

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq)
{
	int i;

	/* Perform collision timestep for cell distribution */
	for(i = 0; i < 19; i++)
	{
		currentCell[i] = currentCell[i] - (currentCell[i] - feq[i])/(*tau);
	}
}

void doCollision(double *collideField, int *flagField, const double * const tau, int xlength)
{
	unsigned long x, y, z;
	double density = 0.0;
	double velocity[3] = {0};
	double feq[19] = {0};

	/* Loop through fluid cells, not boundary nodes */
	for(z = 0; z <= (xlength+1); z++)
	{
		for(y = 0; y <= (xlength+1); y++)
		{
			for(x = 0; x <= (xlength+1); x++)
			{
				const unsigned long idx = z*(xlength+2)*(xlength+2) + y*(xlength+2) + x;

				if(flagField[idx] == FLUID)
				{
					computeDensity(&(collideField[19*idx]), &density);
					computeVelocity(&(collideField[19*idx]), &density, velocity);
					computeFeq(&density, velocity, feq);
					computePostCollisionDistributions(&(collideField[19*idx]), tau, feq);
				}
			}
		}
	}
}

