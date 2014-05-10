#include "boundary.h"
#include "math.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include <stdio.h>

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength)
{
	unsigned long idx = 0;
	unsigned long temp_idx  = 0;
	unsigned long idx_inv  = 0;
	unsigned long x, y, z, i, x_inv, y_inv, z_inv;
	double density = 0.0;
	double cuDotProd = 0.0;
	const double cs = 1.0/sqrt(3.0);

	/* for z walls */
	for(z = 0; z <= (xlength + 1); z += (xlength + 1))
	{
		for(y = 1; y <= xlength; y++)
		{
			for(x = 1; x <= xlength; x++)
			{
				temp_idx = z*(xlength+2)*(xlength+2) + y*(xlength+2) + x;

				if(z == 0)
				{
					idx = (z+1)*(xlength+2)*(xlength+2) + y*(xlength+2) + x;
				}
				else
				{
					idx = (z-1)*(xlength+2)*(xlength+2) + y*(xlength+2) + x;
				}

				for(i = 0; i < 19; i++)
				{
					if ((z == 0 && LATTICEVELOCITIES[i][2] == -1) || (z == (xlength + 1) && LATTICEVELOCITIES[i][2] == 1))
					{
						x_inv = x + LATTICEVELOCITIES[i][0];
						y_inv = y + LATTICEVELOCITIES[i][1];

						idx_inv = z*(xlength+2)*(xlength+2) + y_inv*(xlength+2) + x_inv;

						if(flagField[temp_idx] == NO_SLIP)
						{
							collideField[19*idx_inv + (18 - i)] = collideField[19*idx + i];
						}
						else if(flagField[temp_idx] == MOVING_WALL)
						{
							computeDensity(&collideField[19*idx], &density);
							cuDotProd = LATTICEVELOCITIES[i][0]*wallVelocity[0] + LATTICEVELOCITIES[i][1]*wallVelocity[1] +
									LATTICEVELOCITIES[i][2]*wallVelocity[2];
							collideField[19*idx_inv + (18 - i)] = collideField[19*idx + i] +
									((2 * LATTICEWEIGHTS[i] * density * cuDotProd) / (cs*cs));
						}
					}
				}
			}
		}
	}

	/* for y walls */
	for(y = 0; y <= (xlength + 1); y += (xlength + 1))
	{
		for(z = 1; z <= xlength; z++)
		{
			for(x = 1; x <= xlength; x++)
			{
				temp_idx = z*(xlength+2)*(xlength+2) + y*(xlength+2) + x;

				if(y == 0)
				{
					idx = z*(xlength+2)*(xlength+2) + (y+1)*(xlength+2) + x;
				}
				else
				{
					idx = z*(xlength+2)*(xlength+2) + (y-1)*(xlength+2) + x;
				}

				for(i = 0; i < 19; i++)
				{
					if ((y == 0 && LATTICEVELOCITIES[i][1] == -1) || (y == (xlength + 1) && LATTICEVELOCITIES[i][1] == 1))
					{
						x_inv = x + LATTICEVELOCITIES[i][0];
						z_inv = z + LATTICEVELOCITIES[i][2];

						idx_inv = z_inv*(xlength+2)*(xlength+2) + y*(xlength+2) + x_inv;

						if(flagField[temp_idx] == NO_SLIP)
						{
							collideField[19*idx_inv + (18 - i)] = collideField[19*idx + i];
						}
						else if(flagField[temp_idx] == MOVING_WALL)
						{
							computeDensity(&collideField[19*idx], &density);
							cuDotProd = LATTICEVELOCITIES[i][0]*wallVelocity[0] + LATTICEVELOCITIES[i][1]*wallVelocity[1] +
									LATTICEVELOCITIES[i][2]*wallVelocity[2];
							collideField[19*idx_inv + (18 - i)] = collideField[19*idx + i] +
									((2 * LATTICEWEIGHTS[i] * density * cuDotProd) / (cs*cs));
						}
					}
				}
			}
		}
	}

	/* for x walls */
	for(x = 0; x <= (xlength + 1); x += (xlength + 1))
	{
		for(z = 1; z <= xlength; z++)
		{
			for(y = 1; y <= xlength; y++)
			{
				temp_idx = z*(xlength+2)*(xlength+2) + y*(xlength+2) + x;

				if(x == 0)
				{
					idx = z*(xlength+2)*(xlength+2) + y*(xlength+2) + x + 1;
				}
				else
				{
					idx = z*(xlength+2)*(xlength+2) + y*(xlength+2) + x - 1;
				}

				for(i = 0; i < 19; i++)
				{
					if ((x == 0 && LATTICEVELOCITIES[i][0] == -1) || (x == (xlength + 1) && LATTICEVELOCITIES[i][0] == 1))
					{
						y_inv = y + LATTICEVELOCITIES[i][1];
						z_inv = z + LATTICEVELOCITIES[i][2];

						idx_inv = z_inv*(xlength+2)*(xlength+2) + y_inv*(xlength+2) + x;

						if(flagField[temp_idx] == NO_SLIP)
						{
							collideField[19*idx_inv + (18 - i)] = collideField[19*idx + i];
						}
						else if(flagField[temp_idx] == MOVING_WALL)
						{
							computeDensity(&collideField[19*idx], &density);
							cuDotProd = LATTICEVELOCITIES[i][0]*wallVelocity[0] + LATTICEVELOCITIES[i][1]*wallVelocity[1] +
									LATTICEVELOCITIES[i][2]*wallVelocity[2];
							collideField[19*idx_inv + (18 - i)] = collideField[19*idx + i] +
									((2 * LATTICEWEIGHTS[i] * density * cuDotProd) / (cs*cs));
						}
					}
				}
			}
		}
	}
}

