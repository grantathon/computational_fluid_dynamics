#include "streaming.h"

void doStreaming(double *collideField, double *streamField, int *flagField, int xlength)
{
	unsigned long x, y, z;
	unsigned long idx0, idx1, idx2, idx3, idx4, idx5, idx6, idx7, idx8, idx10, idx11, idx12,
		idx13, idx14, idx15, idx16, idx17, idx18;

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
					idx0 = (z-1)*(xlength+2)*(xlength+2) + (y-1)*(xlength+2) + x;
					idx1 = (z-1)*(xlength+2)*(xlength+2) + y*(xlength+2) + x - 1;
					idx2 = (z-1)*(xlength+2)*(xlength+2) + y*(xlength+2) + x;
					idx3 = (z-1)*(xlength+2)*(xlength+2) + y*(xlength+2) + x + 1;
					idx4 = (z-1)*(xlength+2)*(xlength+2) + (y+1)*(xlength+2) + x;
					idx5 = z*(xlength+2)*(xlength+2) + (y-1)*(xlength+2) + x - 1;
					idx6 = z*(xlength+2)*(xlength+2) + (y-1)*(xlength+2) + x;
					idx7 = z*(xlength+2)*(xlength+2) + (y-1)*(xlength+2) + x + 1;
					idx8 = z*(xlength+2)*(xlength+2) + y*(xlength+2) + x - 1;
					idx10 = z*(xlength+2)*(xlength+2) + y*(xlength+2) + x + 1;
					idx11 = z*(xlength+2)*(xlength+2) + (y+1)*(xlength+2) + x - 1;
					idx12 = z*(xlength+2)*(xlength+2) + (y+1)*(xlength+2) + x;
					idx13 = z*(xlength+2)*(xlength+2) + (y+1)*(xlength+2) + x + 1;
					idx14 = (z+1)*(xlength+2)*(xlength+2) + (y-1)*(xlength+2) + x;
					idx15 = (z+1)*(xlength+2)*(xlength+2) + y*(xlength+2) + x - 1;
					idx16 = (z+1)*(xlength+2)*(xlength+2) + y*(xlength+2) + x;
					idx17 = (z+1)*(xlength+2)*(xlength+2) + y*(xlength+2) + x + 1;
					idx18 = (z+1)*(xlength+2)*(xlength+2) + (y+1)*(xlength+2) + x;

					streamField[19*idx] 	= collideField[19*idx0+18]; 	/* (0, -1, -1) */
					streamField[19*idx+1] 	= collideField[19*idx1+17]; 	/* (-1, 0, -1) */
					streamField[19*idx+2] 	= collideField[19*idx2+16]; 	/* (0, 0, -1) */
					streamField[19*idx+3] 	= collideField[19*idx3+15]; 	/* (1, 0, -1) */
					streamField[19*idx+4] 	= collideField[19*idx4+14]; 	/* (0, 1, -1) */
					streamField[19*idx+5] 	= collideField[19*idx5+13]; 	/* (-1, -1, 0) */
					streamField[19*idx+6] 	= collideField[19*idx6+12]; 	/* (0, -1, 0) */
					streamField[19*idx+7] 	= collideField[19*idx7+11]; 	/* (1, -1, 0) */
					streamField[19*idx+8] 	= collideField[19*idx8+10]; 	/* (-1, 0, 0) */
					streamField[19*idx+9] 	= collideField[19*idx+9]; 	/* (0, 0, 0) */
					streamField[19*idx+10] = collideField[19*idx10+8]; 	/* (1, 0, 0) */
					streamField[19*idx+11] = collideField[19*idx11+7]; 	/* (-1, 1, 0) */
					streamField[19*idx+12] = collideField[19*idx12+6]; 	/* (0, 1, 0) */
					streamField[19*idx+13] = collideField[19*idx13+5]; 	/* (1, 1, 0) */
					streamField[19*idx+14] = collideField[19*idx14+4]; 	/* (0, -1, 1) */
					streamField[19*idx+15] = collideField[19*idx15+3]; 	/* (-1, 0, 1) */
					streamField[19*idx+16] = collideField[19*idx16+2]; 	/* (0, 0, 1) */
					streamField[19*idx+17] = collideField[19*idx17+1]; 	/* (1, 0, 1) */
					streamField[19*idx+18] = collideField[19*idx18]; 	/* (0, 1, 1) */
				}
			}
		}
	}
}

