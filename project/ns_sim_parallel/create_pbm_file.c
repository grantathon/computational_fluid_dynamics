#include "string.h"
#include <stdio.h>
#include <stdlib.h>


void make_pbm_file(const char *filename, int imax, int jmax)
{

  FILE *input = NULL;
  int i, j;

  /* Initialize reading/writing of output file */
  input = fopen(filename, "w");
  if(input == 0)
  {
    char szBuff[80];
    sprintf(szBuff, "Can not append or create file %s", filename);
  }
/*P1
# This describes a step (block) that will be located within our (dim_x+2)*(dim_y+2) domain
imax jmax
*/

  fprintf(input, "P1\n");
  fprintf(input, "# This describes a step (block) that will be located within our (dim_x+2)*(dim_y+2) domain\n");
  fprintf(input, "%i %i\n", imax, jmax);
  fprintf(input, "\n");

  for(j = 1; j < jmax+1; j++) {
    for(i = 1; i < imax+1; i++) {
      if ((j > jmax/2) && (i <= jmax/2))
      {
        fprintf(input, "%i ", 0);
      }
      else
      {
        fprintf(input, "%i ", 1);
      }

      if ((i % imax) == 0)
      {
        fprintf(input, "\n");
      }
    }
  }

  fclose(input);
}

int main(int argc, char *argv[])
{
	/* Input file with user parameters */
	const char *problem = "flow_over_a_step";
	char *problemDataFile_pbm = (char*)malloc(strlen(problem) + 4);

	/* Make pbm file  */
	strcpy(problemDataFile_pbm, problem);
	strcat(problemDataFile_pbm, ".pbm");

	int imax = 0;
	int jmax = 0;

	/* Read passed Reynolds number and unique Monte Carlo ID */
	if(argc == 3)
	{
		imax = atoi(argv[1]);
		jmax = atoi(argv[2]);

		if(imax <= 0)
		{
			printf("imax is less than zero! Setting to default imax=100\n");
			imax = 100;
		}

		if(jmax < 0)
		{
			printf("jmax is less than zero! Setting to default jmax=20\n");
			jmax = 20;
		}

		printf("imax: %i\n", imax);
		printf("jmax: %i\n", jmax);
	}
	else
	{
		printf("not enough input arguments! Exiting.\n");
		return 0;
	}
	


	/* Read parameters from DAT file, store locally, and check for potential error */
	make_pbm_file((const char *)problemDataFile_pbm, imax, jmax);
	printf("PBM file successfully created!\n");
	return 1;
}
