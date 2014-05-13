#include "helper.h"
#include "init.h"
#include "ns_definitions.h"
#include <string.h>

int read_parameters(const char * szFileName,
					double *Re,                /* reynolds number   */
                    double *UI,                /* velocity x-direction */
                    double *VI,                /* velocity y-direction */
                    double *PI,                /* pressure */
                    double *GX,                /* gravitation x-direction */
                    double *GY,                /* gravitation y-direction */
                    double *t_end,             /* end time */
                    double *xlength,           /* length of the domain x-dir.*/
                    double *ylength,           /* length of the domain y-dir.*/
                    double *dt,                /* time step */
                    double *dx,                /* length of a cell x-dir. */
                    double *dy,                /* length of a cell y-dir. */
                    int  *imax,                /* number of cells x-direction*/
                    int  *jmax,                /* number of cells y-direction*/
                    double *alpha,             /* uppwind differencing factor*/
                    double *omg,               /* relaxation factor */
                    double *tau,               /* safety factor for time step*/
                    int  *itermax,             /* max. number of iterations  */
		                               	   	   /* for pressure per time step */
                    double *eps,               /* accuracy bound for pressure*/
                    double *dt_value,			/* time for output */
                    int *wl,			/* domain boundary conditions */
                    int *wr,			/* for right, left, top and */
                    int *wt,			/* bottom surfaces. */
                    int *wb
)
{
	/* Read domain boundary conditions*/
	READ_INT   ( szFileName, *wl );
	READ_INT   ( szFileName, *wr );
	READ_INT   ( szFileName, *wt );
	READ_INT   ( szFileName, *wb );

	/* Read all parameters */
	READ_DOUBLE( szFileName, *xlength );
	READ_DOUBLE( szFileName, *ylength );

	READ_DOUBLE( szFileName, *Re    );
	READ_DOUBLE( szFileName, *t_end );
	READ_DOUBLE( szFileName, *dt    );

	READ_INT   ( szFileName, *imax );
	READ_INT   ( szFileName, *jmax );

	READ_DOUBLE( szFileName, *omg   );
	READ_DOUBLE( szFileName, *eps   );
	READ_DOUBLE( szFileName, *tau   );
	READ_DOUBLE( szFileName, *alpha );

	READ_INT   ( szFileName, *itermax );
	READ_DOUBLE( szFileName, *dt_value );

	READ_DOUBLE( szFileName, *UI );
	READ_DOUBLE( szFileName, *VI );
	READ_DOUBLE( szFileName, *GX );
	READ_DOUBLE( szFileName, *GY );
	READ_DOUBLE( szFileName, *PI );

	*dx = *xlength / (double)(*imax);
	*dy = *ylength / (double)(*jmax);

	return 1;
}

void read_special_BC(const char * szFileName,
		double *UI,	/* Read inlet velocity BC*/
		double *VI,
		double *delta_p
)
{
	READ_DOUBLE( szFileName, *UI );
	READ_DOUBLE( szFileName, *VI );

	if (strcmp(szFileName, "plane_shear_flow.dat") == 0)
	{
		READ_DOUBLE( szFileName, *delta_p );
	}

}

/*
 * Initialize the U, V, and P grids to what UI, VU, and PI
 * are set to, respectively.
 */
void init_uvp(
	double UI,
	double VI,
	double PI,
	int imax,
	int jmax,
	double ***U,
	double ***V,
	double ***P
)
{
	/* Horizontal velocity */
	*U = matrix(0, imax+1, 0, jmax+1);
	init_matrix(*U, 0, imax+1, 0, jmax+1, UI);

	/* Vertical velocity */
	*V = matrix(0, imax+1, 0, jmax+1);
	init_matrix(*V, 0, imax+1, 0, jmax+1, VI);

	/* Pressure */
	*P = matrix(0, imax+1, 0, jmax+1);
	init_matrix(*P, 0, imax+1, 0, jmax+1, PI);
}

void init_flag(const char *problem, int imax, int jmax, int ***Flag)
{
	/* initialize Flag as a matrix of integers */
	*Flag = imatrix(0, imax + 1, 0, jmax + 1);

	/* TODO: Define any objects within our domain here */

	/* start with problem b) */
	if(strcmp(problem, "plane_shear_flow") == 0)
	{
		init_imatrix(*Flag, 0, imax + 1, 0, jmax + 1, C_F);
	}
	else if (strcmp(problem, "Karman_vortex") == 0)
	{

	}
	else if(strcmp(problem, "flow_over_a_step") == 0)
	{
		
	}
}



