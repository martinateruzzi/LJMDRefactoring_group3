/* 
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "prototypes.h"
/* helper function: apply minimum image convention */
static double pbc(double x, const double boxby2)
{
    while (x >  boxby2) x -= 2.0*boxby2;
    while (x < -boxby2) x += 2.0*boxby2;
    return x;
}

/* compute forces */
void force(mdsys_t *sys) 
{
    	double epot = 0.0;
	double sigma6, c6, c12, rcsq;
	
	sigma6 = sys->sigma * sys->sigma * sys->sigma * sys->sigma * sys->sigma * sys->sigma;
    	c6 =4.0 * sys->epsilon * sigma6;
    	c12 =4.0 * sys->epsilon * sigma6 * sigma6;
    	rcsq = sys->rcut * sys->rcut;


#pragma omp parallel reduction(+:epot)
{
	int nthreads = omp_get_num_threads();
	int tid=omp_get_thread_num();
	int N = nthreads * sys->natoms;
	
 	double *fx, *fy, *fz;
	
	printf( "check ");	

   	
	double ffac, rsq;
    	double rx,ry,rz;
    	int i,j;
    /* zero energy and forces */

//	printf("%d", tid);
 	fx = sys->fx + (tid*sys->natoms); 	azzero( fx, N );
	fy = sys->fy + (tid*sys->natoms); 	azzero( fy, N );
	fz = sys->fz + (tid*sys->natoms); 	azzero( fz, N );
	

	for(i = tid; i < (sys->natoms)-1 ; i += nthreads) {
		for(j=i+1; j < (sys->natoms); ++j) {

            	/* get distance between particle i and j */
            		rx = pbc( sys->rx[i] - sys->rx[j], 0.5 * sys->box );
            		ry = pbc( sys->ry[i] - sys->ry[j], 0.5 * sys->box );
            		rz = pbc( sys->rz[i] - sys->rz[j], 0.5*sys->box);
            		rsq = rx*rx + ry*ry + rz*rz;
               	/* compute force and energy if within cutoff */
            		if (rsq < rcsq) {
				double r6, rinv; 
				rinv = 1.0/rsq; 
				r6 = rinv * rinv * rinv;
				ffac = ( 12.0 * c12 * r6 - 6.0 * c6 ) * r6 * rinv;
				epot += r6 * ( c12 * r6 - c6 );

                		fx[ i ] += rx * ffac; 		//sys->fx[i] += rx*ffac;
                		fy[ i ] += ry * ffac;		//sys->fy[i] += ry*ffac;
                		fz[ i ] += rz * ffac;		//sys->fz[i] += rz*ffac;

                		fx[ j ] -= rx * ffac;		//sys->fx[j] -= rx*ffac;
                		fy[ j ] -= ry * ffac;		//sys->fy[j] -= ry*ffac;
                		fz[ j ] -= rz * ffac; 		//sys->fz[j] -= rz*ffac;
            		}
        	}	 
    	}
	i = 1 + ( sys->natoms / nthreads );
//	printf("numth: %d ", nthreads);
//	printf("tid: %d\n", tid);

	#pragma omp barrier
	int fromidx = tid * i;
	int toidx = fromidx + i;
	if (toidx > sys->natoms) toidx = sys->natoms;
	for (i=1; i < nthreads; ++i) {
		int offs = i * sys->natoms;
		for (int j=fromidx; j < toidx; ++j) {
			sys->fx[ j ] += sys->fx[ offs+j ];
			sys->fy[ j ] += sys->fy[ offs+j ];
			sys->fz[ j ] += sys->fz[ offs+j ];
			}
		}

}
 
 sys->epot = epot;
// printf("\n I am epot  %f, and %f \n",epot , sys->epot);
}


