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
    double ffac, sigma6, c6, c12, rcsq, rsq;
    double rx,ry,rz;
    int i,j;

    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    sigma6 = sys->sigma*sys->sigma*sys->sigma*sys->sigma*sys->sigma*sys->sigma;
    c6 =4.0*sys->epsilon*sigma6;
    c12 =4.0*sys->epsilon*sigma6*sigma6;
    rcsq = sys->rcut * sys->rcut;

    for(i=0; i < (sys->natoms)-1 ; ++i) {
        for(j=i+1; j < (sys->natoms); ++j) {

            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
            rsq = rx*rx + ry*ry + rz*rz;
      
            /* compute force and energy if within cutoff */
            if (rsq < rcsq) {

		double r6,rinv; rinv=1.0/rsq; r6=rinv*rinv*rinv;
		ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
		sys->epot += r6*(c12*r6 - c6);

                sys->fx[i] += rx*ffac;
                sys->fy[i] += ry*ffac;
                sys->fz[i] += rz*ffac;

                sys->fx[j] -= rx*ffac;
                sys->fy[j] -= ry*ffac;
                sys->fz[j] -= rz*ffac;
            }
        }
    }
}
