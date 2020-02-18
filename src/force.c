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
  int i,j,ii;

  /* zero energy and forces */
  double epot=0.0;
  azzero(sys->cx,sys->natoms);
  azzero(sys->cy,sys->natoms);
  azzero(sys->cz,sys->natoms);

  sigma6 = sys->sigma*sys->sigma*sys->sigma*sys->sigma*sys->sigma*sys->sigma;
  c6 =4.0*sys->epsilon*sigma6;
  c12 =4.0*sys->epsilon*sigma6*sigma6;
  rcsq = sys->rcut * sys->rcut;
    
  MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);
  MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);
  MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, sys->mpicomm);

  //for(i=0; i < (sys->natoms)-1 ; i+=sys->nsize) {
   // ii=i+sys->mpirank;
   // if (ii >= (sys->natoms - 1)) break;
  for(i=sys->mpirank; i < sys->natoms; i+=sys->nprocs) {  
    for(j=i+1; j < (sys->natoms); ++j) {

      /* get distance between particle i and j */
      rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
      ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
      rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
      rsq = rx*rx + ry*ry + rz*rz;
      
      /* compute force and energy if within cutoff */
      if (rsq < rcsq) {

	double r6,rinv;
	rinv=1.0/rsq;
	r6=rinv*rinv*rinv;
	ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
	epot += r6*(c12*r6 - c6); /*change epot is fo every prcessor, then we will do the reduce*/

	sys->cx[i] += rx*ffac;
	sys->cy[i] += ry*ffac;
	sys->cz[i] += rz*ffac;

	sys->cx[j] -= rx*ffac;
	sys->cy[j] -= ry*ffac;
	sys->cz[j] -= rz*ffac;
      }
    }
  }
    MPI_Reduce(sys->cx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
    MPI_Reduce(sys->cy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
    MPI_Reduce(sys->cz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);
    MPI_Reduce(&epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, sys->mpicomm);

}
