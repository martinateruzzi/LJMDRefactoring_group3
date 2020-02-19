/* 
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#ifndef DATASTRUCTURE_H
#define DATASTRUCTURE_H

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/* structure to hold the complete information 
 * about the MD system */
struct _mdsys {
  int natoms,nfi,nsteps;
  double dt, mass, epsilon, sigma, box, rcut;
  double ekin, epot, temp;
  double *rx, *ry, *rz;
  double *vx, *vy, *vz;
  double *fx, *fy, *fz;
  double *cx, *cy, *cz;
  int nsize, mpirank,nprocs, nthreads;
  MPI_Comm mpicomm;
};
typedef struct _mdsys mdsys_t;

#endif
