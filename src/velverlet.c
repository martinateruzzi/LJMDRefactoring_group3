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

/* a few physical constants */
const double kboltz=0.0019872067;     /* boltzman constant in kcal/mol/K */
const double mvsq2e=2390.05736153349; /* m*v^2 in kcal/mol */

/* first part: propagate velocities by half and positions by full step */
void update_velocities_positions(mdsys_t *sys)
{
    int i;

    /* first part: propagate velocities by half and positions by full step */
    for (i=0; i<sys->natoms; ++i) {
	//printf("loop %d",i);
	//printf("natoms %d",sys->natoms);
	//printf("for loop velverlet ok %d total",i);
        sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        //printf("for loop velverlet ok %d fx",i);
	sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        //printf("for loop velverlet ok %d fy",i);
	sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass;
        //printf("for loop velverlet ok %d fz",i);
	sys->rx[i] += sys->dt*sys->vx[i];
        //printf("for loop velverlet ok %d vx",i);
	sys->ry[i] += sys->dt*sys->vy[i];
        //printf("for loop velverlet ok %d vy",i);
	sys->rz[i] += sys->dt*sys->vz[i];
	//printf("for loop velverlet ok %d vz",i);
    }

}

/* second part: propagate velocities by another half step */
void update_velocities(mdsys_t *sys)
{
    int i;

    /* second part: propagate velocities by another half step */
    for (i=0; i<sys->natoms; ++i) {
	   // printf("second loop velverlet %d",i);
        sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass;
    }

}

/* velocity verlet */
void velverlet(mdsys_t *sys)
{

	/* first part: propagate velocities by half and positions by full step */
	update_velocities_positions(sys);

	/* compute forces and potential energy */
	force(sys);

	/* second part: propagate velocities by another half step */
	update_velocities(sys);
}

/* compute kinetic energy */
void ekin(mdsys_t *sys)
{   
    int i;

    sys->ekin=0.0;
    for (i=0; i<sys->natoms; ++i) {
        sys->ekin += 0.5*mvsq2e*sys->mass*(sys->vx[i]*sys->vx[i] + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
    }
    sys->temp = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}
