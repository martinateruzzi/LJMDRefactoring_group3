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

/* generic file- or pathname buffer length */
#define BLEN 200

/* structure to hold the complete information 
 * about the MD system */
struct _mdsys {
    int natoms,nfi,nsteps;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
};
typedef struct _mdsys mdsys_t;

/**UTILITIES**/

/* helper function: read a line and then return
   the first string with whitespace stripped off */
int get_a_line(FILE *fp, char *buf);

/* helper function: zero out an array */
void azzero(double *d, const int n);

/**FORCE**/

/* compute forces */
void force(mdsys_t *sys);

/**VELVERLET**/

/* first part: propagate velocities by half and positions by full step */
void update_velocities_positions(mdsys_t *sys);

/* second part: propagate velocities by another half step */
void update_velocities(mdsys_t *sys);

/* velocity verlet */
void velverlet(mdsys_t *sys);

/* compute kinetic energy */
void ekin(mdsys_t *sys);

/**OUTPUT**/

/* append data to output. */
void output(mdsys_t *sys, FILE *erg, FILE *traj);
