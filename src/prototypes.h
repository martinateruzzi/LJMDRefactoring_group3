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

#include "datastructures.h"

/* generic file- or pathname buffer length */
#define BLEN 200

/**UTILITIES**/

/* helper function: read a line and then return
   the first string with whitespace stripped off */
int get_a_line(FILE *fp, char *buf);

/* helper function: zero out an array */
void azzero(double *d, const int n);

/* helper function: apply minimum image convention */
double pbc(double x, const double boxby2);

/**FORCE**/

/* compute forces */
void force(mdsys_t *sys);

/**VELVERLET**/

/* velocity verlet */
void velverlet(mdsys_t *sys);

/* compute kinetic energy */
void ekin(mdsys_t *sys);

/**OUTPUT**/

/* append data to output. */
void output(mdsys_t *sys, FILE *erg, FILE *traj);
