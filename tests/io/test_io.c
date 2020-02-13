#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "prototypes.h"

int main (int argc, char **argv )
{
 
  int nprint, i;
  char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
  FILE *fp,*traj,*erg;
  mdsys_t sys;

  /* TEST SYS*/
  double check_array1[]= {108.0,39.948,0.2379,3.405,8.5,17.1580,10000.0,5.0,100.0};
  int flag = 0;  
  if(get_a_line(stdin,line)) return 1;
  sys.natoms=atof(line);
  if(check_array1[0]!=sys.natoms)
    {
      printf("Reading inputs error natoms");
      flag = 1;
    }
  
  if(get_a_line(stdin,line)) return 1;
  sys.mass=atof(line);
  if(check_array1[1]!=sys.mass)
    {
      printf("Reading inputs error mass");
      flag = 1;
    }
 
  if(get_a_line(stdin,line)) return 1;
  sys.epsilon=atof(line);
  if(check_array1[2]!=sys.epsilon)
    {
      printf("Reading inputs error epsilon");
      flag = 1;
    }
  if(get_a_line(stdin,line)) return 1;
  sys.sigma=atof(line);
  if(check_array1[3]!=sys.sigma)
    {
      printf("Reading inputs error sigma");
      flag = 1;
    }
  if(get_a_line(stdin,line)) return 1;
  sys.rcut=atof(line);
  if(check_array1[4]!=sys.rcut)
    {
      printf("Reading inputs error rcut");
      flag = 1;
    }
  get_a_line(stdin,line);
  sys.box=atof(line);
  if(check_array1[5]!=sys.box)
    {
      printf("Reading inputs error box");
      flag = 1;
    }

  if(get_a_line(stdin,restfile)) return 1;
  if(get_a_line(stdin,trajfile)) return 1;
  if(get_a_line(stdin,ergfile)) return 1;	
  if(get_a_line(stdin,line)) return 1;
  sys.nsteps=atoi(line);
  if(check_array1[6]!=sys.nsteps)
    {
      printf("Reading inputs error nsteps");
      flag = 1;
    }
  if(get_a_line(stdin,line)) return 1;
  sys.dt=atof(line);
  if(check_array1[7]!=sys.dt)
    {
      printf("Reading inputs error dt");
      flag = 1;
    }
  if(get_a_line(stdin,line)) return 1;
  nprint=atoi(line);
  if(check_array1[8]!=nprint)
    {
      printf("Reading inputs error nprint");
      flag = 1;
    }


  sys.rx=(double *)malloc(sys.natoms*sizeof(double));
  sys.ry=(double *)malloc(sys.natoms*sizeof(double));
  sys.rz=(double *)malloc(sys.natoms*sizeof(double));
  sys.vx=(double *)malloc(sys.natoms*sizeof(double));
  sys.vy=(double *)malloc(sys.natoms*sizeof(double));
  sys.vz=(double *)malloc(sys.natoms*sizeof(double));
  sys.fx=(double *)malloc(sys.natoms*sizeof(double));
  sys.fy=(double *)malloc(sys.natoms*sizeof(double));
  sys.fz=(double *)malloc(sys.natoms*sizeof(double));


  /* TEST restart */
  fp=fopen(restfile,"r");

  double checkposition[] = {-6.65121935702062,1.49356278110318,20.6354475761263};
  double checkvelocity[] = {-7.6107091738481837e-05,6.5742688995404971e-04,-2.9040212591161766e-04};
  
 
  if(fp) {
    for (i=0; i<sys.natoms; ++i) {
      fscanf(fp,"%lf%lf%lf",sys.rx+i, sys.ry+i, sys.rz+i);
    }

    if(sys.rz[107]!=checkposition[2])
      {
	printf("wrong position check");
      flag = 1;
      }
	      
    for (i=0; i<sys.natoms; ++i) {
      fscanf(fp,"%lf%lf%lf",sys.vx+i, sys.vy+i, sys.vz+i);
    }

    if(sys.vz[107]!=checkvelocity[2])
      {
	printf("wrong velocity check");
      flag = 1;
      }
	    
    fclose(fp);
  } else {
    perror("cannot read restart file");
    return 3;
  }
  if (flag == 0)
  {
	  printf("No input Error!\n");
  }
 
}
