from ctypes import *
import numpy as np
import sys as _sys

""" ljmd.py: Wrapper for ljmd.c code
This is a wrapper around Leonard-Jones Molecular Dynamics codes 
written in c. for a description on how to use the wrapper 
refer to the project website at:
https://github.com/saliei/LJMDRefactoring_group3
"""
__author__      =   "Saeid Aliei"
__copyright__   =   "Copyright 2020, ICTP-MHPC"
__licence__     =   "GPL"
__version__     =   "1.0.0"
__maintainer__  =   "Saeid Aliei"
__email__       =   "saeidaliei2019@gmail.com"
__date__        =   "Dev"

class Mdsys(Structure):
    """
    Wrapper around _mdsys structure.
    
    Attributes
    ----------
    natoms: int
        number of atoms present in the system
    nfi: int
        step of the simulation
    nsteps: int
        number of steps to evolve the system
    dt: float
        time step of simulation
    mass: float
        mass of each molecule in AMU units
    epsilon: float
        constant in kcal/mol units
    rcut: float
        potential reach, after this force is zero, in units of Angstrom
    ekin: float
        kinetic energy of the system
    epot: float
        potential energy of the system
    temp: float
        temperature of the system
    rx: pointer(double)
        pointer to the array of x positions of the system
    ry: pointer(double)
        pointer to the array of y positions of the system
    rz: pointer(double)
        pointer to the array of z positions of the system
    vx: pointer(double)
        pointer to the array of x velocities of the system
    vy: pointer(double)
        pointer to the array of y velocities of the system
    vz: pointer(double)
        pointer to the array of z velocities of the system
    fx: pointer(double)
        pointer to the array of x velocities of the system
    fy: pointer(double)
        pointer to the array of y velocities of the system
    fz: pointer(double)
        pointer to the array of z velocities of the system
    """
    _fields_ = [
        ("natoms", c_int),
        ("nfi", c_int),
        ("nsteps", c_int),
        ("dt", c_double),
        ("mass", c_double),
        ("epsilon", c_double),
        ("sigma", c_double),
        ("box", c_double),
        ("rcut", c_double),
        ("ekin", c_double),
        ("epot", c_double),
        ("temp", c_double),
        ("rx", POINTER(c_double)),
        ("ry", POINTER(c_double)),
        ("rz", POINTER(c_double)),
        ("vx", POINTER(c_double)),
        ("vy", POINTER(c_double)),
        ("vz", POINTER(c_double)),
        ("fx", POINTER(c_double)),
        ("fy", POINTER(c_double)),
        ("fz", POINTER(c_double)),
        ]

class LJMD:
    """
    LJMD object to define a system of particles.
    """
    def __init__(self, initfile):
        self.initfile = initfile
        self.loadinit()
        self.sysinit()

    def loadinit(self):
        try:
            with open("examples/"+self.initfile, 'r') as file:
                self.args = [line.split("#","1")[0].rstrip() for line in file]
        except Exception as err:
            print("Error reading init file: {}".format(str(err)))
            _sys.exit(1)
        self.natoms = int(self.args[0])
        self.mass = float(self.args[1])
        self.epsilon = float(self.args[2])
        self.sigma = float(self.args[3])
        self.rcut = float(self.args[4])
        self.box = float(self.args[5])
        self.restfile = str(self.args[6])
        self.trajfile = str(self.args[7])
        self.ergfile = str(self.args[8])
        self.nsteps = int(self.args[9])
        self.dt = float(self.args[10])
        self.nprint = int(self.args[11])
        print("Loading init file.\n")

    def sysinit(self):
        self.sys = Mdsys(
            natoms = self.natoms,
            nsteps = self.nsteps,
            dt = self.dt,
            mass = self.mass,
            epsilon = self.epsilon,
            sigma = self.sigma,
            box = self.box,
            rcut = self.rcut,            
            )
        try:
            _initarr = np.loadtxt("examples/"+self.restfile, dtype=c_double)
        except Exception as err:
            print("Error reading restart file: {}".format(str(err)))
            _sys.exit(1)
        _dblp = POINTER(c_double)
        _rx = _initarr[:,0][:self.natoms]
        _ry = _initarr[:,1][:self.natoms]
        _rz = _initarr[:,2][:self.natoms]
        _vx = _initarr[:,0][self.natoms:]
        _vy = _initarr[:,1][self.natoms:]
        _vz = _initarr[:,2][self.natoms:]
        _fx = np.zeros(self.natoms, dtype=c_double)
        _fy = np.zeros(self.natoms, dtype=c_double)
        _fz = np.zeros(self.natoms, dtype=c_double)
        self.sys.rx = _rx.ctypes.data_as(_dblp)
        self.sys.ry = _ry.ctypes.data_as(_dblp)
        self.sys.rz = _rz.ctypes.data_as(_dblp)
        self.sys.vx = _vx.ctypes.data_as(_dblp)
        self.sys.vy = _vy.ctypes.data_as(_dblp)
        self.sys.vz = _vz.ctypes.data_as(_dblp)
        self.sys.fx = _fx.ctypes.data_as(_dblp)
        self.sys.fy = _fy.ctypes.data_as(_dblp)
        self.sys.fz = _fz.ctypes.data_as(_dblp)
        self.sys.nfi = 0
        self.force()
        self.ekin()
        print("System initialized.\n")

    def output(self):
        _ergstr = "{:8} {:20.8} {:20.8} {:20.8} {:20.8}\n".format(
            self.sys.nfi, self.sys.temp, self.sys.ekin, self.sys.epot, 
            self.sys.ekin+self.sys.epot)
        print(_ergstr)
        try:
            with open("refrences/"+self.ergfile, "a") as file:
                file.write(_ergstr)
        except Exception as err:
            print("Error writing to file: {}".format(str(err)))

        try:
            with open("refrences/"+self.trajfile, "a") as file:
                file.write("nfil = {:8} \t etot = {:.8}\n".format(
                    self.sys.nfi, self.sys.ekin+self.sys.epot))
                for i in range(self.natoms):
                    _trajstr = "Ar  {:20.8} {:20.8} {:20.8}\n".format(
                        self.sys.rx[i], self.sys.ry[i], self.sys.rz[i])
                    file.write(_trajstr)
        except Exception as err:
            print("Error writing to file: {}".format(str(err)))

    def loadengine(self):
        self._dll = CDLL("../libljmd.so")
        self._dll.force.argtypes = [POINTER(self.sys)]
        self._dll.force.restype = None
        self._dll.ekin.argtypes = [POINTER(self.sys)]
        self._dll.ekin.restype = None
        self._dll.velverlet.argtypes = [POINTER(self.sys)]
        self._dll.velverlet.restype = None

    def force(self):
        self._dll.force(byref(self.sys))

    def ekin(self):
        self._dll.force(byref(self.sys))

    def velverlet(self):
        self._dll.velverlet(byref(self.sys))