 from ctypes import *
import numpy as np
import sys as _sys

class Mdsys(Structure):
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
        ("ekinetic", c_double),
        ("epotential", c_double),
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
    
    def __init__(self, initfile):
        self.initfile = initfile
        self.loadinit()
        self.rvfinit()
        self.loadengine()
    
    def loadinit(self):
        try:
            self.args = []
            with open(self.initfile, 'r') as file:
                for line in file:
                    line = line.split("#", 1)[0]
                    line = line.rstrip()
                    self.args.append(line)
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
        self.nfi = 0

    def reloadinit(self, initfile):
        self.initfile = initfile
        self.loadinit()

    def rvfinit(self):
        try:
            _initarr = np.loadtxt(self.restfile, dtype=c_double)
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
        self.rx = _rx.ctypes.data_as(_dblp)
        self.ry = _ry.ctypes.data_as(_dblp)
        self.rz = _rz.ctypes.data_as(_dblp)
        self.vx = _vx.ctypes.data_as(_dblp)
        self.vy = _vy.ctypes.data_as(_dblp)
        self.vz = _vz.ctypes.data_as(_dblp)
        self.fx = _fx.ctypes.data_as(_dblp)
        self.fy = _fy.ctypes.data_as(_dblp)
        self.fz = _fz.ctypes.data_as(_dblp)

    def output(self):
        _ergstr = "{:8} {:20.8} {:20.8} {:20.8} {:20.8}\n".format(
            self.nfi, self.temp, self.ekinetic, self.epotential, 
            self.ekinetic+self.epotential)
        print(_ergstr)
        try:
            with open(self.ergfile, "a") as file:
                file.write(_ergstr)
        except Exception as err:
            print("Error writing to file: {}".format(str(err)))

        try:
            with open(self.trajfile, "a") as file:
                file.write("nfil = {:8} \t etot = {:.8}\n".format(
                    self.nfi, self.ekinetic+self.epotential))
                for i in range(self.natoms):
                    _trajstr = "Ar  {:20.8} {:20.8} {:20.8}\n".format(
                        self.rx[i], self.ry[i], self.rz[i])
                    file.write(_trajstr)
        except Exception as err:
            print("Error writing to file: {}".format(str(err)))

    def loadengine(self):
        self._dll = CDLL("../libljmd.so")
        self._dll.force.argtypes = [POINTER(Mdsys)]
        self._dll.force.restype = None
        self._dll.ekin.argtypes = [POINTER(Mdsys)]
        self._dll.ekin.restype = None
        self._dll.velverlet.argtypes = [POINTER(Mdsys)]
        self._dll.velverlet.restype = None

    def force(self):
        self._dll.force(byref(self))

    def ekin(self):
        self._dll.ekin(byref(self))

    def velverlet(self):
        self._dll.velverlet(byref(self))