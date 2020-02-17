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
            self.nfi, self.temp, self.ekin, self.epot, self.ekin+self.epot)
        print(_ergstr)
        try:
            with open(self.ergfile, "a") as file:
                file.write(_ergstr)
        except Exception as err:
            print("Error reading file: {}".format(str(err)))

        try:
            with open(self.trajfile, "a") as file:
                file.write("nfil = {:8} \t etot = {:.8}".format(
                    self.nfi, self.ekin+self.epot))
                

    def loadengine(self):
        self._dll = CDLL("./libljmd.so")
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


sys = Mdsys(initfile="../examples/argon_108.inp")
sys.loadengine()

# sys = Mdsys()
print(sys.natoms)
print(sys.force())
# print(sys.ekin())
"""
libkinetic = CDLL("./libkinetic.so")
ekin = libkinetic.ekin
ekin.argtypes = [POINTER(Mdsys)]
ekin.restype = None

libforce = CDLL("./libforce.so")
force = libforce.force
force.argtypes = [POINTER(Mdsys)]
force.restype = None

natoms = 3
mass = 39.948
epsilon = 0.2379
sigma = 3.405
rcut = 8.5
box = 17.1580
rx = np.array([6.67103294321331, -10.6146871435653, 12.6336939877734], dtype=c_double)
ry = np.array([1.06574058650169, -3.33432278188177, -2.59038677851747], dtype=c_double)
rz = np.array([-1.78412295775301, -16.5259458407765, 4.61680014503288], dtype=c_double)
vx = np.array([-1.5643224621482283e-03, 4.8497508563925346e-04, -4.3352481732883966e-04], dtype=c_double)
vy = np.array([4.1676710257651452e-04, 2.2858522230176587e-05, -6.1985040462745732e-04], dtype=c_double)
vz = np.array([-7.5611349562333923e-04, 4.0710138209103827e-04, -4.6520198934056357e-04], dtype=c_double)
nsteps = 100
dt = 5
nprint = 5
fx = np.zeros(natoms, dtype=c_double)
fy = np.zeros(natoms, dtype=c_double)
fz = np.zeros(natoms, dtype=c_double)
nfi = 0

sys = Mdsys()
sys.natoms = natoms
sys.mass = mass
sys.epsilon = epsilon
sys.sigma = sigma
sys.rcut = rcut
sys.box = box
dblp = POINTER(c_double)
rxp = rx.ctypes.data_as(dblp)
ryp = ry.ctypes.data_as(dblp)
rzp = rz.ctypes.data_as(dblp)
vxp = vx.ctypes.data_as(dblp)
vyp = vy.ctypes.data_as(dblp)
vzp = vz.ctypes.data_as(dblp)
fxp = fx.ctypes.data_as(dblp)
fyp = fy.ctypes.data_as(dblp)
fzp = fz.ctypes.data_as(dblp)
sys.rx = rxp
sys.ry = ryp
sys.rz = rzp
sys.vx = vxp
sys.vy = vyp
sys.vz = vzp
sys.fx = fxp
sys.fy = fyp
sys.fz = fzp
sys.nfi = nfi
sys.nsteps = nsteps
sys.dt = dt

ekin(byref(sys))
print(sys.ekin)
print(sys.temp)
for i in range(natoms):
    print(sys.fx[i], sys.fy[i], sys.fz[i])
force(byref(sys))
for i in range(natoms):
    print(sys.fx[i], sys.fy[i], sys.fz[i])
nprint = int()
"""