import ctypes as C

clib = C.CDLL("./libgeom.so")
clib.area.argtypes = [C.Structure]
clib.area.restype = C.c_float

class Rectangle(C.Structure):
    _fields_ = [
            ("width", C.c_float),
            ("height", C.c_float)
            ]
    def __init__(self, width, height):
        self.width = width
        self.height = height

    def area(self):
        return clib.area(self)
