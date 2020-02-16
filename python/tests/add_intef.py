import numpy as np
import ctypes

math = ctypes.CDLL("./libmymath.so")

math.add_float.restype = ctypes.c_float
math.add_float.argtypes = [ctypes.c_float, ctypes.c_float]
print(math.add_int(4, 5))
print(math.add_float(10.0, 20.0))

a = ctypes.c_float(3.0)
b = ctypes.c_float(7.0)
res = ctypes.c_float()
math.add_float_ref(ctypes.byref(a), ctypes.byref(b), ctypes.byref(res))
print(res.value) 

i = ctypes.pointer(a)
j = ctypes.pointer(b)
k = ctypes.pointer(res)
math.add_float_ref(i, j, k)
print(k.contents)

a = (ctypes.c_int * 3)(2, 4, 7)
b= (ctypes.c_int * 3)(3, 1, -2)
c = (ctypes.c_int * 3)(0, 0, 0)
n = ctypes.c_int(3)
math.add_int_array(a, b, c, n)
print(c[0], c[1], c[2])

a = np.array([10, 20, 30], dtype=ctypes.c_int)
b = np.array([30, 20, 10], dtype=ctypes.c_int)
res = np.zeros(3, dtype=ctypes.c_int)
n = ctypes.c_int(3)
intp = ctypes.POINTER(ctypes.c_int)
i = a.ctypes.data_as(intp)
j = b.ctypes.data_as(intp)
k = res.ctypes.data_as(intp)
math.add_int_array(i, j, k, n)
print(res)
