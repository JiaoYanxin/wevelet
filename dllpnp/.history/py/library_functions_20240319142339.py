#c++中对应函数声明import ctypes
from structures import Coeff
import ctypes
# 加载动态链接库
libpnp = ctypes.CDLL('../build/lib/libpnp_lib.so')


# 设置函数原型
libpnp.PnPDataFidelity.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), 
                                   ctypes.POINTER(Coeff), ctypes.c_int, ctypes.c_int, ctypes.c_int, 
                                   ctypes.c_int, ctypes.c_float, ctypes.POINTER(ctypes.c_float), 
                                   ctypes.POINTER(ctypes.c_float)]
# 设置函数返回类型
libpnp.PnPDataFidelity.restype = None


libpnp.TranslateAngleToCoefficients.argtypes = [ctypes.POINTER(ctypes.c_float), 
                                                 ctypes.POINTER(Coeff), 
                                                 ctypes.c_int]
libpnp.TranslateAngleToCoefficients.restype = None

libpnp.DecorateCoefficients.argtypes=[ctypes.POINTER(Coeff),
                                      ctypes.c_float,ctypes.c_float,ctypes.c_float,
                                      ctypes.c_int]
libpnp.DecorateCoefficients.restype = None

libpnp.PnPDataFidelity.argtypes = [ctypes.POINTER(ctypes.c_float),  # proj
                                   ctypes.POINTER(ctypes.c_float),  # voldata
                                   ctypes.POINTER(Coeff),           # coeffv
                                   ctypes.c_int,                    # projsZ
                                   ctypes.c_int,                    # width
                                   ctypes.c_int,                    # length
                                   ctypes.c_int,                    # height
                                   ctypes.c_float,                  # lamb
                                   ctypes.POINTER(ctypes.c_float),  # u_k
                                   ctypes.POINTER(ctypes.c_float)]  # d_k

# 设置函数返回类型
libpnp.PnPDataFidelity.restype = None

libpnp.ATOM.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p,
                        ctypes.c_int, ctypes.c_float,
                        ctypes.POINTER(ctypes.c_float),  # u_k
                        ctypes.POINTER(ctypes.c_float)]  # d_k

# 设置返回类型
# 你的函数返回void，所以这里是None
libpnp.ATOM.restype = None