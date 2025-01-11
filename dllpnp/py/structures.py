#c++中对应结构体
import ctypes


class Point2D:#point = Point2D()的默认值为0
    def __init__(self, x=0, y=0):#类的构造器 
        self.x = x
        self.y = y

    def __repr__(self):#用于打印 print(point)  # 输出: Point2D(5, 10)
        return f"Point2D(x={self.x}, y={self.y})"
    
class Point3D:
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return f"Point3D(x={self.x}, y={self.y}, z={self.z})"

class Slice: #slice_instance = Slice(10, 20)
    def __init__(self, width, height, data=None, x=0, y=0):
        self.coord = Point2D(x, y)
        self.width = width
        self.height = height

        if data is None:#初始化
            self.data = [0.0] * (width * height)#浮点数列表
            self.external = False
        else:
            self.data = data
            self.external = True

    def set_coord(self, x, y):
        self.coord.x = x
        self.coord.y = y

    # 在Python中，通常不需要显式的析构函数，
    # 因为Python有垃圾回收机制来处理内存管理。
# 在 coeff_definitions.py 中
import ctypes

class Coeff(ctypes.Structure):
    _fields_ = [("a", ctypes.c_double * 10), 
                ("b", ctypes.c_double * 10)]