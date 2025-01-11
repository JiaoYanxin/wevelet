#point = Point3DF(1, 2, 3)
class Point3DF:#定义类
    def __init__(self, x, y, z):#构造器，在创建新的 Point3DF 类实例时自动调用
        self.x = x
        self.y = y
        self.z = z

class Point3D:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

class Coeff:
    def __init__(self, a, b):
        self.a = a  # Assuming a is a list [a0, a1, a2, a3]
        self.b = b  # Assuming b is a list [b0, b1, b2, b3]

class Weight:
    def __init__(self, x_min=0, y_min=0, x_min_del=0, y_min_del=0):
        self.x_min = x_min
        self.y_min = y_min
        self.x_min_del = x_min_del
        self.y_min_del = y_min_del

def val_coef(origin, coord, coeff):
    X = coord.x - origin.x
    Y = coord.y - origin.y
    Z = coord.z - origin.z

    n = [0, 0]
    n[0] = coeff.a[0] + coeff.a[1] * X + coeff.a[2] * Y + coeff.a[3] * Z
    n[1] = coeff.b[0] + coeff.b[1] * X + coeff.b[2] * Y + coeff.b[3] * Z

    x = n[0] + origin.x
    y = n[1] + origin.y

    weight = Weight()
    weight.x_min = int(x)  # floor
    weight.y_min = int(y)  # floor

    weight.x_min_del = x - weight.x_min
    weight.y_min_del = y - weight.y_min

    return weight
