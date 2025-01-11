
import numpy as np
#point = Point3DF(1, 2, 3)
class Point3D:#定义类
    def __init__(self, x, y, z):#构造器，在创建新的 Point3DF 类实例时自动调用
        self.x = x
        self.y = y
        self.z = z

#struct Coeff
# { // 20个双精度数，前10个是a后10个是b

# 	double a[10];
# 	double b[10];
# };
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


#void Reproject_admm_htb(const Point3DF &origin, int width, int length, int height, const Coeff &coeff, float *htb, const float *slcdata)
def htb(origin, width, length, height, coeff, htb, slcdata):
    volsize = width * length * height
    for z in range(height):
        for y in range(length):
            for x in range(width):
                coord = Point3D(x, y, z)
                wt = val_coef(origin, coord, coeff)  # Assuming val_coef returns a Weight object

                index = x + y * width + z * width * length
                if 0 <= wt.x_min < width and 0 <= wt.y_min < length:
                    n = wt.x_min + wt.y_min * width
                    htb[index] += (1 - wt.x_min_del) * (1 - wt.y_min_del) * slcdata[n]

                if 0 <= wt.x_min + 1 < width and 0 <= wt.y_min < length:
                    n = wt.x_min + 1 + wt.y_min * width
                    htb[index] += wt.x_min_del * (1 - wt.y_min_del) * slcdata[n]

                if 0 <= wt.x_min < width and 0 <= wt.y_min + 1 < length:
                    n = wt.x_min + (wt.y_min + 1) * width
                    htb[index] += (1 - wt.x_min_del) * wt.y_min_del * slcdata[n]

                if 0 <= wt.x_min + 1 < width and 0 <= wt.y_min + 1 < length:
                    n = (wt.x_min + 1) + (wt.y_min + 1) * width
                    htb[index] += wt.x_min_del * wt.y_min_del * slcdata[n]


#ATb+mu*(uk-dk)
# def ATbmuIT(atb, atb_lt, uk, dk, width, length, height, mu):
#     volsize = width * length * height
#     # 确保 atb, atb_lt, uk, dk 是 NumPy 数组
#     # atb = np.asarray(atb).reshape(height, length, width)
#     # atb_lt = np.asarray(atb_lt).reshape(height, length, width)
#     # uk = np.asarray(uk).reshape(height, length, width)
#     # dk = np.asarray(dk).reshape(height, length, width)

#     for z in range(height):
#         for y in range(length):
#             for x in range(width):
#                 atb_lt[z, y, x] = (uk[z, y, x] - dk[z, y, x]) * mu + atb[z, y, x]
def ATbmuIT(atb, atb_lt, uk, dk, mu):
    atb_lt[:] = atb + mu * (uk - dk)


def Reproject_admm_atax(origin, width, length, height, x0, coeffv, atax, proj):
    volsize = width * length * height
    ax = np.zeros((length, width))
    w = np.zeros((length, width))
    
    for idx in range(proj):
        ax.fill(0)
        w.fill(0)

        for z in range(height):
            for y in range(length):
                for x in range(width):
                    coord = Point3D(x, y, z)
                    wt = val_coef(origin, coord, coeffv[idx])  # Assuming val_coef returns a Weight object
                   # First condition
                    if 0 <= wt.x_min < width and 0 <= wt.y_min < length:
                        n = wt.x_min + wt.y_min * width
                        ax[wt.y_min, wt.x_min] += (1 - wt.x_min_del) * (1 - wt.y_min_del) *  x0[z * width * length + y * width + x]
                        w[wt.y_min, wt.x_min] += (1 - wt.x_min_del) * (1 - wt.y_min_del)

                    # Second condition
                    if 0 <= wt.x_min + 1 < width and 0 <= wt.y_min < length:
                        n = wt.x_min + 1 + wt.y_min * width
                        ax[wt.y_min, wt.x_min] += wt.x_min_del * (1 - wt.y_min_del) *  x0[z * width * length + y * width + x]
                        w[wt.y_min, wt.x_min] += wt.x_min_del * (1 - wt.y_min_del)

                    # Third condition
                    if 0 <= wt.x_min < width and 0 <= wt.y_min + 1 < length:
                        n = wt.x_min + (wt.y_min + 1) * width
                        ax[wt.y_min, wt.x_min] += (1 - wt.x_min_del) * wt.y_min_del *  x0[z * width * length + y * width + x]
                        w[wt.y_min, wt.x_min] += (1 - wt.x_min_del) * wt.y_min_del

                    # Fourth condition
                    if 0 <= wt.x_min + 1 < width and 0 <= wt.y_min + 1 < length:
                        n = (wt.x_min + 1) + (wt.y_min + 1) * width
                        ax[wt.y_min, wt.x_min] += wt.x_min_del * wt.y_min_del *  x0[z * width * length + y * width + x]
                        w[wt.y_min, wt.x_min] += wt.x_min_del * wt.y_min_del

        # Apply weighted sum to atax
        for z in range(height):
            for y in range(length):
                for x in range(width):
                    coord = Point3D(x, y, z)
                    wt = val_coef(origin, coord, coeffv[idx])
                    # First condition
                    if 0 <= wt.x_min < width and 0 <= wt.y_min < length:
                        n = wt.x_min + wt.y_min * width
                        atax[z * width * length  + y * width + x] += ((1 - wt.x_min_del) * (1 - wt.y_min_del) * ax[wt.y_min, wt.x_min] / w[wt.y_min, wt.x_min]
                                                        if w[wt.y_min, wt.x_min] != 0 else 0)

                    # Second condition
                    if 0 <= wt.x_min + 1 < width and 0 <= wt.y_min < length:
                        n = wt.x_min + 1 + wt.y_min * width
                        atax[z * width * length  + y * width + x] += (wt.x_min_del * (1 - wt.y_min_del) * ax[wt.y_min, wt.x_min + 1] / w[wt.y_min, wt.x_min + 1]
                                                        if w[wt.y_min, wt.x_min + 1] != 0 else 0)

                    # Third condition
                    if 0 <= wt.x_min < width and 0 <= wt.y_min + 1 < length:
                        n = wt.x_min + (wt.y_min + 1) * width
                        atax[z * width * length  + y * width + x] += ((1 - wt.x_min_del) * wt.y_min_del * ax[wt.y_min + 1, wt.x_min] / w[wt.y_min + 1, wt.x_min]
                                                        if w[wt.y_min + 1, wt.x_min] != 0 else 0)

                    # Fourth condition
                    if 0 <= wt.x_min + 1 < width and 0 <= wt.y_min + 1 < length:
                        n = (wt.x_min + 1) + (wt.y_min + 1) * width
                        atax[z * width * length  + y * width + x] += (wt.x_min_del * wt.y_min_del * ax[wt.y_min + 1, wt.x_min + 1] / w[wt.y_min + 1, wt.x_min + 1]
                                                        if w[wt.y_min + 1, wt.x_min + 1] != 0 else 0)