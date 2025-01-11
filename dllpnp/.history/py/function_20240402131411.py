
import numpy as np
import numpy as np
from scipy.sparse.linalg import cg
import math
def d2r(degrees):
    # Convert degrees to radians
    return degrees * math.pi / 180

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
    def __init__(self):
        # Initialize a and b as arrays of zeroes
        self.a = np.zeros(10, dtype=float)
        self.b = np.zeros(10, dtype=float)

class Weight:
    def __init__(self, x_min=0, y_min=0, x_min_del=0, y_min_del=0):
        self.x_min = x_min
        self.y_min = y_min
        self.x_min_del = x_min_del
        self.y_min_del = y_min_del

def translate_angle_to_coefficients(angles):
    coeffs = [Coeff() for _ in range(len(angles))]

    for i, angle in enumerate(angles):
        beta = d2r(angle)

        coeffs[i].a[0] = 0
        coeffs[i].a[1] = math.cos(beta)
        coeffs[i].a[2] = 0
        coeffs[i].a[3] = -math.sin(beta)
        coeffs[i].b[0] = 0
        coeffs[i].b[1] = 0
        coeffs[i].b[2] = 1
        coeffs[i].b[3] = 0

    return coeffs

def decorate_coefficients(coeffs, pitch_angle, offset, zshift, anglesize):
    alpha = -math.radians(pitch_angle)
    beta = math.radians(offset)
    t = -zshift

    ca, sa = math.cos(alpha), math.sin(alpha)
    cb, sb = math.cos(beta), math.sin(beta)
    ca2, sa2 = ca * ca, sa * sa
    cb2, sb2 = cb * cb, sb * sb

    for i in range(anglesize):
        a_orig = coeffs[i].a.copy()
        b_orig = coeffs[i].b.copy()

        # Calculations for coefficients a
        coeffs[i].a[0] = a_orig[0]
        coeffs[i].a[1] = a_orig[2] * sa * sb + a_orig[3] * ca * sb + a_orig[1] * cb
        coeffs[i].a[2] = a_orig[2] * ca - a_orig[3] * sa
        coeffs[i].a[3] = a_orig[3] * ca * cb + a_orig[2] * cb * sa - a_orig[1] * sb
        coeffs[i].a[4] = a_orig[4] * ca * cb - a_orig[5] * cb * sa - a_orig[6] * sa2 * sb - 2 * a_orig[9] * ca * sa * sb + a_orig[6] * ca2 * sb + 2 * a_orig[8] * ca * sa * sb
        coeffs[i].a[5] = a_orig[4] * cb2 * sa - a_orig[5] * ca * sb2 + a_orig[5] * ca * cb2 - 2 * a_orig[7] * cb * sb - a_orig[4] * sa * sb2 + 2 * a_orig[9] * ca2 * cb * sb + 2 * a_orig[8] * cb * sa2 * sb + 2 * a_orig[6] * ca * cb * sa * sb
        coeffs[i].a[6] = 2 * a_orig[8] * ca * cb * sa - a_orig[4] * ca * sb + a_orig[5] * sa * sb + a_orig[6] * ca2 * cb - a_orig[6] * cb * sa2 - 2 * a_orig[9] * ca * cb * sa
        coeffs[i].a[7] = a_orig[7] * cb2 + a_orig[5] * ca * cb * sb + a_orig[9] * ca2 * sb2 + a_orig[8] * sa2 * sb2 + a_orig[4] * cb * sa * sb + a_orig[6] * ca * sa * sb2
        coeffs[i].a[8] = a_orig[8] * ca2 + a_orig[9] * sa2 - a_orig[6] * ca * sa
        coeffs[i].a[9] = a_orig[7] * sb2 + a_orig[9] * ca2 * cb2 + a_orig[8] * cb2 * sa2 - a_orig[4] * cb * sa * sb + a_orig[6] * ca * cb2 * sa - a_orig[5] * ca * cb * sb
 # Calculations for coefficients b
        coeffs[i].b[0] = b_orig[0]
        coeffs[i].b[1] = b_orig[2] * sa * sb + b_orig[3] * ca * sb + b_orig[1] * cb
        coeffs[i].b[2] = b_orig[2] * ca - b_orig[3] * sa
        coeffs[i].b[3] = b_orig[3] * ca * cb + b_orig[2] * cb * sa - b_orig[1] * sb
        coeffs[i].b[4] = b_orig[4] * ca * cb - b_orig[5] * cb * sa - b_orig[6] * sa2 * sb - 2 * b_orig[9] * ca * sa * sb + b_orig[6] * ca2 * sb + 2 * b_orig[8] * ca * sa * sb
        coeffs[i].b[5] = b_orig[4] * cb2 * sa - b_orig[5] * ca * sb2 + b_orig[5] * ca * cb2 - 2 * b_orig[7] * cb * sb - b_orig[4] * sa * sb2 + 2 * b_orig[9] * ca2 * cb * sb + 2 * b_orig[8] * cb * sa2 * sb + 2 * b_orig[6] * ca * cb * sa * sb
        coeffs[i].b[6] = 2 * b_orig[8] * ca * cb * sa - b_orig[4] * ca * sb + b_orig[5] * sa * sb + b_orig[6] * ca2 * cb - b_orig[6] * cb * sa2 - 2 * b_orig[9] * ca * cb * sa
        coeffs[i].b[7] = b_orig[7] * cb2 + b_orig[5] * ca * cb * sb + b_orig[9] * ca2 * sb2 + b_orig[8] * sa2 * sb2 + b_orig[4] * cb * sa * sb + b_orig[6] * ca * sa * sb2
        coeffs[i].b[8] = b_orig[8] * ca2 + b_orig[9] * sa2 - b_orig[6] * ca * sa
        coeffs[i].b[9] = b_orig[7] * sb2 + b_orig[9] * ca2 * cb2 + b_orig[8] * cb2 * sa2 - b_orig[4] * cb * sa * sb + b_orig[6] * ca * cb2 * sa - b_orig[5] * ca * cb * sb

        # ... (z_shift adjustments if applicable)

        
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
                        
def ATAmuITI(vol, ata, mu):
    # Perform the operation in a vectorized manner
    ltl = ata + mu * vol
    return ltl
