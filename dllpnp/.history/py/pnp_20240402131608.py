import function as f
import numpy as np

#vol是z*x*y的
def pnp(vol):
    #计算光线
    f.translate_angle_to_coefficients(angles.data(), params.data(), angles.size());

    #先计算htb
    length = vol.shape[0]   # '深度'或第一个维度
    height = vol.shape[1]  # '高度'或第二个维度
    width = vol.shape[2]   # '宽度'或第三个维度
    ori=f.Point3D(width*0.5,length*0.5,height*0.5)
    htb(ori, width, length, height, coeff, htb, slcdata)
