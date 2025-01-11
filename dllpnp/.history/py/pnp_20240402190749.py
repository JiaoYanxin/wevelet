import function as f
import numpy as np

#vol是z*x*y的
#angle np、列表都行
def pnp(vol, angles,projs):
    #计算光线
    params=f.translate_angle_to_coefficients(angles);
    f.decorate_coefficients(params, 0, 0, 0,  len(angles));
    #先计算htb
    height = vol.shape[0]   # '深度'或第一个维度
    length = vol.shape[1]  # '高度'或第二个维度
    width = vol.shape[2]   # '宽度'或第三个维度
    ori=f.Point3D(width*0.5,length*0.5,height*0.5)
    htb=np.zeros(vol)
    for idx in range(height):
        proj = projs[idx]  # 提取当前切片
        f.htb(ori,width,length,height,params[idx],htb,proj)
    print("**  PnP with Data Fidelity item iteration,iter **\n");
    #共轭梯度法的初值就是vol.data
    
    htb(ori, width, length, height, coeff, htb, slcdata)
