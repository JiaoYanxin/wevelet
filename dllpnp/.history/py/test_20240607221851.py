# mpirun -n 1 python ./py/test.py
#参数1：x_k-1
#参数2：u_k
#参数3：d_k
#参数4：输出文件x_k
import ctypes
from library_functions import libpnp
from structures import Coeff
import numpy as np
import mrcfile
import pybind11
print(pybind11.get_include())
import sys
# from mpi4py import MPI

def ReadData(x_k_1_path,u_k_path,d_k_path):
    # tilt_path='/data/重建数据集/deepdewadge/model0/noisy/proj_noisy.mrc'
    # angle_path = '/data/重建数据集/deepdewadge/model0/proj.tlt'
    # with mrcfile.open(tilt_path,  permissive=True) as mrc:# 打开MRC文件
    #     projdata = mrc.data#z x y
    # projdata_np = np.array(projdata, dtype=np.float32)#转为np
    # projdata_flattened = projdata_np.reshape(-1)# 重塑为一维数组
    # print("projection shape: ",projdata.shape)
    # z_count = projdata.shape[0]
    # width= projdata.shape[1]
    # length= projdata.shape[2]
    # print("z_count",z_count)
    # print("width",width)
    # print("length",length)
    # print("height",height)

    # angles = []
    # with mrcfile.open(angle_path, 'r') as anglefile:# 打开MRC文件
    #     for line in anglefile:
    #     # 将每行转换为浮点数并添加到列表中
    #         angles.append(float(line.strip()))
    # angles_array = np.array(angles)


    #x_k_path = sys.argv[4]
    with mrcfile.open(u_k_path,  permissive=True) as ukmrc:
        uk =  ukmrc.data#z x y
    uk_np = np.array(uk, dtype=np.float32)#转为np
    uk_flattened = uk_np.reshape(-1)# 重塑为一维数组

    with mrcfile.open(d_k_path, permissive=True) as dkmrc:# 打开MRC文件
        dk = dkmrc.data#z x y
    dk_np = np.array(dk, dtype=np.float32)#转为np
    dk_flattened = dk_np.reshape(-1)# 重塑为一维数组

    with mrcfile.open(x_k_1_path, permissive=True) as volmrc:# 打开MRC文件
        xk = volmrc.data#z x y
    xk_np = np.array(xk, dtype=np.float32)#转为np
    xk_flattened = xk_np.reshape(-1)# 重塑为一维数组

    return uk_flattened,dk_flattened,xk_flattened

# #读入数据
# iternum=1#迭代次数

# height=200#thickness
# lamb=0.4



#projdata_flattened_ctypes = projdata_flattened.ctypes.data_as(ctypes.POINTER(ctypes.c_float))#创建一个指向原始 NumPy 数组数据的指针
#voldata_flattened_ctypes = voldata_flattened.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
# uk_flattened_ctypes = uk_flattened.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
# dk_flattened_ctypes = dk_flattened.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
#开始mpi
# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# size = comm.Get_size()
#print(MPI.Is_initialized())  # 应该输出 True
# for i in range(iternum):
#     tilt_path_bytes = tilt_path.encode('utf-8')
#     angle_path_bytes = angle_path.encode('utf-8')
#     out_path_bytes = x_k_path.encode('utf-8')
#     # 这个程序会写出vol结果 不会改变u_k d_k
#     libpnp.ATOM(voldata_flattened_ctypes,tilt_path_bytes,angle_path_bytes,out_path_bytes,38,0.4,uk_flattened_ctypes,dk_flattened_ctypes)

#MPI.Finalize()#结束mpi
# 将一维数组重塑为三维数组
#voldata = voldata_flattened.reshape((height, width, length))
#projdata_np = projdata_flattened.reshape((len(angles), width, length))
# 将该数组保存为 MRC 文件
#with mrcfile.new(out_path, overwrite=True) as mrc:
#    mrc.set_data(voldata.astype(np.float32))  # 确保数据类型为 float32
if __name__ == '__main__':
    if len(sys.argv) < 5:
        print("使用方法: python script.py x_k_path u_k_path d_k_path")
        sys.exit(1)
    x_k_1_path = sys.argv[1]
    u_k_path = sys.argv[2]
    d_k_path = sys.argv[3]
    x_k_path = sys.argv[4]
    #iternum=1#迭代次数

    height=200#thickness
    lamb=0.4
    tilt_path='/home/jiaoyx/cryorecnet/preprocessed_test.ali'
    tilt_path_bytes = tilt_path.encode('utf-8')
    angle_path = '/home/jiaoyx/cryorecnet/proj.tlt'
    angle_path_bytes = angle_path.encode('utf-8')
    
    out_path_bytes = x_k_path.encode('utf-8')

    uk_flattened,dk_flattened,voldata_flattened = ReadData(x_k_1_path,u_k_path,d_k_path)
    uk_flattened_ctypes = uk_flattened.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    dk_flattened_ctypes = dk_flattened.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    voldata_flattened_ctypes = voldata_flattened.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    libpnp.ATOM(voldata_flattened_ctypes,tilt_path_bytes,angle_path_bytes,out_path_bytes,height,lamb,uk_flattened_ctypes,dk_flattened_ctypes)




