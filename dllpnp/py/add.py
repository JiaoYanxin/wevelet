#参数1：x_k
#参数1：u_k(网络输出结果)
#参数3：d_k1
#参数4：d_k(输出结果)
import mrcfile
import numpy as np
import sys
# 检查是否提供了足够的参数
if len(sys.argv) < 3:
    print("使用方法: python script.py x_k_path u_k_path d_k_path")
    sys.exit(1)

# 获取参数
x_k_path = sys.argv[1]
u_k_path = sys.argv[2]

out_path = sys.argv[3]

#第一次 dk-1=0
with mrcfile.open(x_k_path, permissive=True) as xkmrc:
    x_k_data = xkmrc.data
with mrcfile.open(u_k_path, permissive=True) as ukmrc:
    u_k_data = ukmrc.data

dk=x_k_data+ u_k_data

#写出
with mrcfile.open(x_k_path, mode='r', permissive=True) as mrc:
    x_k_data = mrc.data.copy()  # 使用 copy 确保从文件中读取数据

# 创建一个新的MRC文件
with mrcfile.new(out_path, overwrite=True) as new_mrc:
    # 使用 x_k 的 header 作为基础
    new_mrc.set_data(dk.astype(np.float32))

    # Copy header information from the original file
    # Instead of directly assigning the header, update relevant fields individually
    original_header = mrc.header
    new_mrc.header.dmin = np.min(dk)
    new_mrc.header.dmax = np.max(dk)
    new_mrc.header.dmean = np.mean(dk)
