#参数1：x_k
#参数1：u_k(网络输出结果)
#参数3：d_k1
#参数4：d_k(输出结果)
import mrcfile
import numpy as np
import sys
# 检查是否提供了足够的参数
if len(sys.argv) < 5:
    print("使用方法: python script.py x_k_path u_k_path d_k_path")
    sys.exit(1)

# 获取参数
x_k_path = sys.argv[1]
u_k_path = sys.argv[2]
d_k_1_path = sys.argv[3]
d_k_path = sys.argv[4]

# 您的脚本逻辑...

#input
#x_k_path="/data/重建数据集/deepdewadge/model0/noisy/reconfull_noisy_r.mrc"
#d_k_1_path=""
#u_k_path="/data/重建数据集/deepdewadge/model0/result/refined_tomogram.mrc"
#output
#d_k_path="/data/重建数据集/deepdewadge/model0/result/0_dk.mrc"

#第一次 dk-1=0
with mrcfile.open(x_k_path, permissive=True) as xkmrc:
    x_k_data = xkmrc.data
with mrcfile.open(d_k_1_path) as dkmrc:
    d_k_1_data = dkmrc.data
with mrcfile.open(u_k_path, permissive=True) as ukmrc:
    u_k_data = ukmrc.data

dk=(x_k_data- u_k_data)*1.5+d_k_1_data

#写出
with mrcfile.open(x_k_path, mode='r', permissive=True) as mrc:
    x_k_data = mrc.data.copy()  # 使用 copy 确保从文件中读取数据

# 创建一个新的MRC文件
with mrcfile.new(d_k_path, overwrite=True) as new_mrc:
    # 使用 x_k 的 header 作为基础
    new_mrc.set_data(dk.astype(np.float32))

    # Copy header information from the original file
    # Instead of directly assigning the header, update relevant fields individually
    original_header = mrc.header
    #new_mrc.header.cella = original_header.cella
   # new_mrc.header.cellb = original_header.cellb
   # new_mrc.header.cellc = original_header.cellc
   # new_mrc.header.mapc = original_header.mapc
   # new_mrc.header.mapr = original_header.mapr
   # new_mrc.header.maps = original_header.maps
    # ... copy other necessary header fields in the same way

    # Update statistical information
    new_mrc.header.dmin = np.min(dk)
    new_mrc.header.dmax = np.max(dk)
    new_mrc.header.dmean = np.mean(dk)
