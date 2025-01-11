import mrcfile
import numpy as np

#input
# 获取参数
x_k_path = sys.argv[1]
d_k_1_path = sys.argv[2]
u_k_path = sys.argv[3]#输出

with mrcfile.open(x_k_path) as mrc:
    x_k_data = mrc.data
with mrcfile.open(d_k_1_path) as mrc:
    d_k_1_data = mrc.data


#第一次没有d_k-1 uk=prox(xk)
#第二次：
uk=x_k_data-d_k_1_data

#写出
# 读取 x_k 和 d_k_1 的数据
with mrcfile.open(x_k_path, mode='r') as mrc:
    x_k_data = mrc.data.copy()  # 使用 copy 确保从文件中读取数据

# 创建一个新的MRC文件
with mrcfile.new(u_k_path, overwrite=True) as new_mrc:
    # 使用 x_k 的 header 作为基础
    new_mrc.set_data(uk.astype(np.float32))
    new_mrc.header = mrc.header  # 复制原始文件的header信息

    # 更新统计信息
    new_mrc.header.dmin = np.min(uk)
    new_mrc.header.dmax = np.max(uk)
    new_mrc.header.dmean = np.mean(uk)

# 这里 new_mrc_path 就是包含处理后数据的新MRC文件的路径