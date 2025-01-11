import mrcfile
import numpy as np

x_k_path = " /home/jiaoyx/cryorecnet/0_xk.mrc"
out_path=" /home/jiaoyx/cryorecnet/zero.mrc"
with mrcfile.open(x_k_path, permissive=True) as mrc:
    x_k_data = mrc.data

out=np.zeros(x_k_data.shape)
with mrcfile.new(out_path, overwrite=True) as new_mrc:
    # 使用 x_k 的 header 作为基础
    new_mrc.set_data(out.astype(np.float32))
    #new_mrc.header = mrc.header  # 复制原始文件的header信息

    # 更新统计信息
    new_mrc.header.dmin = 0
    new_mrc.header.dmax = 0
    new_mrc.header.dmean = 0
