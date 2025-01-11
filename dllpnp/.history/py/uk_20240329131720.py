import mrcfile
import numpy as np

#input
x_k_path=""
d_k_1_path=""
#output
u_k_path=""

with mrcfile.open(x_k_path) as mrc:
    x_k_data = mrc.data
with mrcfile.open(d_k_1_path) as mrc:
    d_k_1_data = mrc.data


#第一次没有d_k-1 uk=prox(xk)
#第二次：
uk=x_k_data-d_k_1_data

#写出