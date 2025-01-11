import mrcfile
import numpy as np


x_k_path=""
d_k_1_path=""
with mrcfile.open(x_k_path) as mrc:
    x_k_data = mrc.data
with mrcfile.open(d_k_1_path) as mrc:
    d_k_1_data = mrc.data


#第一次没有d_k-1 uk=prox(xk)
uk=x_k_data-d_k_1_data

#写出