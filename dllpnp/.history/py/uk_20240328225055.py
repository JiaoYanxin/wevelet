import mrcfile
import numpy as np


u_k-1_path=""
with mrcfile.open(u_k-1_path) as mrc:
    uk_data = mrc.data
