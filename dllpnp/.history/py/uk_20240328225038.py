import mrcfile
import numpy as np


uk_path=""
with mrcfile.open(uk_path) as mrc:
    uk_data = mrc.data
