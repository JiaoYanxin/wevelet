#######################################
from sqlalchemy import true
import mrc
# import argparse
import numpy as np
import seaborn as sns
import pandas as pd
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from past.utils import old_div
from builtins import range
from EMAN2 import *
import sys
from math import *
from PIL import Image
# get a test model
a = test_image_3d(1)
# alternatively load yours from disk
a = EMData("/home/jiaoyx/singlexsart/pnp/data/BBa/result/ICON_omit33/mask_ICONreconstruction_r.mrc")
# 指定投影角度（这里假设角度为45度）
projection_angle = 1.0

# 创建一个旋转变换对象，将投影角度应用于 z 轴
transform = Transform({"type": "xyz", "ytilt": projection_angle})
# make a projection along the z axis
proj = a.project("standard", transform)
# Another way of making a projection..
save_path = "/home/jiaoyx/singlexsart/python/FRC-main/fsc/proj_image.mrc"
proj.write_image(save_path)
display(proj)

