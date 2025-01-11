#!/opt/anaconda3-2019/bin/python
#######################################
# Author：xueyangcs
# Personal page: https://gitee.com/xueyangcs
# Last update: 2022-09-29 20:31:07
# Function:
# 用于比较两个mrc文件的每层差距
# Usage:
# 输入两个mrc的文件名，脚本将自动比较这两个文件的第0层并输出平均差值、最大差值、最小差值
# 随后脚本将询问比较第几层图像，并输出对应层的平均差值、最大差值、最小差值

#jyx:

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

# parser = argparse.ArgumentParser(description='Compare two tomogram')
# parser.add_argument('-i', '--input_tomo', type=str,
#                     default='BBb.rec', dest='input_tomo')
# parser.add_argument('-b', '--base_tomo', type=str,
#                     default='BBb.rec', dest='base_tomo')


def calc_error(slice1, slice2):
    diffMat = np.maximum(slice1-slice2, slice2-slice1)
    avgerror = np.sum(diffMat)/(slice1.shape[0]*slice1.shape[1])
    maxerror = np.max(diffMat)
    minerror = np.min(diffMat)
    return diffMat, avgerror, maxerror, minerror


def read_mrc_1slice(input_tomo, sliceIdx):
    with open(input_tomo, 'rb') as f:
        content = f.read()
        header1, data_start1, data_size = mrc.get_header(content)
        sliceSize = header1.nx * header1.ny * data_size
        f.seek(data_start1+sliceIdx * sliceSize)
        content = f.read(sliceSize)
        array = mrc.load_data(header1, content)
        maxdepth = header1.nz # 层数
    return array, maxdepth


def compare_tomo(input_tomo, base_tomo):
    compareSliceIdx = 0
    slice1, maxdepth = read_mrc_1slice(input_tomo, compareSliceIdx)
    slice2, maxdepth = read_mrc_1slice(base_tomo, compareSliceIdx)
    diffMat, avgerror, maxerror, minerror = calc_error(slice1, slice2)

    print("0 slice avg error=%.8f" % avgerror, " max error=%.8f" %
          maxerror, " min error=", minerror)
    while true:
        compareSliceIdx = int(input("compareSliceIdx="))
        if compareSliceIdx == -1:
            return
        if compareSliceIdx >= maxdepth:
            print("Out of depth!")
            continue
        slice1, maxdepth = read_mrc_1slice(input_tomo, compareSliceIdx)
        slice2, maxdepth = read_mrc_1slice(base_tomo, compareSliceIdx)
        diffMat, avgerror, maxerror, minerror = calc_error(slice1, slice2)
        print(compareSliceIdx, "slice avg error=%.8f" % avgerror, " max error=%.8f" %
              maxerror, " min error=", minerror)


def compare_tomo_jyx(input_tomo, base_tomo):
    compareSliceIdx = 0
    slice1, maxdepth = read_mrc_1slice(input_tomo, compareSliceIdx)
    slice2, maxdepth = read_mrc_1slice(base_tomo, compareSliceIdx)
    diffMat, avgerror, maxerror, minerror = calc_error(slice1, slice2)

    print("0 slice avg error=%.8f" % avgerror, " max error=%.8f" %
          maxerror, " min error=", minerror)
    while true:
        compareSliceIdx = int(input("compareSliceIdx="))
        if compareSliceIdx == -1:
            return
        if compareSliceIdx >= maxdepth:
            print("Out of depth!")
            continue
        slice1, maxdepth = read_mrc_1slice(input_tomo, compareSliceIdx)
        slice2, maxdepth = read_mrc_1slice(base_tomo, compareSliceIdx)
        diffMat, avgerror, maxerror, minerror = calc_error(slice1, slice2)

        sns.heatmap(diffMat, center=0, cmap=sns.diverging_palette(20,220,n=200))
        plt.show()
        print(compareSliceIdx, "slice avg error=%.8f" % avgerror, " max error=%.8f" %
              maxerror, " min error=", minerror)








# 输入两个参数分别是需要比较的两个文件
if __name__ == '__main__':
    base_tomo = sys.argv[2]
    tomo = sys.argv[1]
    #compare_tomo(input_tomo=tomo, base_tomo=base_tomo)
    compare_tomo_jyx(input_tomo=tomo, base_tomo=base_tomo)

