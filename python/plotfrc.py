import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 文件路径
path1 = "/data/论文比较/frc/BBb/SIRT_BBb/imod.txt"
path4 = "/data/论文比较/frc/BBb/SIRT_BBb/Ours(slice-based)).txt"
path3 = "/data/论文比较/frc/BBb/SIRT_BBb/Ours(unified memory).txt"
path2 = "/data/论文比较/frc/BBb/SIRT_BBb/tomo3d.txt"
path5 = "/data/重建数据集/isonet/result/4.txt"
file_paths = [path1,path4,path3,path2]

data = {}

# 读取每个文件的数据+
for path in file_paths:
    file_name = path.split('/')[-1]
    data[file_name] = pd.read_csv(path, delimiter='\t', header=None)

plt.figure(figsize=(15, 10))

target_frc_value = 0.143

# for file_name, df in data.items():
#     # 使用线性插值找到目标FRC值对应的分数
#     fraction_at_target_frc = np.interp(target_frc_value, df[0], df[1])
#     print(f"{file_name}: Fraction at FRC 0.143 is {fraction_at_target_frc:.3f}")

# 绘制每个数据集的FRC曲线
for file_name, df in data.items():
    plt.plot(df[0], df[1], label=file_name.replace('.txt', ''))


plt.axhline(y=target_frc_value, color='orange', linestyle='--')
plt.xlabel('Fraction of Nyquist')
plt.ylabel('Fourier Ring Correlation')
plt.legend()
plt.grid(True)

# 保存为矢量图（SVG 格式）
plt.savefig('frc_data_comparison.svg', format='svg')
plt.show()