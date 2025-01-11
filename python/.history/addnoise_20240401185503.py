import numpy as np
import mrcfile  # 假设你有一个用于读写 MRC 文件的库
import mrc, argparse, os

# 读取原始 MRC 文件
def read_mrc_file(file_path):
    with mrcfile.open(file_path, permissive=True) as mrc:
        data = mrc.data
    return data


# 写入带噪声的 MRC 文件
def write_mrc_file(input_file_path,file_path, data):
    # 转换数据类型为int16

    with open(input_file_path, 'rb') as f:
        content = f.read()
    tomo, header, _ = mrc.parse(content)

    header = header._replace(mode=2)  # 32-bit real
    header = header._replace(amin=tomo.min())
    header = header._replace(amax=tomo.max())
    header = header._replace(amean=tomo.mean())
  #  data_int16 = data.astype(np.int16)
   # with mrcfile.new(file_path, overwrite=True) as mrc:
  #      mrc.set_data(data_int16)
    with open(file_path, 'wb') as f:
        mrc.write(f, data)


# 添加高斯噪声
def add_gaussian_noise(images, SNR):
    images = images.astype(np.float32)
    # 计算每个图像的平均功率
    sig_power = np.sum(np.abs(images) ** 2, axis=(1, 2)) / (images.shape[1] * images.shape[2])

    # 根据给定的信噪比计算每个图像的噪声功率
    noise_power_db = 10 * np.log10(sig_power) - SNR

    # 计算每个图像的噪声功率
    noise_power = 10 ** (noise_power_db / 10)

    # 生成噪声并添加到每个图像上
    noisy_images = np.empty_like(images)
    for i in range(images.shape[0]):
        noise = np.random.normal(scale=np.sqrt(noise_power[i]), size=images[i].shape)
        noisy_images[i] = images[i] + noise

    return noisy_images


# 读取原始 MRC 文件
input_file_path="/data/重建数据集/deepdewadge/model0/grandmodel_noisefree.mrc"
output_file_path="/data/重建数据集/deepdewadge/model0/grandmodel_noise.mrc"
original_mrc_data = read_mrc_file(input_file_path)

# 添加高斯噪声
SNR = 0.5  # 信噪比为 0.5
noisy_mrc_data = add_gaussian_noise(original_mrc_data, SNR)

# 写入带噪声的 MRC 文件
write_mrc_file(input_file_path, output_file_path, noisy_mrc_data)
