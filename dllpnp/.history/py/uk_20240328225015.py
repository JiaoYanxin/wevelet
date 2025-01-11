import mrcfile

uk_path=""
with mrcfile.open(uk_path) as mrc:# 打开MRC文件
    uk_data = mrc.data#z x y