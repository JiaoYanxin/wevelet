#!/bin/bash

# 检查是否有足够的参数传递给脚本
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 input_file [optional_second_parameter]"
    exit 1
fi

# 使用 $1 来获取命令行传递的第一个参数（输入文件）
input=$1


#提取出文件名 比如/home/a.txt 提取出a
file_name_without_extension=$(basename -s .txt "$input")

# 使用 $2 作为可选的第二个参数（未在你的示例中使用，但可保留）
fscname=$2

# 对输入文件进行重命名，添加前缀 "re_"
rename="./data/BBa/result/pnp/re_$file_name_without_extension"

# 打印新的文件名用于调试
echo "Re_named file: $rename"

# 你的 xyzproj 命令（根据需要调整）
xyzproj -input $input -output $rename -angles 1,1,0 -axis y -ray

# 打印完成信息
echo "处理完成，输出文件为：$re_name"

#计算fsc 需要eman2环境
e2proc3d.py /home/jiaoyx/singlexsart/pnp/data/BBa/omit34/34.mrc /home/jiaoyx/singlexsart/python/FRC-main/fsc/$2.txt --calcfsc $rename

echo "FSC计算完成，输出文件为：$2.txt"
