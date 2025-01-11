import mrc, argparse, os
import numpy as np

parser = argparse.ArgumentParser(description='Compute rotated tomogram')
parser.add_argument('-i', '--input_tomo', type=str, default='BBb.rec', dest='input_tomo')
parser.add_argument('-o', '--output_tomo', type=str, default='BBb.rec', dest='output_tomo')


def rotate_tomo(input_tomo, output_tomo):

    with open(input_tomo, 'rb') as f:
        content = f.read()
    tomo, header, _ = mrc.parse(content)
    tomo = tomo.astype(np.float32)
#旋转计算
    tomo = tomo.transpose(1,0,2)

    header = header._replace(mode=2)  # 32-bit real
    header = header._replace(amin=tomo.min())
    header = header._replace(amax=tomo.max())
    header = header._replace(amean=tomo.mean())

    with open(output_tomo, 'wb') as f:
        mrc.write(f, tomo)

if __name__ == '__main__':
    args = parser.parse_args()
    output_tomo = args.output_tomo
    tomo = args.input_tomo
    rotate_tomo(input_tomo=tomo, output_tomo=output_tomo)

