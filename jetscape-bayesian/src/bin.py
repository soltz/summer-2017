from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import argparse
import numpy as np

# Set up argument parser
parser = argparse.ArgumentParser(
    description='Program to bin data by multiplicity percentile.',
    usage='%(prog)s FILENAME [options]')
parser.add_argument('filename', metavar='FILENAME', type=str,
                    help='filename with xjg quench mult data')
parser.add_argument('-n', '--num_bins', type=int, default=5,
                    metavar='NUMBER_OF_BINS',
                    help='number of bins into which to divide data')

# Parse arguments
args = parser.parse_args()

# Set up filename patterns
filename = args.filename
prefix = '.'.join(filename.split('.')[0:-1])
ext = filename.split('.')[-1]

# Load and sort data
data = np.loadtxt(filename)
data = data.tolist()
sorted_data = sorted(data, key=lambda x: x[1])

# Prepare for quick binning
size = len(sorted_data)
inc = 100 / args.num_bins
low = 0
high = inc
i = 0

while low < 100:
    low_i = int(low * size / 100)
    high_i = int(high * size / 100)
    sorted_data[low_i:high_i] = sorted(sorted_data[low_i:high_i], key=lambda x: x[0])
    with open('{}_{}to{}.{}'.format(prefix, low, high, ext), 'w') as f:
        while (i / size * 100 < high):
            f.write('{} {} {}\n'.format(sorted_data[i][0],
                                       int(sorted_data[i][1]),
                                       sorted_data[i][2]))
            i += 1
    low = high
    high += inc
