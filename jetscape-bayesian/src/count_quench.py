from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import argparse
import numpy as np

parser = argparse.ArgumentParser(
    description='Program to count events for each quench factor.',
    usage='%(prog)s FILENAME')
parser.add_argument('filename', metavar='FILENAME', type=str,
                    help='filename with xjg quench mult data')

# Parse arguments
args = parser.parse_args()

# Set up filename
filename = args.filename

# Load data
data = np.loadtxt(filename)
data = data.tolist()

# Count events for each quench
quench_dict = {}
for x in data:
    if x[0] in quench_dict:
        quench_dict[x[0]] += 1
    else:
        quench_dict[x[0]] = 1

# Print final dict
print(quench_dict)
