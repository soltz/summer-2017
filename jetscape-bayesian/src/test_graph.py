from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# data_src = './data/Pyhon_nobg_data.txt'
# name = 'nobg_py'
# mult_i = 1
# xgj_i = 2
# quench_i = 3

# data_src = './data/C_nobg_data.txt'
# name = 'nobg'
# data_src = './data/C_bg_data.txt'
# name = 'bg'
# data_src = './data/C_bgsub_data.txt'
# name = 'bgsub'
data_src = './data/C_nobg_orig_data.txt'
name = 'nobg_orig'
mult_i = 1
xgj_i = 2
quench_i = 0

data = np.loadtxt(data_src)

filter_xgj_0 = True

bin_bounds = [0, 50, 100, 150, 200]
# bin_bounds = [0, 50, 100, 150, 200, 300, 400, 500, 600, 800, 1000]

bins = []
for i in range(len(bin_bounds) - 1):
    bins.append([])

for i in data:
    for j in range(len(bins)):
        if bin_bounds[j] < i[mult_i] <= bin_bounds[j + 1]:
            if not(filter_xgj_0 and i[xgj_i] == 0):
                bins[j].append(i)

bins = [np.array(i) for i in bins]
xgj_edges = list(np.arange(0, 2.1, 0.1))
quench_edges = list(np.arange(0.125, 1.225, 0.05))

for i in range(len(bins)):
    quench = bins[i][:,quench_i]
    xgj = bins[i][:,xgj_i]
    # plt.hist2d(quench, xgj, bins=40, norm=LogNorm())
    # plt.colorbar()
    # plt.show()
    H, xedges, yedges = np.histogram2d(quench, xgj, bins=(quench_edges, xgj_edges))
    H = H.T
    X, Y = np.meshgrid(xedges, yedges)
    plt.pcolormesh(X, Y, H, norm=LogNorm())
    plt.colorbar()
    plt.xlabel('quench')
    plt.ylabel('xgj')
    plt.title('2-D histogram with multiplicity between {} and {}'.format(bin_bounds[i], bin_bounds[i + 1]))
    plt.savefig('graphs/{}_xgj_quench_2dhist_mult{}_{}.pdf'.format(name, bin_bounds[i], bin_bounds[i + 1]))
    plt.clf()
    # plt.show()

    xgj_1_0 = []
    xgj_0_2 = []
    xgj_0_5 = []

    for j in range(len(quench)):
        if quench[j] == 1.0:
            xgj_1_0.append(xgj[j])
        elif quench[j] == 0.2:
            # pass
            xgj_0_2.append(xgj[j])
        elif quench[j] == 0.5:
            xgj_0_5.append(xgj[j])


    # Quench 1.0 histogram
    plt.hist(xgj_1_0, 40, normed=1)
    plt.xlabel('xgj')
    plt.ylabel('Probability')
    plt.title('Histogram for quench=1.0 and multiplicity between {} and {}'.format(bin_bounds[i], bin_bounds[i + 1]))
    plt.axis([0, 2.0, 0, 5.0])
    plt.savefig('graphs/{}_xgj_hist_quench1.0_mult{}_{}.pdf'.format(name, bin_bounds[i], bin_bounds[i + 1]))
    plt.clf()
    # plt.show()

    # Quench 0.2 histogram
    plt.hist(xgj_0_2, 40, normed=1)
    plt.xlabel('xgj')
    plt.ylabel('Probability')
    plt.title('Histogram for quench=0.2 and multiplicity between {} and {}'.format(bin_bounds[i], bin_bounds[i + 1]))
    plt.axis([0, 2.0, 0, 5.0])
    plt.savefig('graphs/{}_xgj_hist_quench0.2_mult{}_{}.pdf'.format(name, bin_bounds[i], bin_bounds[i + 1]))
    plt.clf()
    # plt.show()

    # Quench 0.5 histogram
    plt.hist(xgj_0_5, 40, normed=1)
    plt.xlabel('xgj')
    plt.ylabel('Probability')
    plt.title('Histogram for quench=0.5 and multiplicity between {} and {}'.format(bin_bounds[i], bin_bounds[i + 1]))
    plt.axis([0, 2.0, 0, 5.0])
    plt.savefig('graphs/{}_xgj_hist_quench0.5_mult{}_{}.pdf'.format(name, bin_bounds[i], bin_bounds[i + 1]))
    plt.clf()
    # plt.show()
