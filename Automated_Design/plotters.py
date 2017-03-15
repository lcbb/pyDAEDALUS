import numpy as np
from matplotlib import pyplot as plt
from os import path


def plot_edge_length_distributions(fname_no_ply,
                                   scale_edge_length_PLY,
                                   rounded_edge_length_PLY,
                                   results_foldername):
    min_len_nt = min(rounded_edge_length_PLY)
    bins_for_hist = range(min_len_nt, max(rounded_edge_length_PLY) + 3, 1)
    bins_for_plotting = [x - 0.25 for x in bins_for_hist[:-1]]

    fig_31 = plt.figure(0, figsize=(8, 6))
    fig_31.clf()
    y31 = np.histogram(scale_edge_length_PLY, bins=bins_for_hist)[0]
    plt.bar(bins_for_plotting, y31, width=0.5)
    plt.xlim((min_len_nt - 0.5, max(rounded_edge_length_PLY) + 1.5))
    plt.ylim((0, max(y31) + 1))
    plt.title('Minimum edge length {} bp'.format(min_len_nt))
    plt.xlabel('Edge length (bp)')
    plt.ylabel('Number of edges')
    plt.xticks(bins_for_hist)

    fig_32 = plt.figure(1, figsize=(8, 6))
    fig_32.clf()
    y32 = np.histogram(rounded_edge_length_PLY, bins=bins_for_hist)[0]
    bins_a = [x - 0.125 for x in bins_for_plotting]
    bins_b = [x + 0.38 for x in bins_for_plotting]
    plt.bar(bins_a, y31, width=0.25, color='b')
    plt.bar(bins_b, y32, width=0.25, color='r')

    plt.xlim((min_len_nt - 0.5, max(rounded_edge_length_PLY) + 1.5))
    plt.ylim((0, max(max(y31) + 1, 1)))
    plt.title('Edges rounded to nearest 10.5 bp')
    plt.xlabel('Edge length (bp)')
    plt.ylabel('Number of edges')
    plt.xticks(bins_for_hist)

    shape_name = path.basename(path.normpath(fname_no_ply))
    shape_name_with_len = shape_name + '_{}_'.format(min_len_nt)
    base_filename = path.join(results_foldername, shape_name_with_len)
    fig_31.savefig(base_filename + 'min_edge_length_dist.png',
                   bbox_inches='tight')
    fig_32.savefig(base_filename + 'edges_rounded_to_10_5.png',
                   bbox_inches='tight')
