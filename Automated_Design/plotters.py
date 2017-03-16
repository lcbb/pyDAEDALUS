import numpy as np
from matplotlib import pyplot as plt
from os import path


def plot_edge_length_distributions(shape_name,
                                   scale_edge_length_PLY,
                                   rounded_edge_length_PLY,
                                   results_foldername):
    min_len_nt = min(rounded_edge_length_PLY)
    max_len_nt = max(rounded_edge_length_PLY)
    bins_for_hist = range(min_len_nt, max_len_nt + 3, 1)
    bins_for_plotting = [x - 0.25 for x in bins_for_hist[:-1]]

    # Everything pertaining to `width_multiplier` is to help with the shapes
    # with a very large difference between the largest edge lengths and the
    # smallest edge lengths.  The general idea is you accept the possibility
    # of overlapping bins (if there happen to be a few very similarly lengthed
    # edges) to make the bins unfaithfully larger (that is, large enough to
    # see on the graph).  This is all because the actual-width (as opposed to
    # the width chosen to show the bins as) of the bins is important to keep
    # at 1.
    bin_width = 0.25  # for rounded bins, so half of non-rounded bins
    bin_offset = float(bin_width)/2
    width_multiplier = (1 + (max_len_nt - min_len_nt) / 75)


    fig_31 = plt.figure(0, figsize=(16, 8))
    fig_31.clf()
    y31 = np.histogram(scale_edge_length_PLY, bins=bins_for_hist)[0]
    plt.bar([x + bin_width for x in bins_for_plotting],
            y31, width=2 * bin_width * width_multiplier)
    plt.xlim((min_len_nt - width_multiplier,
              max_len_nt + 1.5 * width_multiplier))
    plt.ylim((0, max(y31) + 1))
    plt.title('Minimum edge length {} bp'.format(min_len_nt))
    plt.xlabel('Edge length (bp)')
    plt.ylabel('Number of edges')
    # plt.xticks(bins_for_hist)

    fig_32 = plt.figure(1, figsize=(16, 8))
    fig_32.clf()
    y32 = np.histogram(rounded_edge_length_PLY, bins=bins_for_hist)[0]
    bins_a = [x - (bin_offset * width_multiplier)
              for x in bins_for_plotting]
    bins_b = [x + bin_width + (bin_offset * width_multiplier)
              for x in bins_for_plotting]
    plt.bar(bins_a, y31, width=bin_width * width_multiplier, color='b')
    plt.bar(bins_b, y32, width=bin_width * width_multiplier, color='r')

    plt.xlim((min_len_nt - width_multiplier,
              max_len_nt + 1.5 * width_multiplier))
    plt.ylim((0, max(max(y31) + 1, 1)))
    plt.title('Edges rounded to nearest 10.5 bp')
    plt.xlabel('Edge length (bp)')
    plt.ylabel('Number of edges')
    # plt.xticks(bins_for_hist)

    shape_name_with_len = shape_name + '_{}_'.format(min_len_nt)
    base_filename = path.join(results_foldername, shape_name_with_len)
    fig_31.savefig(base_filename + 'min_edge_length_dist.png',
                   bbox_inches='tight')
    fig_32.savefig(base_filename + 'edges_rounded_to_10_5.png',
                   bbox_inches='tight')
