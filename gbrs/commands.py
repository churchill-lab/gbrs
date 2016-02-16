import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
from itertools import combinations_with_replacement


def plot(gpbfile, outfile, sample_name=''):

    # TODO: Get the following from kwarg
    grid_size = 42586
    xt_max = 4501
    xt_size = 475
    width = 0.01
    chrs = map(str, np.arange(19) + 1) + ['X', 'Y', 'MT']
    haplotypes = ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')
    hid = dict(zip(haplotypes, np.arange(8)))
    cc_colors = ('yellow', 'gray', 'salmon', 'blue', 'dodgerblue', 'forestgreen', 'red', 'darkviolet')
    genotypes = np.array([h1+h2 for h1, h2 in combinations_with_replacement(haplotypes, 2)])

    #
    # Main body
    #
    genoprob = np.load(gpbfile)
    fig = pyplot.figure()
    fig.set_size_inches((16, 16))
    ax = fig.add_subplot(111)
    ax.set_xlim(0, 4500*width+width)
    ax.set_ylim(1, 95)
    num_recomb_total = 0
    for cid in xrange(len(chrs)):
        c = chrs[cid]
        if c in genoprob.files:  # Skip drawing Y chromosome if the sample is female
            genotype_calls = genotypes[genoprob[c].argmax(axis=0)]
            hap = []
            col1 = []
            col2 = []
            oldcol1 = 'NA'
            oldcol2 = 'NA'
            num_recomb = 0
            num_genes_in_chr = len(genotype_calls)
            for i in xrange(num_genes_in_chr):
                hap.append((i*width, width))
                c1 = cc_colors[hid[genotype_calls[i][0]]]
                c2 = cc_colors[hid[genotype_calls[i][1]]]
                if i > 0:
                    if c1 == c2:
                        if col1[-1] != col2[-1]: # When homozygous region starts, remember the most recent het
                            oldcol1 = col1[-1]
                            oldcol2 = col2[-1]
                    else:
                        if col1[-1] == col2[-1]: # When heterozygous region starts
                            if c1 == oldcol2 or c2 == oldcol1:
                                c1, c2 = c2, c1
                        elif c1 == col2[-1] or c2 == col1[-1]:
                            c1, c2 = c2, c1
                    if c1 != col1[-1] or c2 != col2[-1]:
                        num_recomb += 1
                col1.append(c1)
                col2.append(c2)
            num_recomb_total += num_recomb
            # plot
            ax.broken_barh(hap, (len(chrs)*4-cid*4+1, 1), facecolors=col1, edgecolor='face')
            ax.broken_barh(hap, (len(chrs)*4-cid*4, 1), facecolors=col2, edgecolor='face')
            ax.text((num_genes_in_chr+50)*width, len(chrs)*4-cid*4+0.5, '(%d)' % num_recomb)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_right()
        ax.get_yaxis().tick_left()
        ax.get_yaxis().set_ticks([])
        ax.set_yticklabels(chrs)
        pyplot.yticks(np.arange(len(chrs)*4+1, 1, -4), fontsize=14)
        ax.set_xticklabels([ '%dM' % xt for xt in np.arange(0, xt_max*grid_size/1000000, xt_size*grid_size/1000000)])
        pyplot.xticks(np.arange(0, xt_max*width, xt_size*width))
        title_txt = 'Genome reconstruction: ' + sample_name
        title_txt += "\n(Total %d recombinations)" % num_recomb_total
        ax.set_title(title_txt, fontsize=18, loc='center')
    fig.savefig(outfile, dpi=300)
    pyplot.close(fig)
