import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
from itertools import combinations_with_replacement
from collections import defaultdict
from scipy.interpolate import interp1d
from genome_info import CHRS, CHRLENS, NUM_CHRS, GENEPOS


def unit_vector(vector):
    return vector / np.linalg.norm(vector)


def get_genotype_probability(aln_profile, aln_specificity, sigma=0.12):
    # 'aln_specificity' should be a set of unit vectors (at least one of the entry is larger than 1.)
    num_haps = len(aln_profile)
    aln_vec = unit_vector(aln_profile)
    genoprob = []
    for i in xrange(num_haps):
        v1 = unit_vector(aln_specificity[i])
        for j in xrange(i, num_haps):
            if j == i:
                genoprob.append(sum(np.power(aln_vec - v1, 2))) # homozygotes
            else:
                v2 = unit_vector(aln_specificity[j])
                geno_vec = unit_vector(v1 + v2)
                # compute directional similarity
                genoprob.append(sum(np.power(aln_vec - geno_vec, 2))) # for heterozygotes
    genoprob = np.exp(np.array(genoprob) / (-2 * sigma * sigma))
    return np.array(genoprob / sum(genoprob))


def reconstruct(**kwargs):
    exprfile = kwargs.get('exprfile')
    avecfile = kwargs.get('avecfile')
    gidfile = kwargs.get('gidfile')
    tprobfile = kwargs.get('tprobfile')
    expr_threshold = kwargs.get('expr_threshold')
    sigma = kwargs.get('sigma')
    outdir = kwargs.get('outdir')

    # Load meta info and alignment specificity
    gid_genome_order = np.load(gidfile)
    avecs = np.load(open(avecfile))

    # Load expression level
    expr = dict()
    with open(exprfile) as fh:
        curline = fh.next()
        haplotypes = curline.rstrip().split("\t")[1:-1]
        num_haps = len(haplotypes)
        genotypes = [h1+h2 for h1, h2 in combinations_with_replacement(haplotypes, 2)]
        num_genos = len(genotypes)
        for curline in fh:
            item = curline.rstrip().split("\t")
            expr[item[0]] = np.array(map(float, item[1:-1]))

    # Get null model probability
    init_vec = []
    for h1, h2 in combinations_with_replacement(haplotypes, 2):
        if h1 == h2:
            init_vec.append(np.log(1.0/(num_haps*num_haps)))
        else:
            init_vec.append(np.log(2.0/(num_haps*num_haps)))
    init_vec = np.array(init_vec)

    # Get initial emission probability
    naiv_avecs = np.eye(num_haps) + (np.ones((num_haps, num_haps)) - np.eye(num_haps)) * 0.0001
    eprob = dict()
    for gid, evec in expr.iteritems():
        if sum(evec) < expr_threshold:
            eprob[gid] = init_vec
        elif gid not in avecs.files:
            eprob[gid] = np.log(get_genotype_probability(evec, naiv_avecs, sigma=0.450) + np.nextafter(0, 1))
        else:
            eprob[gid] = np.log(get_genotype_probability(evec, avecs[gid], sigma=sigma) + np.nextafter(0, 1))

    # Load transition probabilities
    tprob = np.load(tprobfile)

    # Get forward probability
    alpha = dict()
    alpha_scaler = dict()
    for c in CHRS:
        if c in tprob.files:
            tprob_c = tprob[c]
            gid_genome_order_c = gid_genome_order[c]
            num_genes_in_chr = len(gid_genome_order_c)
            alpha_c = np.zeros((num_genos, num_genes_in_chr))
            alpha_scaler_c = np.zeros(num_genes_in_chr)
            alpha_c[:, 0] = init_vec + eprob[gid_genome_order_c[0]]
            normalizer = np.log(sum(np.exp(alpha_c[:, 0])))
            alpha_c[:, 0] -= normalizer # normalization
            alpha_scaler_c[0] = -normalizer
            for i in xrange(1, num_genes_in_chr):
                alpha_c[:, i] = np.log(np.exp(alpha_c[:, i-1] + tprob_c[i-1]).sum(axis=1) + np.nextafter(0, 1)) + eprob[gid_genome_order_c[i]]
                normalizer = np.log(sum(np.exp(alpha_c[:, i])))
                alpha_c[:, i] -= normalizer  # normalization
                alpha_scaler_c[i] = -normalizer
            alpha[c] = alpha_c
            alpha_scaler[c] = alpha_scaler_c

    # Get backward probability
    beta = dict()
    for c in CHRS:
        if c in tprob.files:
            tprob_c = tprob[c]
            gid_genome_order_c = gid_genome_order[c]
            num_genes_in_chr = len(gid_genome_order_c)
            beta_c = np.zeros((num_genos, num_genes_in_chr))
            beta_c[:, -1] = alpha_scaler[c][-1]  #init_vec + eprob[gid_genome_order_c[-1]]
            for i in xrange(num_genes_in_chr-2, -1, -1):
                beta_c[:, i] = np.log(np.exp(tprob_c[i].transpose() +
                                             beta_c[:, i+1] +
                                             eprob[gid_genome_order_c[i+1]] +
                                             alpha_scaler[c][i]).sum(axis=1))
            beta[c] = beta_c

    # Get forward-backward probability
    gamma = dict()
    for c in CHRS:
        if c in tprob.files:
            gamma_c = np.exp(alpha[c] + beta[c])
            normalizer = gamma_c.sum(axis=0)
            gamma[c] = gamma_c / normalizer
    np.savez_compressed(os.path.join(outdir, 'gbrs.gamma.npz'), **gamma)

    # Run Viterbi
    delta = dict()
    for c in CHRS:
        if c in tprob.files:
            tprob_c = tprob[c]
            gid_genome_order_c = gid_genome_order[c]
            num_genes_in_chr = len(gid_genome_order_c)
            delta_c = np.zeros((num_genos, num_genes_in_chr))
            delta_c[:, 0] = init_vec + eprob[gid_genome_order_c[0]]
            for i in xrange(1, num_genes_in_chr):
                delta_c[:, i] = (delta_c[:, i-1] + tprob_c[i-1]).max(axis=1) + eprob[gid_genome_order_c[i]]
            delta[c] = delta_c
    viterbi_states = defaultdict(list)
    for c in CHRS:
        if c in tprob.files:
            tprob_c = tprob[c]
            gid_genome_order_c = gid_genome_order[c]
            num_genes_in_chr = len(gid_genome_order_c)
            sid = delta[c][:, num_genes_in_chr-1].argmax()
            viterbi_states[c].append(genotypes[sid])
            for i in reversed(xrange(num_genes_in_chr-1)):
                sid = (delta[c][:, i] + tprob_c[i][sid]).argmax()
                viterbi_states[c].append(genotypes[sid])
            viterbi_states[c].reverse()
    viterbi_states = dict(viterbi_states)
    np.savez_compressed(os.path.join(outdir, 'gbrs.genotypes.npz'), **viterbi_states)


def interpolate(**kwargs):
    gridfile = kwargs.get('gridfile')
    probfile = kwargs.get('probfile')
    outfile = kwargs.get('outfile')

    x_grid = defaultdict(list)
    with open(gridfile) as fh:
        fh.next()  # skip header
        for curline in fh:
            item = curline.rstrip().split("\t")
            x_grid[item[0]].append(float(item[1]))
    x_grid = dict(x_grid)

    x_grid_complete = dict()
    for c in CHRS:
        if c in GENEPOS.files:
            x = [float(coord) for m, coord in GENEPOS[c]]
            x_min = min(x_grid[c][0]-100.0, 0.0)
            x_max = max(x_grid[c][-1]+100.0, CHRLENS[c])
            x = np.append([x_min], x)
            x = np.append(x, [x_max])
            x_grid_complete[c] = x

    gamma_gene = np.load(probfile)
    gene_model_chr = {}
    gene_intrp_chr = {}
    for c in CHRS:
        if c in gamma_gene.files:
            gamma_gene_c = gamma_gene[c].transpose()
            y = np.append([gamma_gene_c[0, :]], gamma_gene_c, axis=0)
            y = np.append(y, [y[-1, :]], axis=0)
            gene_model_chr[c] = interp1d(x_grid_complete[c], y, axis=0)
            gene_intrp_chr[c] = gene_model_chr[c](x_grid[c])
    np.savez_compressed(outfile, **gene_intrp_chr)


def plot(**kwargs):
    gpbfile = kwargs.get('gpbfile')
    outfile = kwargs.get('outfile')
    sample_name = kwargs.get('sample_name')
    grid_size = kwargs.get('grid_size')
    xt_max = kwargs.get('xt_max')
    xt_size = kwargs.get('xt_size')
    width = kwargs.get('width')

    # TODO: Get the following from kwarg (or think about more elegant generalization)
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
    for cid in xrange(len(CHRS)):
        c = CHRS[cid]
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
            ax.broken_barh(hap, (NUM_CHRS*4-cid*4+1, 1), facecolors=col1, edgecolor='face')
            ax.broken_barh(hap, (NUM_CHRS*4-cid*4, 1), facecolors=col2, edgecolor='face')
            ax.text((num_genes_in_chr+50)*width, NUM_CHRS*4-cid*4+0.5, '(%d)' % num_recomb)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_right()
        ax.get_yaxis().tick_left()
        ax.get_yaxis().set_ticks([])
        ax.set_yticklabels(CHRS)
        pyplot.yticks(np.arange(NUM_CHRS*4+1, 1, -4), fontsize=14)
        ax.set_xticklabels([ '%dM' % xt for xt in np.arange(0, xt_max*grid_size/1000000, xt_size*grid_size/1000000)])
        pyplot.xticks(np.arange(0, xt_max*width, xt_size*width))
        title_txt = 'Genome reconstruction: ' + sample_name
        title_txt += "\n(Total %d recombinations)" % num_recomb_total
        ax.set_title(title_txt, fontsize=18, loc='center')
    fig.savefig(outfile, dpi=300)
    pyplot.close(fig)
