import os
import sys
import numpy as np
from itertools import combinations_with_replacement, product
from collections import defaultdict, OrderedDict
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot


DATA_DIR = os.getenv('GBRS_DATA', '.')


def get_chromosome_info(caller="gbrs"):
    faifile = os.path.join(DATA_DIR, "ref.fa.fai")
    try:
        chrlens = OrderedDict(np.loadtxt(faifile, usecols=(0, 1), dtype="|S8,<i4"))
    except:
        print(f"[{caller}] Make sure if $GBRS_DATA is set correctly. Currently it is {DATA_DIR}", file=sys.stderr)
        raise
    else:
        return chrlens


def get_founder_info(caller="gbrs"):
    fcofile = os.path.join(DATA_DIR, "founder.hexcolor.info")
    try:
        fcolors = OrderedDict(np.loadtxt(fcofile, usecols=(0, 1), dtype="string", delimiter="\t", comments=None))
    except:
        print(f"[{caller}] Make sure if $GBRS_DATA is set correctly. Currently it is {DATA_DIR}", file=sys.stderr)
        raise
    else:
        return fcolors


def unit_vector(vector):
    if sum(vector) > 1e-6:
        return vector / np.linalg.norm(vector)
    else:
        return vector


def print_vecs(vecs, format_str="%10.1f", show_sum=False):
    for i in range(vecs.shape[0]):
        v = vecs[i]
        print(" ".join( format_str % elem for elem in v ))
        if show_sum:
            print("\t=>", format_str % sum(v))
        else:
            print


def get_genotype_probability(aln_profile, aln_specificity, sigma=0.12):
    # 'aln_specificity' should be a set of unit vectors (at least one of the entry is larger than 1.)
    num_haps = len(aln_profile)
    aln_vec = unit_vector(aln_profile)
    genoprob = []
    for i in range(num_haps):
        v1 = unit_vector(aln_specificity[i])
        for j in range(i, num_haps):
            if j == i:
                genoprob.append(sum(np.power(aln_vec - v1, 2))) # homozygotes
            else:
                v2 = unit_vector(aln_specificity[j])
                geno_vec = unit_vector(v1 + v2)
                # compute directional similarity
                genoprob.append(sum(np.power(aln_vec - geno_vec, 2))) # for heterozygotes
    genoprob = np.exp(np.array(genoprob) / (-2 * sigma * sigma))
    return np.array(genoprob / sum(genoprob))


def ris_step(gen_left, gen_right, rec_frac, haps=('A', 'B'), gamma_scale=0.1, is_x_chr=False, forward_direction=True):
    """
    Log transition probability for RIL by sib-mating
    Originally part of r/qtl2 designed/coded by Karl Broman (http://kbroman.org/qtl2/)
    Ported to python by Karl Broman (https://gist.github.com/kbroman/14984b40b0eab71e51891aceaabec850)
    Extended to open the possibility of heterogyzosity by KB Choi

    :param gen_left: left genotype
    :param gen_right: right genotype
    :param rec_frac: interval distance (cM)
    :param haps: list of parental strain
    :param gamma_scale: amount we allow heterozygosity
    :param is_x_chr: whether it is 'X' chromosome
    :param forward_direction: direction of intervals
    :return: log_e transition probability
    """
    it = combinations_with_replacement(haps, 2)
    diplotype = ["%s%s" % (ht1, ht2) for ht1, ht2 in it]

    if is_x_chr:
        R = (2*rec_frac)/(1.0 + 4.0*rec_frac)
        gamma = R * gamma_scale
        if forward_direction:
            if gen_left == diplotype[0]:
                if gen_right == diplotype[0]:
                    return np.log(1.0-R) - np.log(1+gamma)
                elif gen_right == diplotype[1]:
                    return np.log(gamma) - np.log(1+gamma)
                elif gen_right == diplotype[2]:
                    return np.log(R) - np.log(1+gamma)
            elif gen_left == diplotype[1]:
                return np.log(1/3.0)
            if gen_left == diplotype[2]:
                if gen_right == diplotype[0]:
                    return np.log(2.0*R) - np.log(1+gamma)
                elif gen_right == diplotype[1]:
                    return np.log(gamma) - np.log(1+gamma)
                elif gen_right == diplotype[2]:
                    return np.log(1.0-2.0*R) - np.log(1+gamma)
        else:  # backward direction
            if gen_left == diplotype[0]:
                if gen_right == diplotype[0]:
                    return np.log(1.0-2.0*R) - np.log(1+gamma)
                elif gen_right == diplotype[1]:
                    return np.log(gamma) - np.log(1+gamma)
                elif gen_right == diplotype[2]:
                    return np.log(2.0*R) - np.log(1+gamma)
            elif gen_left == diplotype[1]:
                return np.log(1/3.0)
            elif gen_left == diplotype[2]:
                if gen_right == diplotype[0]:
                    return np.log(R) - np.log(1+gamma)
                elif gen_right == diplotype[1]:
                    return np.log(gamma) - np.log(1+gamma)
                elif gen_right == diplotype[2]:
                    return np.log(1.0-R) - np.log(1+gamma)

    else:  # autosome
        R = 4.0*rec_frac/(1+6.0*rec_frac)
        gamma = R * gamma_scale
        if gen_left == diplotype[0]:
            if gen_right == diplotype[0]:
                return np.log(1.0-R) - np.log(1+gamma)
            elif gen_right == diplotype[1]:
                return np.log(gamma) - np.log(1+gamma)
            elif gen_right == diplotype[2]:
                return np.log(R) - np.log(1+gamma)
        elif gen_left == diplotype[1]:
            return np.log(1/3.0)
        elif gen_left == diplotype[2]:
            if gen_right == diplotype[0]:
                return np.log(R) - np.log(1+gamma)
            elif gen_right == diplotype[1]:
                return np.log(gamma) - np.log(1+gamma)
            elif gen_right == diplotype[2]:
                return np.log(1.0-R) - np.log(1+gamma)


def f2_step(gen_left, gen_right, rec_frac, is_x_chr=False, forward_direction=True):
    return NotImplementedError


def cc_step(gen_left, gen_right, rec_frac, is_x_chr=False, forward_direction=True):
    return NotImplementedError


def do_step(gen_left, gen_right, rec_frac, is_x_chr=False, forward_direction=True):
    return NotImplementedError


def get_transition_prob(**kwargs):
    haplotype = kwargs.get("haps").split(",")
    num_haplotypes = len(haplotype)
    it = combinations_with_replacement(haplotype, 2)
    diplotype = ["%s%s" % (ht1, ht2) for ht1, ht2 in it]
    num_diplotypes = len(diplotype)
    diplotype_index = np.arange(num_diplotypes)

    locs_by_chro = defaultdict(list)
    gpos_by_chro = defaultdict(list)
    with open(kwargs.get("mkrfile")) as fh:
        for curline in fh:
            item = curline.rstrip().split('\t')
            locs_by_chro[item[1]].append((item[0], float(item[3])))
            gpos_by_chro[item[1]].append((item[0], int(item[2])))
    locs_by_chro = dict(locs_by_chro)
    gpos_by_chro = dict(gpos_by_chro)
    np.savez_compressed(os.path.join(DATA_DIR, "ref.gene_pos.ordered.npz"), **gpos_by_chro)

    mating_scheme = kwargs.get("mating_scheme")
    if mating_scheme == "RI":
        stepfunc = ris_step
    elif mating_scheme == "F2":
        stepfunc = f2_step
    elif mating_scheme == "CC":
        stepfunc = cc_step
    elif mating_scheme == "DO":
        stepfunc = do_step

    gamma_scale = kwargs.get("gamma_scale")
    epsilon = kwargs.get("epsilon")
    tprob = dict()
    for c in locs_by_chro.keys():
        is_x_chr = (c == "X")
        pdiff = np.diff(np.array([e[1] for e in locs_by_chro[c]]))
        pdiff[pdiff < epsilon] = epsilon
        ndiff = len(pdiff)
        tprob[c] = np.ndarray(shape=(ndiff, num_diplotypes, num_diplotypes), dtype=float)
        for dt1id, dt2id in list(product(diplotype_index, repeat=2)):
            dt1 = diplotype[dt1id]
            dt2 = diplotype[dt2id]
            for i, d in enumerate(pdiff):
                tprob[c][i][dt1id, dt2id] = stepfunc(dt1, dt2, d, gamma_scale=gamma_scale, haps=haplotype, is_x_chr=is_x_chr)
    np.savez_compressed(os.path.join(DATA_DIR, kwargs.get("outfile")), **tprob)


def get_alignment_spec(**kwargs):
    strains = kwargs.get("haps").split(",")
    num_strains = len(strains)

    gname = np.loadtxt(os.path.join(DATA_DIR, "ref.gene2transcripts.tsv"), usecols=(0,), dtype="string")
    num_genes = len(gname)
    gid = dict(zip(gname, np.arange(num_genes)))

    flist = defaultdict(list)
    with open(kwargs.get("smpfile")) as fh:
        for curline in fh:
            item = curline.rstrip().split("\t")
            flist[item[0]].append(item[1])
    flist = dict(flist)

    dset = dict()
    for st in strains:
        dmat_strain = np.zeros((num_genes, num_strains))
        for tpmfile in flist[st]:
            dmat_sample = np.zeros((num_genes, num_strains))
            if not os.path.isfile(tpmfile):
                print(f"File {tpmfile} does not exist.")
                continue
            with open(tpmfile) as fh:
                fh.readline()  # header
                for curline in fh:
                    item = curline.rstrip().split("\t")
                    if item[0] in gid:
                        row = gid[item[0]]
                        dmat_sample[row, :] = map(float, item[1:(num_strains+1)])
            dmat_strain += dmat_sample
        dset[st] = dmat_strain / len(flist[st])
    min_expr = kwargs.get("min_expr")
    axes = dict()
    ases = dict()
    avecs = dict()
    for g in gname:
        axes[g] = np.zeros((num_strains, num_strains))
        ases[g] = np.zeros((1, num_strains))
        good = np.zeros(num_strains)
        for i, st in enumerate(strains):
            v = dset[st][gid[g], :]
            axes[g][i, :] = v
            ases[g][0, i] = sum(v)
            if sum(v) > min_expr:
                good[i] = 1.0
        if sum(good) > 0:  # At least one strain expresses
            avecs[g] = np.zeros((num_strains, num_strains))
            for i in range(num_strains):
                avecs[g][i, :] = unit_vector(axes[g][i, :])
    np.savez_compressed(os.path.join(DATA_DIR, "axes.npz"), **axes)
    np.savez_compressed(os.path.join(DATA_DIR, "ases.npz"), **ases)
    np.savez_compressed(os.path.join(DATA_DIR, "avecs.npz"), **avecs)


def reconstruct(**kwargs):
    chrlens = get_chromosome_info(caller="gbrs::reconstruct")
    chrs = chrlens.keys()

    outbase = kwargs.get("outbase")
    if outbase is None:
        out_gtype = "gbrs.reconstructed.genotypes.tsv"
        out_gprob = "gbrs.reconstructed.genoprobs.npz"
    else:
        out_gtype = outbase + ".genotypes.tsv"
        out_gprob = outbase + ".genoprobs.npz"
    out_gtype_ordered = os.path.splitext(out_gtype)[0] + ".npz"

    tprobfile = kwargs.get("tprobfile")
    exprfile = kwargs.get("exprfile")
    expr_threshold = kwargs.get("expr_threshold")
    sigma = kwargs.get("sigma")

    # Load alignment specificity
    avecfile = kwargs.get("avecfile")
    if avecfile is None:
        avecfile = os.path.join(DATA_DIR, "avecs.npz")
    avecs = np.load(avecfile)

    # Load meta info
    gposfile = kwargs.get("gposfile")
    if gposfile is None:
        gposfile = os.path.join(DATA_DIR, "ref.gene_pos.ordered.npz")
    gene_pos = np.load(gposfile)
    gid_genome_order = dict.fromkeys(gene_pos.files)
    for c in gene_pos.files:
        gid_genome_order[c] = np.array([g for g, p in gene_pos[c]])

    # Load expression level
    expr = dict()
    with open(exprfile) as fh:
        curline = fh.readline()
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
    for gid, evec in expr.items():
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
    for c in chrs:
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
            for i in range(1, num_genes_in_chr):
                alpha_c[:, i] = np.log(np.exp(alpha_c[:, i-1] + tprob_c[i-1]).sum(axis=1) + np.nextafter(0, 1)) + \
                                eprob[gid_genome_order_c[i]]
                normalizer = np.log(sum(np.exp(alpha_c[:, i])))
                alpha_c[:, i] -= normalizer  # normalization
                alpha_scaler_c[i] = -normalizer
            alpha[c] = alpha_c
            alpha_scaler[c] = alpha_scaler_c

    # Get backward probability
    beta = dict()
    for c in chrs:
        if c in tprob.files:
            tprob_c = tprob[c]
            gid_genome_order_c = gid_genome_order[c]
            num_genes_in_chr = len(gid_genome_order_c)
            beta_c = np.zeros((num_genos, num_genes_in_chr))
            beta_c[:, -1] = alpha_scaler[c][-1]  #init_vec + eprob[gid_genome_order_c[-1]]
            for i in range(num_genes_in_chr-2, -1, -1):
                beta_c[:, i] = np.log(np.exp(tprob_c[i].transpose() +
                                             beta_c[:, i+1] +
                                             eprob[gid_genome_order_c[i+1]] +
                                             alpha_scaler[c][i]).sum(axis=1))
            beta[c] = beta_c

    # Get forward-backward probability
    gamma = dict()
    for c in chrs:
        if c in tprob.files:
            gamma_c = np.exp(alpha[c] + beta[c])
            normalizer = gamma_c.sum(axis=0)
            gamma[c] = gamma_c / normalizer
    np.savez_compressed(out_gprob, **gamma)

    # Run Viterbi
    delta = dict()
    for c in chrs:
        if c in tprob.files:
            tprob_c = tprob[c]
            gid_genome_order_c = gid_genome_order[c]
            num_genes_in_chr = len(gid_genome_order_c)
            delta_c = np.zeros((num_genos, num_genes_in_chr))
            delta_c[:, 0] = init_vec + eprob[gid_genome_order_c[0]]
            for i in range(1, num_genes_in_chr):
                delta_c[:, i] = (delta_c[:, i-1] + tprob_c[i-1]).max(axis=1) + eprob[gid_genome_order_c[i]]
            delta[c] = delta_c
    viterbi_states = defaultdict(list)
    gtcall_g = dict()
    for c in chrs:
        if c in tprob.files:
            tprob_c = tprob[c]
            gid_genome_order_c = gid_genome_order[c]
            num_genes_in_chr = len(gid_genome_order_c)
            sid = delta[c][:, num_genes_in_chr-1].argmax()
            viterbi_states[c].append(genotypes[sid])
            for i in reversed(range(num_genes_in_chr-1)):
                sid = (delta[c][:, i] + tprob_c[i][sid]).argmax()
                viterbi_states[c].append(genotypes[sid])
                gtcall_g[gid_genome_order_c[i]] = genotypes[sid]
            viterbi_states[c].reverse()
    viterbi_states = dict(viterbi_states)
    np.savez_compressed(out_gtype_ordered, **viterbi_states)
    with open(out_gtype, "w") as fhout:
        fhout.write("#Gene_ID\tDiplotype\n")
        for g in sorted(gtcall_g.keys()):
            fhout.write(f"{g}\t{gtcall_g[g]}\n")


def interpolate(**kwargs):
    chrlens = get_chromosome_info(caller="gbrs::interpolate")
    chrs = chrlens.keys()

    probfile = kwargs.get("probfile")

    gposfile = kwargs.get("gposfile")
    if gposfile is None:  # if gposfile is not specified
        gposfile = os.path.join(DATA_DIR, "ref.gene_pos.ordered.npz")
        try:
            x_gene = np.load(gposfile)
        except:
            print(f"[gbrs::interpolate] Please make sure if $GBRS_DATA is set correctly: {DATA_DIR}", file=sys.stderr)
            raise
        else:
            pass
    else:
        x_gene = np.load(gposfile)

    gridfile = kwargs.get("gridfile")
    if gridfile is None:
        gridfile = os.path.join(DATA_DIR, "ref.genome_grid.64k.txt")

    outfile = kwargs.get("outfile")
    if outfile is None:
        outfile = "gbrs.interpolated." + os.path.basename(probfile)

    x_grid = defaultdict(list)
    with open(gridfile) as fh:
        fh.readline()  # skip header (Assuming there is just one line of header)
        for curline in fh:
            item = curline.rstrip().split("\t")
            x_grid[item[1]].append(float(item[3]))  # x_grid[chr] = [...positions in cM...]
    x_grid = dict(x_grid)

    x_gene_extended = dict()  # Adding end points in case we have to extrapolate at the 1st or last grid
    for c in x_grid.keys():
        if c in x_gene.files:
            x = [float(coord) for m, coord in x_gene[c]]
            #x_min = min(x_grid[c][0]-100.0, 0.0)
            #x_max = max(x_grid[c][-1]+1.0, chrlens[c])
            #x = np.append([x_min], x)
            #x = np.append(x, [x_max])
            x = np.append([0.0], x)
            x = np.append(x, [x_grid[c][-1]+1.0])  # Do we have chromosome length in cM?
            x_gene_extended[c] = x

    gamma_gene = np.load(probfile)
    gene_model_chr = dict()
    gene_intrp_chr = dict()
    for c in x_grid.keys():
        if c in gamma_gene.files:
            gamma_gene_c = gamma_gene[c]
            y = np.hstack((gamma_gene_c[:, 0][:, np.newaxis], gamma_gene_c))
            y = np.hstack((y, y[:, -1][:, np.newaxis]))
            gene_model_chr[c] = interp1d(x_gene_extended[c], y, axis=1)
            gene_intrp_chr[c] = gene_model_chr[c](x_grid[c])
    np.savez_compressed(outfile, **gene_intrp_chr)


def combine(**kwargs):
    raise NotImplementedError


def plot(**kwargs):
    chrlens = get_chromosome_info(caller="gbrs::plot")
    chrs = chrlens.keys()
    num_chrs = len(chrs)

    gpbfile = kwargs.get("gpbfile")
    outfile = kwargs.get("outfile")
    if outfile is None:
        outfile = "gbrs.plotted." + os.path.splitext(os.path.basename(gpbfile))[0] + ".pdf"

    sample_name = kwargs.get("sample_name")
    grid_size = kwargs.get("grid_size")
    xt_max = kwargs.get("xt_max")
    xt_size = kwargs.get("xt_size")
    width = kwargs.get("width")

    hcolors = get_founder_info(caller="gbrs::plot")
    haplotypes = hcolors.keys()
    hid = dict(zip(haplotypes, np.arange(8)))
    genotypes = np.array([h1+h2 for h1, h2 in combinations_with_replacement(haplotypes, 2)])

    #
    # Main body
    #
    genoprob = np.load(gpbfile)
    fig = pyplot.figure()
    fig.set_size_inches((16, 16))
    ax = fig.add_subplot(111)
    ax.set_xlim(0, xt_max*width+width)
    ax.set_ylim(1, 95)
    num_recomb_total = 0
    for cid, c in enumerate(chrs):
        if c in genoprob.files:  # Skip drawing Y chromosome if the sample is female
            genotype_calls = genotypes[genoprob[c].argmax(axis=0)]
            hap = []
            col1 = []
            col2 = []
            oldcol1 = "NA"
            oldcol2 = "NA"
            num_recomb = 0
            num_genes_in_chr = len(genotype_calls)
            for i in range(num_genes_in_chr):
                hap.append((i*width, width))
                c1 = hcolors[genotype_calls[i][0]]
                c2 = hcolors[genotype_calls[i][1]]
                if i > 0:
                    if c1 == c2:
                        if col1[-1] != col2[-1]:  # When homozygous region starts, remember the most recent het
                            oldcol1 = col1[-1]
                            oldcol2 = col2[-1]
                    else:
                        if col1[-1] == col2[-1]:  # When heterozygous region starts
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
            ax.broken_barh(hap, (num_chrs*4-cid*4+1, 1), facecolors=col1, edgecolor="face")
            ax.broken_barh(hap, (num_chrs*4-cid*4, 1), facecolors=col2, edgecolor="face")
            ax.text((num_genes_in_chr+50)*width, num_chrs*4-cid*4+0.5, f"({num_recomb})")
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_right()
        ax.get_yaxis().tick_left()
        ax.get_yaxis().set_ticks([])
        ax.set_yticklabels(chrs)
        pyplot.yticks(np.arange(num_chrs*4+1, 1, -4), fontsize=14)
        ax.set_xticklabels(["%dcM" % xt for xt in np.arange(0, xt_max*grid_size/100, xt_size*grid_size/100)])
        pyplot.xticks(np.arange(0, xt_max*width, xt_size*width))
        title_txt = f"Genome reconstruction: {sample_name}"
        title_txt += f"\n(Total {num_recomb_total} recombinations)"
        ax.set_title(title_txt, fontsize=18, loc='center')
    fig.savefig(outfile, dpi=300)
    pyplot.close(fig)
