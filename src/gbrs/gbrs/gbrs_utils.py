# standard library imports
from collections import defaultdict, OrderedDict
from itertools import combinations_with_replacement, product
import logging
import os
import re

# 3rd party library imports
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as pyplot
import numpy as np

# local library imports
from gbrs import utils

matplotlib.use('Agg')

DATA_DIR = os.getenv('GBRS_DATA', '.')

logging.getLogger('matplotlib').setLevel(logging.WARNING)
logger = utils.get_logger('gbrs')


def get_chromosome_info() -> OrderedDict[str, int]:
    fai_file = os.path.join(DATA_DIR, 'ref.fa.fai')
    try:
        chr_lens = OrderedDict(
            np.loadtxt(fai_file, usecols=(0, 1), dtype='|S8,<i4')
        )

        # convert from bytes to string
        chr_lens = OrderedDict({k.decode(): v for k, v in chr_lens.items()})
        return chr_lens
    except FileNotFoundError:
        raise ValueError(
            'Make sure if $GBRS_DATA is set correctly, and that "ref.fa.fai" is in that directory. Currently it is: '
            f'{DATA_DIR}'
        )


def get_founder_info():
    color_file = os.path.join(DATA_DIR, 'founder.hexcolor.info')

    try:
        founder_colors = OrderedDict(
            np.loadtxt(
                color_file,
                usecols=(0, 1),
                dtype='str',
                delimiter='\t',
                comments=None,
            )
        )

        return founder_colors
    except FileNotFoundError:
        raise ValueError(
            f'Make sure if $GBRS_DATA is set correctly, and that "founder.hexcolor.info" is in that directory. Currently it is: '
            f'{DATA_DIR}'
        )


def unit_vector(vector):
    if sum(vector) > 1e-6:
        return vector / np.linalg.norm(vector)
    else:
        return vector


def print_vecs(vecs, format_str: str = '%10.1f', show_sum: bool = False) -> None:
    for i in range(vecs.shape[0]):
        v = vecs[i]
        print(' '.join(format_str % elem for elem in v))
        if show_sum:
            print('\t=>', format_str % sum(v))
        else:
            print()


def get_genotype_probability(aln_profile, aln_specificity, sigma: float = 0.12):
    # 'aln_specificity' should be a set of unit vectors (at least one of the entry is larger than 1.)
    num_haps = len(aln_profile)
    aln_vec = unit_vector(aln_profile)
    genoprob = []
    for i in range(num_haps):
        v1 = unit_vector(aln_specificity[i])
        for j in range(i, num_haps):
            if j == i:
                genoprob.append(sum(np.power(aln_vec - v1, 2)))   # homozygotes
            else:
                v2 = unit_vector(aln_specificity[j])
                geno_vec = unit_vector(v1 + v2)
                # compute directional similarity
                genoprob.append(
                    sum(np.power(aln_vec - geno_vec, 2))
                )   # for heterozygotes
    genoprob = np.exp(np.array(genoprob) / (-2 * sigma * sigma))
    return np.array(genoprob / sum(genoprob))


def ris_step(
    gen_left,
    gen_right,
    rec_frac,
    haps=('A', 'B'),
    gamma_scale=0.1,
    is_x_chr=False,
    forward_direction=True,
):
    """
    Log transition probability for RIL by sib-mating
    Originally part of r/qtl2 designed/coded by Karl Broman (http://kbroman.org/qtl2/)
    Ported to python by Karl Broman (https://gist.github.com/kbroman/14984b40b0eab71e51891aceaabec850)
    Extended to open the possibility of heterogyzosity by KB Choi

    Args:
        gen_left: left genotype
        gen_right: right genotype
        rec_frac: interval distance (cM)
        haps: list of parental strain
        gamma_scale: amount we allow heterozygosity
        is_x_chr: whether it is 'X' chromosome
        forward_direction: direction of intervals

    Returns:
        log_e transition probability
    """
    it = combinations_with_replacement(haps, 2)
    diplotype = [f'{ht1}{ht2}' for ht1, ht2 in it]

    if is_x_chr:
        R = (2 * rec_frac) / (1.0 + 4.0 * rec_frac)
        gamma = R * gamma_scale
        if forward_direction:
            if gen_left == diplotype[0]:
                if gen_right == diplotype[0]:
                    return np.log(1.0 - R) - np.log(1 + gamma)
                elif gen_right == diplotype[1]:
                    return np.log(gamma) - np.log(1 + gamma)
                elif gen_right == diplotype[2]:
                    return np.log(R) - np.log(1 + gamma)
            elif gen_left == diplotype[1]:
                return np.log(1 / 3.0)
            if gen_left == diplotype[2]:
                if gen_right == diplotype[0]:
                    return np.log(2.0 * R) - np.log(1 + gamma)
                elif gen_right == diplotype[1]:
                    return np.log(gamma) - np.log(1 + gamma)
                elif gen_right == diplotype[2]:
                    return np.log(1.0 - 2.0 * R) - np.log(1 + gamma)
        else:  # backward direction
            if gen_left == diplotype[0]:
                if gen_right == diplotype[0]:
                    return np.log(1.0 - 2.0 * R) - np.log(1 + gamma)
                elif gen_right == diplotype[1]:
                    return np.log(gamma) - np.log(1 + gamma)
                elif gen_right == diplotype[2]:
                    return np.log(2.0 * R) - np.log(1 + gamma)
            elif gen_left == diplotype[1]:
                return np.log(1 / 3.0)
            elif gen_left == diplotype[2]:
                if gen_right == diplotype[0]:
                    return np.log(R) - np.log(1 + gamma)
                elif gen_right == diplotype[1]:
                    return np.log(gamma) - np.log(1 + gamma)
                elif gen_right == diplotype[2]:
                    return np.log(1.0 - R) - np.log(1 + gamma)

    else:  # autosome
        R = 4.0 * rec_frac / (1 + 6.0 * rec_frac)
        gamma = R * gamma_scale
        if gen_left == diplotype[0]:
            if gen_right == diplotype[0]:
                return np.log(1.0 - R) - np.log(1 + gamma)
            elif gen_right == diplotype[1]:
                return np.log(gamma) - np.log(1 + gamma)
            elif gen_right == diplotype[2]:
                return np.log(R) - np.log(1 + gamma)
        elif gen_left == diplotype[1]:
            return np.log(1 / 3.0)
        elif gen_left == diplotype[2]:
            if gen_right == diplotype[0]:
                return np.log(R) - np.log(1 + gamma)
            elif gen_right == diplotype[1]:
                return np.log(gamma) - np.log(1 + gamma)
            elif gen_right == diplotype[2]:
                return np.log(1.0 - R) - np.log(1 + gamma)


def f2_step(
    gen_left, gen_right, rec_frac, is_x_chr=False, forward_direction=True
):
    return NotImplementedError


def cc_step(
    gen_left, gen_right, rec_frac, is_x_chr=False, forward_direction=True
):
    return NotImplementedError


def do_step(
    gen_left, gen_right, rec_frac, is_x_chr=False, forward_direction=True
):
    return NotImplementedError


def get_transition_prob(
    marker_file: str,
    haplotypes: str = 'A,B',
    mating_scheme: str = 'RI',
    gamma_scale: float = 0.01,
    epsilon: float = 0.000001,
    output_file: str = 'tranprob.npz'
) -> None:
    """
    Calculate the transition probabilities.

    Args:
        marker_file: marker file
        haplotypes: the haplotypes
        mating_scheme: mating scheme to use
        gamma_scale: scale for gamma
        epsilon: epsilon value
        output_file: output file
    """
    logger.info(f'Marker File: {marker_file}')
    logger.info(f'Haplotypes: {haplotypes}')
    logger.info(f'Mating Scheme: {mating_scheme}')
    logger.info(f'Gama Scale: {gamma_scale}')
    logger.info(f'Epsilon: {epsilon}')
    logger.info(f'Output File: {output_file}')

    haplotypes = haplotypes.split(',')
    it = combinations_with_replacement(haplotypes, 2)
    diplotype = [f'{ht1}{ht2}' for ht1, ht2 in it]
    num_diplotypes = len(diplotype)
    diplotype_index = np.arange(num_diplotypes)

    locs_by_chro = defaultdict(list)
    gpos_by_chro = defaultdict(list)

    logger.info(f'Loading marker file: {marker_file}')
    with open(marker_file) as fh:
        for line in fh:
            item = line.rstrip().split('\t')
            locs_by_chro[item[1]].append((item[0], float(item[3])))
            gpos_by_chro[item[1]].append((item[0], int(item[2])))
    locs_by_chro = dict(locs_by_chro)
    gpos_by_chro = dict(gpos_by_chro)

    logger.info(f'Saving {os.path.join(DATA_DIR, "ref.gene_pos.ordered.npz")}')
    np.savez_compressed(
        os.path.join(DATA_DIR, 'ref.gene_pos.ordered.npz'), **gpos_by_chro
    )

    if mating_scheme == 'RI':
        step_func = ris_step
    elif mating_scheme == 'F2':
        step_func = f2_step
    elif mating_scheme == 'CC':
        step_func = cc_step
    elif mating_scheme == 'DO':
        step_func = do_step
    else:
        raise ValueError(f'Unknown mating scheme: {mating_scheme}')

    tprob = dict()
    for c in locs_by_chro.keys():
        logger.debug(f'Working on {c}')
        is_x_chr = c == 'X'
        pdiff = np.diff(np.array([e[1] for e in locs_by_chro[c]]))
        pdiff[pdiff < epsilon] = epsilon
        ndiff = len(pdiff)
        tprob[c] = np.ndarray(
            shape=(ndiff, num_diplotypes, num_diplotypes), dtype=float
        )
        for dt1id, dt2id in list(product(diplotype_index, repeat=2)):
            dt1 = diplotype[dt1id]
            dt2 = diplotype[dt2id]
            for i, d in enumerate(pdiff):
                tprob[c][i][dt1id, dt2id] = step_func(
                    dt1,
                    dt2,
                    d,
                    gamma_scale=gamma_scale,
                    haps=haplotypes,
                    is_x_chr=is_x_chr,
                )

    logger.info(f'Saving {os.path.join(DATA_DIR, output_file)}')
    np.savez_compressed(os.path.join(DATA_DIR, output_file), **tprob)
    logger.info('Done')


def get_alignment_spec(
    sample_file: str,
    haplotypes: list[str],
    gene2transcript: str,
    out_dir: str,
    min_expr: float = 2.0
) -> None:
    """

    Args:
        sample_file:
        haplotypes:
        min_expr:

    Returns:

    """
    logger.info(f'Sample File: {sample_file}')
    logger.info(f'Haplotypes: {haplotypes}')
    logger.info(f'Min Expression: {min_expr}')
    logger.info(f'Gene2Transcript: {gene2transcript}')

    num_strains = len(haplotypes)

    gname = np.loadtxt(
        gene2transcript,
        usecols=(0,),
        dtype='str'
    )
    num_genes = len(gname)
    gid = dict(zip(gname, np.arange(num_genes)))

    logger.info(f'Loading {sample_file}')
    flist = defaultdict(list)
    with open(sample_file) as fh:
        for line in fh:
            item = line.rstrip().split('\t')
            flist[item[0]].append(item[1])
    flist = dict(flist)

    dset = dict()
    for st in haplotypes:
        dmat_strain = np.zeros((num_genes, num_strains))
        for tpmfile in flist[st]:
            logger.debug(f'Working on {tpmfile}')
            dmat_sample = np.zeros((num_genes, num_strains))
            if not os.path.isfile(tpmfile):
                print(f'File {tpmfile} does not exist.')
                continue
            with open(tpmfile) as fh:
                fh.readline()  # header
                for curline in fh:
                    item = curline.rstrip().split('\t')
                    if item[0] in gid:
                        row = gid[item[0]]
                        dmat_sample[row, :] = map(
                            float, item[1: (num_strains + 1)]
                        )
            dmat_strain += dmat_sample
        dset[st] = dmat_strain / len(flist[st])

    axes = dict()
    ases = dict()
    avecs = dict()
    for g in gname:
        axes[g] = np.zeros((num_strains, num_strains))
        ases[g] = np.zeros((1, num_strains))
        good = np.zeros(num_strains)
        for i, st in enumerate(haplotypes):
            v = dset[st][gid[g], :]
            axes[g][i, :] = v
            ases[g][0, i] = sum(v)
            if sum(v) > min_expr:
                good[i] = 1.0
        if sum(good) > 0:  # At least one strain expresses
            avecs[g] = np.zeros((num_strains, num_strains))
            for i in range(num_strains):
                avecs[g][i, :] = unit_vector(axes[g][i, :])

    logger.info(f'Saving {os.path.join(out_dir, "axes.npz")}')
    np.savez_compressed(os.path.join(out_dir, 'axes.npz'), **axes)
    logger.info(f'Saving {os.path.join(out_dir, "ases.npz")}')
    np.savez_compressed(os.path.join(out_dir, 'ases.npz'), **ases)
    logger.info(f'Saving {os.path.join(out_dir, "aves.npz")}')
    np.savez_compressed(os.path.join(out_dir, 'avecs.npz'), **avecs)


def reconstruct(
    expression_file: str,
    tprob_file: str,
    avec_file: str = None,
    gpos_file: str = None,
    expr_threshold: float = 1.5,
    sigma: float = 0.12,
    outbase: str = None,
) -> None:
    """
    Reconstruct the genome based upon gene-level TPM quantities.

    Args:
        expression_file: file containing gene-level TPM quantities
        tprob_file: transition probabilities file
        avec_file: alignment specificity file
        gpos_file: meta information for genes (chrom, id, location)
        expr_threshold:
        sigma:
        outbase: base output
    """
    if outbase is None:
        out_gtype = 'gbrs.reconstructed.genotypes.tsv'
        out_gprob = 'gbrs.reconstructed.genoprobs.npz'
    else:
        out_gtype = f'{outbase}.genotypes.tsv'
        out_gprob = f'{outbase}.genoprobs.npz'

    out_gtype_ordered = f'{os.path.splitext(out_gtype)[0]}.npz'

    if avec_file is None:
        avec_file = os.path.join(DATA_DIR, 'avecs.npz')

    if gpos_file is None:
        gpos_file = os.path.join(DATA_DIR, 'ref.gene_pos.ordered.npz')

    logger.info(f'Expression File: {expression_file}')
    logger.info(f'Transition Probabilities File: {tprob_file}')
    logger.info(f'Alignment Specificity File: {avec_file}')
    logger.info(f'Gene Position File: {gpos_file}')
    logger.info(f'Expression Threshold: {expr_threshold}')
    logger.info(f'Sigma: {sigma}')
    logger.info(f'Outbase: {outbase}')

    logger.info('Loading chromosome information')
    chrlens = get_chromosome_info()
    chrs = chrlens.keys()

    # Load alignment specificity
    logger.info(f'Loading alignment specificity: {avec_file}')
    avecs = np.load(avec_file)

    # Load meta info
    logger.info(f'Loading gene meta data: {gpos_file}')
    gene_pos = np.load(gpos_file)

    # gid_genome_order is a dictionary
    # with keys being chromosomes
    # and values being a list of gene ids and their position
    gid_genome_order = dict.fromkeys(gene_pos.files)
    for c in gene_pos.files:
        try:
            gid_genome_order[c] = np.array([g.decode() for g, p in gene_pos[c]])
        except Exception as e:
            logger.error(f"Error occurred: {e}")
            gid_genome_order[c] = np.array([g for g, p in gene_pos[c]])

    # Load expression level
    logger.info(f'Loading expression level data: {expression_file}')
    expr = dict()
    with open(expression_file) as fh:
        curline = fh.readline()
        haplotypes = curline.rstrip().split('\t')[1:-1]
        num_haps = len(haplotypes)
        genotypes = [
            h1 + h2 for h1, h2 in combinations_with_replacement(haplotypes, 2)
        ]
        num_genos = len(genotypes)
        for curline in fh:
            item = curline.rstrip().split('\t')
            expr[item[0]] = np.array(list(map(float, item[1:-1])))

    # Get null model probability
    logger.debug('Get null model probability')
    init_vec = []
    for h1, h2 in combinations_with_replacement(haplotypes, 2):
        if h1 == h2:
            init_vec.append(np.log(1.0 / (num_haps * num_haps)))
        else:
            init_vec.append(np.log(2.0 / (num_haps * num_haps)))
    init_vec = np.array(init_vec)

    # Get initial emission probability
    logger.debug('Get initial emission probability')
    naiv_avecs = (
        np.eye(num_haps)
        + (np.ones((num_haps, num_haps)) - np.eye(num_haps)) * 0.0001
    )
    eprob = dict()
    for gid, evec in expr.items():
        if sum(evec) < expr_threshold:
            eprob[gid] = init_vec
        elif gid not in avecs.files:
            eprob[gid] = np.log(
                get_genotype_probability(evec, naiv_avecs, sigma=0.450)
                + np.nextafter(0, 1)
            )
        else:
            eprob[gid] = np.log(
                get_genotype_probability(evec, avecs[gid], sigma=sigma)
                + np.nextafter(0, 1)
            )

    # Load transition probabilities
    logger.info(f'Loading transition probabilities: {tprob_file}')
    tprob = np.load(tprob_file)

    # Get forward probability
    logger.info('Getting forward probability')
    alpha = dict()
    alpha_scaler = dict()
    for c in chrs:
        if c in tprob.files:
            logger.debug(f'Working on {c}')
            tprob_c = tprob[c]
            gid_genome_order_c = gid_genome_order[c]
            num_genes_in_chr = len(gid_genome_order_c)
            alpha_c = np.zeros((num_genos, num_genes_in_chr))
            alpha_scaler_c = np.zeros(num_genes_in_chr)
            alpha_c[:, 0] = init_vec + eprob[gid_genome_order_c[0]]
            normalizer = np.log(sum(np.exp(alpha_c[:, 0])))
            alpha_c[:, 0] -= normalizer   # normalization
            alpha_scaler_c[0] = -normalizer
            for i in range(1, num_genes_in_chr):
                alpha_c[:, i] = (
                    np.log(
                        np.exp(alpha_c[:, i - 1] + tprob_c[i - 1]).sum(axis=1)
                        + np.nextafter(0, 1)
                    )
                    + eprob[gid_genome_order_c[i]]
                )
                normalizer = np.log(sum(np.exp(alpha_c[:, i])))
                alpha_c[:, i] -= normalizer  # normalization
                alpha_scaler_c[i] = -normalizer
            alpha[c] = alpha_c
            alpha_scaler[c] = alpha_scaler_c

    # Get backward probability
    logger.info('Getting backward probability')
    beta = dict()
    for c in chrs:
        if c in tprob.files:
            logger.debug(f'Working on {c}')
            tprob_c = tprob[c]
            gid_genome_order_c = gid_genome_order[c]
            num_genes_in_chr = len(gid_genome_order_c)
            beta_c = np.zeros((num_genos, num_genes_in_chr))
            beta_c[:, -1] = alpha_scaler[c][
                -1
            ]  # init_vec + eprob[gid_genome_order_c[-1]]
            for i in range(num_genes_in_chr - 2, -1, -1):
                beta_c[:, i] = np.log(
                    np.exp(
                        tprob_c[i].transpose()
                        + beta_c[:, i + 1]
                        + eprob[gid_genome_order_c[i + 1]]
                        + alpha_scaler[c][i]
                    ).sum(axis=1)
                )
            beta[c] = beta_c

    # Get forward-backward probability
    logger.info('Getting forward-backward probability')
    gamma = dict()
    for c in chrs:
        if c in tprob.files:
            logger.debug(f'Working on {c}')
            gamma_c = np.exp(alpha[c] + beta[c])
            normalizer = gamma_c.sum(axis=0)
            gamma[c] = gamma_c / normalizer

    logger.info(f'Saving Reconstructed Genotype Probabilities: {out_gprob}')
    np.savez_compressed(out_gprob, **gamma)

    # Run Viterbi
    logger.info('Running Viterbi')
    delta = dict()
    for c in chrs:
        if c in tprob.files:
            tprob_c = tprob[c]
            gid_genome_order_c = gid_genome_order[c]
            num_genes_in_chr = len(gid_genome_order_c)
            delta_c = np.zeros((num_genos, num_genes_in_chr))
            delta_c[:, 0] = init_vec + eprob[gid_genome_order_c[0]]
            for i in range(1, num_genes_in_chr):
                delta_c[:, i] = (delta_c[:, i - 1] + tprob_c[i - 1]).max(
                    axis=1
                ) + eprob[gid_genome_order_c[i]]
            delta[c] = delta_c
    viterbi_states = defaultdict(list)
    gtcall_g = dict()
    for c in chrs:
        if c in tprob.files:
            tprob_c = tprob[c]
            gid_genome_order_c = gid_genome_order[c]
            num_genes_in_chr = len(gid_genome_order_c)
            sid = delta[c][:, num_genes_in_chr - 1].argmax()
            viterbi_states[c].append(genotypes[sid])
            if (num_genes_in_chr > len(tprob_c)):
                num_genes_in_chr = len(tprob_c)
                # Above avoids IndexError for legacy GRCm38 ref files.
                # In those transprob files: num_genes_in_chr > tprob_c by 1.
            for i in reversed(range(num_genes_in_chr)):
                sid = (delta[c][:, i] + tprob_c[i][sid]).argmax()
                viterbi_states[c].append(genotypes[sid])
                gtcall_g[gid_genome_order_c[i]] = genotypes[sid]
            viterbi_states[c].reverse()

    viterbi_states = dict(viterbi_states)

    logger.info(f'Saving Reconstructed Genotypes: {out_gtype}')
    with open(out_gtype, 'w') as fhout:
        fhout.write('#Gene_ID\tDiplotype\n')
        for g in sorted(gtcall_g.keys()):
            fhout.write(f'{g}\t{gtcall_g[g]}\n')

    logger.info(f'Saving Reconstructed Ordered Genotypes: {out_gtype_ordered}')
    np.savez_compressed(out_gtype_ordered, **viterbi_states)
    logger.info('Done')


def interpolate(
    genoprob_file: str,
    grid_file: str = None,
    gpos_file: str = None,
    output_file: str = None,
) -> None:
    """
    Args:
        genoprob_file: the EMASE genotype probability file
        grid_file: grid file (i.e, ref.genome_grid.64k.txt)
        gpos_file: meta information for genes (chrom, id, location)
        output_file: the PDF file
    """
    if gpos_file is None:
        gpos_file = os.path.join(DATA_DIR, 'ref.gene_pos.ordered.npz')
        try:
            x_gene = np.load(gpos_file)
        except Exception as e:
            logger.error(
                f'{e}: Please make sure if $GBRS_DATA is set correctly: {DATA_DIR}'
            )
            raise
        else:
            pass
    else:
        x_gene = np.load(gpos_file)

    if grid_file is None:
        grid_file = os.path.join(DATA_DIR, 'ref.genome_grid.64k.txt')

    if output_file is None:
        output_file = f'gbrs.interpolated.{os.path.basename(genoprob_file)}'

    logger.info(f'Genotype Probability File: {genoprob_file}')
    logger.info(f'Grid File: {grid_file}')
    logger.info(f'Gene Position File: {gpos_file}')
    logger.info(f'Output File: {output_file}')

    logger.info('Loading chromosome information')
    chrlens = get_chromosome_info()
    chrs = chrlens.keys()
    logger.info(f'Chrom Lens: {chrlens}')
    logger.info(f'Chroms: {chrs}')

    logger.info(f'Loading grid file: {grid_file}')
    x_grid = defaultdict(list)
    with open(grid_file) as fh:
        fh.readline()  # skip header (Assuming there is just one line of header)
        for line in fh:
            item = line.rstrip().split('\t')
            x_grid[item[1]].append(
                float(item[3])
            )  # x_grid[chr] = [...positions in cM...]
    x_grid = dict(x_grid)

    x_gene_extended = (
        dict()
    )  # Adding end points in case we have to extrapolate at the 1st or last grid
    for c in x_grid.keys():
        if c in x_gene.files:
            logger.debug(f'Working on {c}')
            x = [float(coord) for m, coord in x_gene[c]]
            # x_min = min(x_grid[c][0]-100.0, 0.0)
            # x_max = max(x_grid[c][-1]+1.0, chrlens[c])
            # x = np.append([x_min], x)
            # x = np.append(x, [x_max])
            x = np.append([0.0], x)
            x = np.append(
                x, [x_grid[c][-1] + 1.0]
            )  # Do we have chromosome length in cM?
            x_gene_extended[c] = x

    logger.info(f'Loading GBRS genotype probability file: {genoprob_file}')
    gamma_gene = np.load(genoprob_file)
    gene_model_chr = dict()
    gene_intrp_chr = dict()
    for c in x_grid.keys():
        if c in gamma_gene.files:
            logger.debug(f'Working on {c}')
            gamma_gene_c = gamma_gene[c]
            y = np.hstack((gamma_gene_c[:, 0][:, np.newaxis], gamma_gene_c))
            y = np.hstack((y, y[:, -1][:, np.newaxis]))
            gene_model_chr[c] = interp1d(x_gene_extended[c], y, axis=1)
            gene_intrp_chr[c] = gene_model_chr[c](x_grid[c])

    logger.info(f'Saving interpolate probability file: {output_file}')
    np.savez_compressed(output_file, **gene_intrp_chr)
    logger.info('Done')


def combine():
    raise NotImplementedError


def plot(
    genoprob_file: str,
    output_file: str = None,
    output_format: str = 'pdf',
    sample_name: str = '',
    grid_size: int = 2,
    xt_max: int = 5000,
    xt_size: int = 500,
    grid_width: float = 0.01,
) -> None:
    """
    Plot a reconstructed genome.

    Args:
        genoprob_file: the EMASE genotype probability file
        output_file: the PDF file
        output_format: The file format, e.g. 'png', 'pdf', 'svg', ... The
            behavior when this is unset is documented under *fname*
        sample_name: the sample name
        grid_size: size of the grid (advanced)
        xt_max: size of xt (advanced)
        xt_size: max xt (advanced)
        grid_width: grid width (advanced)
    """
    if output_file is None:
        output_file = os.path.splitext(os.path.basename(genoprob_file))[0]
        output_file = f'gbrs.plotted.{output_file}.{output_format}'

    logger.info(f'Genotype Probabilities File: {genoprob_file}')
    logger.info(f'Output File: {output_file}')
    logger.info(f'Output File: {output_format}')
    logger.info(f'Sample Name: {sample_name}')
    logger.info(f'Grid Size: {grid_size}')
    logger.info(f'XT Max: {xt_max}')
    logger.info(f'XT Size: {xt_size}')
    logger.info(f'Grid Width: {grid_width}')

    logger.info('Loading chromosome information')
    chrlens = get_chromosome_info()

    logger.info('Loading founder colors')
    hcolors = get_founder_info()
    haplotypes = hcolors.keys()
    hid = dict(zip(haplotypes, np.arange(8)))
    logger.info(f'Haplotype IDs: {hid}')

    genotypes = np.array(
        [h1 + h2 for h1, h2 in combinations_with_replacement(haplotypes, 2)]
    )

    #
    # Main body
    #
    logger.info(f'Loading GBRS genotype probability file: {genoprob_file}')
    genoprob = np.load(genoprob_file)

    def natural_sort(list):

        def convert(text):
            return int(text) if text.isdigit() else text.lower()

        def alphanum_key(key):
            return [convert(c) for c in re.split('([0-9]+)', key)]

        return sorted(list, key=alphanum_key)

    chrs = [value for value in chrlens.keys() if value in genoprob.files]
    # intersection of ref.fa.fai chroms and those present in genoprob.

    chrs = natural_sort(chrs)
    # natural sorting of the chrom list for display.

    num_chrs = len(chrs)

    fig = pyplot.figure()
    fig.set_size_inches((16, 16))
    ax = fig.add_subplot(111)
    ax.set_xlim(0, xt_max * grid_width + grid_width)
    ax.set_ylim(1, 95)
    num_recomb_total = 0
    for cid, c in enumerate(chrs):
        if (
            c in genoprob.files
        ):  # Skip drawing Y chromosome if the sample is female
            logger.debug(f'Working on {c}')
            genotype_calls = genotypes[genoprob[c].argmax(axis=0)]
            hap = []
            col1 = []
            col2 = []
            oldcol1 = 'NA'
            oldcol2 = 'NA'
            num_recomb = 0
            num_genes_in_chr = len(genotype_calls)
            for i in range(num_genes_in_chr):
                hap.append((i * grid_width, grid_width))
                c1 = hcolors[genotype_calls[i][0]]
                c2 = hcolors[genotype_calls[i][1]]
                if i > 0:
                    if c1 == c2:
                        if (
                            col1[-1] != col2[-1]
                        ):  # When homozygous region starts, remember the
                            # most recent het
                            oldcol1 = col1[-1]
                            oldcol2 = col2[-1]
                    else:
                        if (
                            col1[-1] == col2[-1]
                        ):  # When heterozygous region starts
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
            ax.broken_barh(
                hap,
                (num_chrs * 4 - cid * 4 + 1, 1),
                facecolors=col1,
                edgecolor='face',
            )
            ax.broken_barh(
                hap,
                (num_chrs * 4 - cid * 4, 1),
                facecolors=col2,
                edgecolor='face',
            )
            ax.text(
                (num_genes_in_chr + 50) * grid_width,
                num_chrs * 4 - cid * 4 + 0.5,
                f'({num_recomb})',
            )
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_right()
        ax.get_yaxis().tick_left()
        pyplot.yticks(
            ticks=np.arange(num_chrs * 4 + 1, 1, -4),
            labels=list(chrs),
            fontsize=14,
        )
        pyplot.xticks(
            ticks=np.arange(0, xt_max * grid_width, xt_size * grid_width),
            labels=[
                '%dcM' % xt
                for xt in np.arange(
                    0, xt_max * grid_size / 100, xt_size * grid_size / 100
                )
            ],
        )
        title_txt = f'Genome reconstruction: {sample_name}'
        title_txt += f'\n(Total {num_recomb_total} recombinations)'
        ax.set_title(title_txt, fontsize=18, loc='center')

    logger.info(f'Saving generated plot: {output_file}')
    fig.savefig(output_file, dpi=600, format=output_format)
    pyplot.close(fig)
    logger.info('Done')


def export(
    genoprob_file: str,
    strains: list[str],
    grid_file: str = None,
    output_file: str = None,
) -> None:
    """
    Export genotypes probability file in in GBRS quant format.

    Args:
        genoprob_file: genotype probability file
        strains: list of strains
        grid_file: grid file (i.e, ref.genome_grid.64k.txt)
        output_file: output file in GBRS quant format
    """
    if grid_file is None:
        grid_file = os.path.join(DATA_DIR, 'ref.genome_grid.64k.txt')

    if output_file is None:
        output_file = f'{os.path.splitext(genoprob_file)[0]}.tsv'

    logger.info(f'Genotype Probabilities File: {genoprob_file}')
    logger.info(f'Strains: {strains}')
    logger.info(f'Grid File: {grid_file}')
    logger.info(f'Output File: {output_file}')

    logger.info('Getting suffices for strains')
    num_strains = len(strains)
    hid = dict(zip(strains, np.arange(num_strains)))
    genotypes = [
        h1 + h2 for h1, h2 in combinations_with_replacement(strains, 2)
    ]
    num_genotypes = len(genotypes)
    convmat = np.zeros((num_genotypes, num_strains))
    for g in range(num_genotypes):
        h1, h2 = genotypes[g]
        convmat[g, hid[h1]] += 1
        convmat[g, hid[h2]] += 1
    convmat *= 0.5

    logger.info(f'Loading grid file: {grid_file}')
    with open(grid_file) as fh:
        next(fh)
        grid = OrderedDict()
        for line in fh:
            item = line.rstrip().split('\t')
            chrom = item[1]
            if chrom in grid:
                grid[chrom] += 1
            else:
                grid[chrom] = 1
    num_grids = sum(grid.values())
    logger.debug(f'Number of grids: {num_grids}')

    logger.info(f'Loading GBRS genotype probability file: {genoprob_file}')
    gprob = np.load(genoprob_file)
    chromosomes = list(grid.keys())
    gprob_mat = gprob[chromosomes[0]].transpose()

    for c in chromosomes[1:]:
        logger.debug(f'Stacking {c}')
        gprob_mat = np.vstack((gprob_mat, gprob[c].transpose()))

    logger.info('Converting genotype probability')
    gprob_mat_converted = np.dot(gprob_mat, convmat)

    logger.info(f'Saving GBRS quant format: {output_file}')
    np.savetxt(
        output_file,
        gprob_mat_converted,
        fmt='%.6f',
        delimiter='\t',
        header='\t'.join(strains),
    )

    logger.info('Done')
