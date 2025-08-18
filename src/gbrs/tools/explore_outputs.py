import numpy as np
import pandas as pd
import sys
import os
import argparse
from collections import defaultdict
from itertools import combinations_with_replacement
from gbrs.gbrs import emase_utils
from gbrs.gbrs import gbrs_utils
from gbrs import utils
from pathlib import Path


def get_diplotype_order(haplotypes=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']):
    """Generate the diplotype order used by GBRS."""
    diplotypes = [h1 + h2 for h1, h2 in combinations_with_replacement(haplotypes, 2)]
    return diplotypes


logger = utils.configure_logging('gbrs', 5)
spacer = '-' * 40


def chromosome_sort_key(chrom: str | int) -> (int, str):
    chrom_str = str(chrom).lower().replace('chr', '')

    if chrom_str == 'x':
        return (23, 'x')
    elif chrom_str == 'y':
        return (24, 'y')
    elif chrom_str in ['m', 'mt']:
        return (25, chrom_str)
    elif chrom_str.isdigit():
        return (int(chrom_str), chrom_str)
    else:
        return (26, chrom_str)  # For anything else


def sort_chromosomes(chromosomes: list[str | int]) -> list[str | int]:
    return sorted(chromosomes, key=chromosome_sort_key)


def analyze_tpm_file(file_path: str | Path) -> None:
    """Debug TPM files (multiway or diploid TPM/read count files)."""
    file_path = Path(file_path).expanduser().resolve()

    if not file_path.exists():
        raise FileNotFoundError(file_path)

    logger.info(f'Anlyzing file: {file_path}')

    try:
        with pd.option_context('display.float_format',
                               '{:,.2f}'.format):  # Thousands separator and 2 decimal places
            # read the file
            df = pd.read_csv(file_path, delimiter='\t')
            columns = list(df.columns)

            logger.info(f'File shape: {df.shape}')
            logger.info(f'Columns: {columns}')

            # determine file type
            if columns[-1].lower() == 'notes':
                file_type = 'DIPLOID'
            else:
                file_type = 'MULTIWAY'

            logger.info(f'File type: {file_type}')

            # basic statistics
            logger.info(f'Basic statistics:')
            logger.info(f'\tTotal rows: {len(df)}')

            # analyze founder columns
            founder_cols = [col for col in columns if
                            col in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']]
            logger.info(f'\tFounder columns: {founder_cols}')

            # expression statistics
            if founder_cols:
                logger.info(f'Expression statistics by founder:')
                for founder in founder_cols:
                    values = df[founder]
                    non_zero = values[values > 0]
                    logger.info(
                        f'\t{founder}: mean={values.mean():.3f}, max={values.max():.3f}, non-zero={len(non_zero)} ({len(non_zero) / len(df) * 100:.1f}%)')

            # total expression
            if 'total' in columns:
                total_expr = df['total']
                logger.info(f'Total expression statistics:')
                logger.info(f'\tMean: {total_expr.mean():.3f}')
                logger.info(f'\tMedian: {total_expr.median():.3f}')
                logger.info(f'\tMax: {total_expr.max():.3f}')
                logger.info(
                    f'\tZero loci: {len(total_expr[total_expr == 0])} ({len(total_expr[total_expr == 0]) / len(df) * 100:.1f}%)')

            # genotype analysis (for diploid files)
            if file_type == 'DIPLOID':
                logger.info(f'Genotype analysis:')
                genotype_counts = df['notes'].value_counts()
                logger.info(f'\tUnique genotypes: {len(genotype_counts)}')
                logger.info(f'\tMost common genotypes (top 10):')
                for genotype, count in genotype_counts.head(10).items():
                    logger.info(f'\t{genotype}: {count} loci ({count / len(df) * 100:.1f}%)')

                # homozygous vs heterozygous
                homozygous = df['notes'].apply(
                    lambda x: x[0] == x[1] if isinstance(x, str) and len(x) == 2 else False)
                heterozygous = ~homozygous
                logger.info(
                    f'\tHomozygous: {homozygous.sum()} loci ({homozygous.sum() / len(df) * 100:.1f}%)')
                logger.info(
                    f'\tHeterozygous: {heterozygous.sum()} loci ({heterozygous.sum() / len(df) * 100:.1f}%)')

            # show first few entries
            logger.info(f'First 5 entries:')
            logger.info(df.head())

            # show last few entries
            logger.info(f'Last 5 entries:')
            logger.info(df.tail())

            # show first 3 loci
            #logger.info(f'First 3 loci:')
            #test_loci = list(df['locus'][:3])
            #for gene in test_loci:
            #    if gene in df['locus'].values:
            #        row = df[df['locus'] == gene].iloc[0]
            #        logger.info(f'Loci {gene}:')
            #        if file_type == 'DIPLOID':
            #            logger.info(f'\tGenotype: {row["notes"]}')
            #        if founder_cols:
            #            logger.info(f'\tExpression by founder:')
            #            founders_over_zero = 0
            #            for founder in founder_cols:
            #                if row[founder] > 0:
            #                    logger.info(f'\t\t{founder}: {row[founder]:.3f}')
            #                    founders_over_zero += 1
            #            if founders_over_zero == 0:
            #                logger.info('\t\tNo values over 0')

    except Exception as e:
        print(f'Error reading {file_path}: {e}')

analyze_expected_read_counts_file = analyze_tpm_file


def analyze_alignment_counts_file(file_path: str | Path) -> None:
    file_path = Path(file_path).expanduser().resolve()

    if not file_path.exists():
        raise FileNotFoundError(file_path)

    logger.info(f'Anlyzing file: {file_path}')

    df = pd.read_csv(file_path, sep='\t')

    # expect columns locus, aln_A…aln_H, uniq_A…uniq_H, locus_uniq
    expected_cols = ['locus']
    expected_cols += [f'aln_{c}' for c in 'ABCDEFGH']
    expected_cols += [f'uniq_{c}' for c in 'ABCDEFGH']
    expected_cols += ['locus_uniq']

    missing = [c for c in expected_cols if c not in df.columns]
    if missing:
        logger.error(f'Missing expected columns: {', '.join(missing)}')

    # show first few entries
    logger.info(f'First 5 entries:')
    logger.info(df.head())

    # show last few entries
    logger.info(f'Last 5 entries:')
    logger.info(df.tail())


def analyze_genoprobs_file(file_path: str | Path) -> None:
    file_path = Path(file_path).expanduser().resolve()

    if not file_path.exists():
        raise FileNotFoundError(file_path)

    logger.info(f'Anlyzing file: {file_path}')

    try:
        diplotypes = get_diplotype_order()
        data = np.load(file_path)

        chromosomes = sort_chromosomes(data.files)
        num_chromosomes = len(data.files)

        logger.info(f'Number of chromosomes: {num_chromosomes}')
        logger.info(f'Chromosomes: {', '.join(chromosomes)}')

        chromosome_details = {}
        total_positions = 0

        logger.info(f'Chromosome Details:')
        logger.info(f'{"Chr":<4} {"Loci":<10} {"Shape":<15}')

        for chrom in chromosomes:
            matrix = data[chrom]
            total_positions += matrix.shape[1]
            logger.info(f'{chrom:<4} {matrix.shape[1]:<10,} {str(matrix.shape):<15} ')

        logger.info(f'Total loci: {total_positions:,}')

        logger.info('Detailed Information:')

        for chrom in chromosomes:
            matrix = data[chrom]
            logger.info(f'Chromosome {chrom}:')
            logger.info(f'  Shape: {matrix.shape}')
            logger.info(f'  Data type: {matrix.dtype}')
            logger.info(f'  Number of loci: {matrix.shape[1]}')
            logger.info(f'  Number of diplotypes: {matrix.shape[0]}')

            # Cceck probability properties
            col_sums = matrix.sum(axis=0)
            logger.info(
                f'  Column sums (should be ~1.0): min={col_sums.min():.6f}, max={col_sums.max():.6f}, mean={col_sums.mean():.6f}')

            # find most likely diplotypes
            max_probs = matrix.max(axis=0)
            max_diplotype_indices = matrix.argmax(axis=0)

            # logger.debug(f'{max_probs=}')
            # logger.debug(f'{max_diplotype_indices=}')
            logger.info(
                f'  Max probabilities: min={max_probs.min():.6f}, max={max_probs.max():.6f}, mean={max_probs.mean():.6f}')

            # show distribution of most likely diplotypes
            unique_vals, first_idx, counts = np.unique(max_diplotype_indices, return_index=True,
                                                       return_counts=True)

            # sort primarily by counts (descending), then by first occurrence
            order = np.lexsort((first_idx, -counts))

            unique_sorted = unique_vals[order]
            counts_sorted = counts[order]

            for i, (diplotype_idx, count) in enumerate(zip(unique_sorted, counts_sorted)):
                # print(f'{i=}, {diplotype_idx=}, {count=}')
                if i < 10:
                    logger.info(
                        f'    {diplotypes[diplotype_idx]}: {count} loci ({count / matrix.shape[1] * 100:.1f}%)')

        # show detailed example for first chromosome
        if data.files:
            first_chrom = sorted(data.files)[0]
            matrix = data[first_chrom]
            logger.info(f'DETAILED EXAMPLE: Chromosome {first_chrom}')
            logger.info(f'First 5 loci, all diplotypes:')

            # show first 5 loci with their probabilities
            for locus_idx in range(min(5, matrix.shape[1])):
                logger.info(f'Locus {locus_idx}:')
                probs = matrix[:, locus_idx]
                max_idx = probs.argmax()
                logger.info(f'\tMost likely: {diplotypes[max_idx]} (prob={probs[max_idx]:.6f})')

                # show top 5 diplotypes
                top_indices = np.argsort(probs)[-5:][::-1]
                for i, idx in enumerate(top_indices):
                    logger.info(f'\t{i + 1}. {diplotypes[idx]}: {probs[idx]:.6f}')

    except Exception as e:
        logger.error(f'Error reading {file_path}: {e}')


def analyze_genotypes_file(file_path: str | Path) -> None:
    file_path = Path(file_path).expanduser().resolve()

    if not file_path.exists():
        raise FileNotFoundError(file_path)

    logger.info(f'Anlyzing file: {file_path}')

    try:
        diplotypes = get_diplotype_order()
        data = np.load(file_path)

        chromosomes = sort_chromosomes(data.files)
        num_chromosomes = len(data.files)

        logger.info(f'Number of chromosomes: {num_chromosomes}')
        logger.info(f'Chromosomes: {', '.join(chromosomes)}')

        chromosome_details = {}
        total_positions = 0
        total_homozygous_count = 0
        total_heterozygous_count = 0
        all_genotypes = []
        total_num_recombinations = 0
        total_recombinations = []

        logger.info('Chromosome Details:')
        logger.info(f'{"Chr":<4} {"Positions":<10} {"Shape":<15} {"Recombs":<5}')

        for chrom in chromosomes:
            genotypes = data[chrom]

            homozygous = np.array([g[0] == g[1] for g in genotypes])
            heterozygous = ~homozygous
            unique_genotypes, counts = np.unique(genotypes, return_counts=True)
            genotypes_list = genotypes.tolist()
            all_genotypes.extend(genotypes_list)

            # count recombinatiosns
            num_recombinations = 0
            recombinations = []
            tmp_g = genotypes_list[0]
            for i, g in enumerate(genotypes_list[1:]):
                if g != tmp_g:
                    num_recombinations += 1
                    recombinations.append(f'({i + 1}) {tmp_g} -> {g}')

                tmp_g = g

            total_num_recombinations += num_recombinations
            total_recombinations.extend(recombinations)

            logger.info(
                f'{chrom:<4} {genotypes.shape[0]:<10,} {str(genotypes.shape):<15} {num_recombinations:<5}')

            # summary['chromosome_details'][chrom] = {
            #    'shape': genotypes.shape,
            #    # position and genes are the same, just for ease of use
            #    'num_positions': genotypes.shape[0],
            #    'num_genes': genotypes.shape[0],
            #    'data_type': str(genotypes.dtype),
            #    'memory_mb': genotypes.nbytes / (1024 * 1024),
            #    'unique_genotypes': unique_genotypes,
            #    'unique_genotypes_counts': counts,
            #    'homozygous_count': homozygous.sum(),
            #    'heterozygous_count': heterozygous.sum(),
            #    'num_recombinations': num_recombinations,
            #    'recombinations': recombinations
            # }

            total_positions += genotypes.shape[0]
            total_homozygous_count += homozygous.sum()
            total_heterozygous_count += heterozygous.sum()

        logger.info(f'Total Positions: {total_positions}')
        logger.info(f'Total Homozygous Count: {total_homozygous_count}')
        logger.info(f'Total Heterozygous Count: {total_heterozygous_count}')
        logger.info(f'Total Recombinations: {total_num_recombinations}')

        # for r in total_recombinations:
        #    logger.info(r)
        # logger.info(f'All Genotypes: {all_genotypes}')

        # overall genotype distribution
        logger.info('\n' + '=' * 80)
        logger.info('OVERALL GENOTYPE DISTRIBUTION')
        logger.info('=' * 80)
        all_genotypes_array = np.array(all_genotypes)
        unique_all, counts_all = np.unique(all_genotypes_array, return_counts=True)
        logger.info(f'Total unique genotypes across all chromosomes: {len(unique_all)}')
        logger.info('Genotype distribution (sorted by frequency):')
        sorted_indices = np.argsort(counts_all)[::-1]
        for i, idx in enumerate(sorted_indices):
            genotype = unique_all[idx]
            count = counts_all[idx]
            percentage = count / len(all_genotypes_array) * 100
            logger.info(f"  {genotype}: {count} genes ({percentage:.1f}%)")
    except FileExistsError as e:
        logger.error(f'Error reading {file_path}: {e}')


def analyze_genotypes_file_tsv(file_path: str | Path) -> None:
    file_path = Path(file_path).expanduser().resolve()

    if not file_path.exists():
        raise FileNotFoundError(file_path)

    logger.info(f'Analyzing file: {file_path}')

    try:
        df = pd.read_csv(file_path, delimiter='\t')

        logger.info(f'File shape: {df.shape}')
        logger.info(f'Columns: {list(df.columns)}')

        # basic statistics
        logger.info(f'Basic statistics:')
        logger.info(f'\tTotal genes: {len(df)}')

        # diplotype distribution
        diplotype_counts = df['Diplotype'].value_counts()
        logger.info(f'Diplotype distribution (top 20):')
        for diplotype, count in diplotype_counts.head(20).items():
            logger.info(f'\t{diplotype}: {count} genes ({count / len(df) * 100:.1f}%)')

        # check for homozygous vs heterozygous
        homozygous = df['Diplotype'].apply(lambda x: x[0] == x[1])
        heterozygous = ~homozygous
        logger.info('Homozygous vs Heterozygous:')
        logger.info(
            f'\tHomozygous: {homozygous.sum()} genes ({homozygous.sum() / len(df) * 100:.1f}%)')
        logger.info(
            f'\tHeterozygous: {heterozygous.sum()} genes ({heterozygous.sum() / len(df) * 100:.1f}%)')

        # show first few entries
        logger.info('First 10 entries:')
        logger.info(df.head(10))

        # check for specific loci
        test_loci = list(df['#Gene_ID'][:3])
        for locus in test_loci:
            if locus in df['#Gene_ID'].values:
                row = df[df['#Gene_ID'] == locus].iloc[0]
                logger.info(f'\tLoci {locus}: {row["Diplotype"]}')

    except FileExistsError as e:
        logger.error(f'Error reading {file_path}: {e}')


