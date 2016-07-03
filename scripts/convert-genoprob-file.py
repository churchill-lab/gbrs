#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import logging
import numpy as np
from collections import OrderedDict
from itertools import combinations_with_replacement

DATA_DIR = os.getenv('GBRS_DATA', '.')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--genoprob-file', action='store', dest='gprobfile', type=str, required=True)
    parser.add_argument('-s', '--strains', action='store', dest='strains', type=str, required=True)
    parser.add_argument('-g', '--grid-file', action='store', dest='gridfile', type=str, default=None)
    parser.add_argument('-f', '--output-format', action='store', dest='outformat', type=str, default='tsv')
    parser.add_argument('-d', '--by-strains', action='store_false')
    parser.add_argument('-c', '--by-chromosomes', action='store_true')
    parser.add_argument("-o", action='store', dest='outfile', type=str, default=None)
    parser.add_argument("-v", "--verbose", help="Toggle DEBUG verbosity", action="store_true")
    return parser.parse_args()


def main(args, loglevel):

    logging.basicConfig(format="[gbrs::%(levelname)s] %(message)s", level=loglevel)
    logging.debug("Conversion requested: %s" % args.gprobfile)

    logging.info("Getting suffices for strains...")
    strains = args.strains.split(',')
    num_strains = len(strains)
    hid = dict(zip(strains, np.arange(num_strains)))
    genotypes = [ h1+h2 for h1, h2 in combinations_with_replacement(strains, 2) ]
    num_genotypes = len(genotypes)
    convmat = np.zeros((num_genotypes, num_strains))
    for g in xrange(num_genotypes):
        h1, h2 = genotypes[g]
        convmat[g, hid[h1]] += 1
        convmat[g, hid[h2]] += 1
    convmat *= 0.5

    if args.gridfile is None:
        args.gridfile = os.path.join(DATA_DIR, 'ref.genome_grid.64k.txt')
    with open(args.gridfile) as fh:
        fh.next()
        grid = OrderedDict()
        for curline in fh:
            item = curline.rstrip().split('\t')
            chro = item[0]
            if grid.has_key(chro):
                grid[chro] += 1
            else:
                grid[chro] = 1
    num_grids = sum(grid.values())
    logging.info("There were %d grids." % num_grids)

    logging.info("Loading Reading in GBRS genotype probability file...")
    gprob = np.load(args.gprobfile)
    chromosomes = grid.keys()
    gprob_mat = gprob[chromosomes[0]].transpose()
    for c in chromosomes[1:]:
        gprob_mat = np.vstack((gprob_mat, gprob[c].transpose()))

    logging.info("Converting genotype probability...")
    gprob_mat_converted = np.dot(gprob_mat, convmat)

    logging.info("Exporting to GBRS quant format...")
#    if args.outfile is None:
#        outfile = os.path.splitext(args.gprobfile)[0] + '.tsv'
#    else:
#        outfile = args.outbase
    if args.outfile is None:
        args.outfile = os.path.splitext(args.gprobfile)[0] + '.tsv'
    np.savetxt(args.outfile, gprob_mat_converted, fmt='%.6f', delimiter='\t', header='\t'.join(strains))

    logging.info("DONE")


if __name__ == '__main__':
    arguments = parse_args()
    # Setup logging
    if arguments.verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO
    sys.exit(main(arguments, log_level))
