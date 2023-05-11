# standard library imports
from collections import defaultdict
from itertools import dropwhile
import gzip
import logging
import os
import string
import subprocess

# 3rd party library imports
logging.getLogger('numexpr').setLevel(logging.WARNING)
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np

# local library imports
from gbrs import utils
from gbrs.emase.AlignmentMatrixFactory import AlignmentMatrixFactory
from gbrs.emase.AlignmentPropertyMatrix import AlignmentPropertyMatrix
from gbrs.emase.EMfactory import EMfactory

logger = logging.getLogger('gbrs')


def bam2emase(
    alignment_file: str,
    haplotypes: list[str],
    locusid_file: str,
    output_file: str = 'alignments.transcriptome.h5',
    delim: str = '_',
    index_dtype: str = 'uint32',
    data_dtype: str = 'uint8'
) -> None:
    """
    Convert BAM file to EMASE format (hdf5).

    Args:
        alignment_file: The BAM file
        haplotypes: tuple or list of haplotypes
        locusid_file: filename for the locus (usually transcripts) info
        output_file: file name of the output file, if not specified one will be
            generated
        delim: delimiter string between locus and haplotype in BAM file
        index_dtype: data type of indices ptr, defaults to uint32
        data_dtype: data type of the stored value, defaults to uint8
    """
    logger.info(f'BAM File: {alignment_file}')
    logger.info(f'Locus ID File: {locusid_file}')
    logger.info(f'Output File: {output_file}')
    logger.info(f'Haplotypes: {haplotypes}')
    logger.info(f'Delimiter: {delim}')
    logger.info(f'Index dtype: {index_dtype}')
    logger.info(f'Data dtype: {data_dtype}')

    logger.info(f'Parsing Locus ID File: {locusid_file}')
    loci = utils.get_names(locusid_file)

    logger.info(f'Parsing BAM File: {alignment_file}')
    alignmat_factory = AlignmentMatrixFactory(alignment_file)
    alignmat_factory.prepare(
        haplotypes, loci, delim=delim, outdir=os.path.dirname(output_file)
    )

    logger.info(f'Saving EMASE Formatted File: {output_file}')
    alignmat_factory.produce(
        output_file, index_dtype=index_dtype, data_dtype=data_dtype
    )
    alignmat_factory.cleanup()
    logger.info('Done')


def combine(
    emase_files: list[str],
    output_file: str,
    comp_lib: str = 'zlib'
) -> None:
    """
    Combine EMASE files

    Args:
        emase_files: list of EMASE files to compress
        output_file: name of the combined EMASE file
        comp_lib: compression library to use
    """
    for x in emase_files:
        logger.info(f'EMASE file: {x}')
    logger.info(f'Output File: {output_file}')
    logger.info(f'Compression Library: {comp_lib}')

    aln_list = list()
    for f in emase_files:
        # TODO: Check if the loci and haplotypes are in the equivalent orders
        logger.info(f'Loading EMASE file: {f}')
        aln_mat = AlignmentPropertyMatrix(h5file=f)
        logger.debug(f'Number Loci: {aln_mat.num_loci}')
        logger.debug(f'Number Haplotypes: {aln_mat.num_haplotypes}')
        logger.debug(f'Number Reads: {aln_mat.num_reads}')
        aln_list.append(aln_mat)

    combined_aln = aln_list[0].copy()
    for i, aln in enumerate(aln_list[1:]):
        logger.info(f'Combining EMASE file: {emase_files[i]}')
        combined_aln = combined_aln.combine(aln)
        logger.debug(f'Combined Number Loci: {combined_aln.num_loci}')
        logger.debug(f'CombinedNumber Haplotypes: {combined_aln.num_haplotypes}')
        logger.debug(f'Combined Number Reads: {combined_aln.num_reads}')

    logger.info(f'Saving EMASE file {output_file}')
    combined_aln.save(h5file=output_file, complib=comp_lib)
    logger.info('Dome')


def count_alignments(
    alignment_file: str,
    group_file: str,
    outbase: str = 'emase'
):
    logger.info(f'Alignment File: {alignment_file}')
    logger.info(f'Group File: {group_file}')
    logger.info(f'Outbase: {outbase}')

    logger.info(f'Loading EMASE file: {alignment_file}')
    aln_mat = AlignmentPropertyMatrix(h5file=alignment_file, grpfile=group_file)
    logger.debug(f'Number Loci: {aln_mat.num_loci}')
    logger.debug(f'Number Haplotypes: {aln_mat.num_haplotypes}')
    logger.debug(f'Number Reads: {aln_mat.num_reads}')

    outfile1 = f'{outbase}.isoforms.alignment_counts'
    logger.info(f'Generating isoform Alignment Counts: {outfile1}')
    aln_mat.report_alignment_counts(filename=outfile1)

    aln_mat._bundle_inline(reset=True)

    outfile2 = f'{outbase}.genes.alignment_counts'
    logger.info(f'Generating gene Alignment Counts: {outfile2}')
    aln_mat.report_alignment_counts(filename=outfile2)

    logger.info('Done')


def get_num_shared_multireads(aln_mat: AlignmentPropertyMatrix) -> np.ndarray:
    hapsum = aln_mat.sum(axis=AlignmentPropertyMatrix.Axis.HAPLOTYPE)
    hapsum.data = np.ones(hapsum.nnz)
    cnt_mat = hapsum.transpose() * hapsum
    return cnt_mat


def count_shared_multireads_pairwise(
    alignment_file: str,
    group_file: str,
    outbase: str = 'emase'
):
    logger.info(f'Alignment File: {alignment_file}')
    logger.info(f'Group File: {group_file}')
    logger.info(f'Outbase: {outbase}')

    logger.info(f'Loading EMASE file: {alignment_file}')
    aln_mat = AlignmentPropertyMatrix(h5file=alignment_file, grpfile=group_file)
    logger.debug(f'Number Loci: {aln_mat.num_loci}')
    logger.debug(f'Number Haplotypes: {aln_mat.num_haplotypes}')
    logger.debug(f'Number Reads: {aln_mat.num_reads}')

    outfile1 = f'{outbase}.isoforms.shared_read_counts'
    logger.info(f'Generating isoform Shared Read Counts: {outfile1}')
    cnt_mat = get_num_shared_multireads(aln_mat)
    np.savez_compressed(outfile1, counts=cnt_mat)

    aln_mat._bundle_inline(reset=True)

    outfile2 = f'{outbase}.isoforms.shared_read_counts'
    logger.info(f'Generating genes Shared Read Counts: {outfile2}')
    cnt_mat = get_num_shared_multireads(aln_mat)
    np.savez_compressed(outfile2, counts=cnt_mat)

    logger.info('Done')


def create_hybrid(
    fasta_list: list[str],
    haplotypes: list[str],
    output_file: str = 'emase.pooled.targets.fa',
    build_bowtie_index: bool = False
):
    out_dir = os.path.dirname(output_file)
    if out_dir != '' and not os.path.exists(out_dir):
        os.mkdir(out_dir)

    for x in fasta_list:
        logger.info(f'Fasta File: {x}')
    logger.info(f'Haplotype List: {haplotypes}')
    logger.info(f'Output File: {output_file}')

    # Get pooled transcriptome
    outbase = os.path.splitext(output_file)[0]
    num_haps = len(fasta_list)
    lenfile = f'{outbase}.info'
    seqout = open(output_file, 'w')
    lenout = open(lenfile, 'w')

    logger.debug('Looping through haplotypes')
    for hid in range(num_haps):
        fasta = fasta_list[hid]
        hapname = haplotypes[hid]
        logger.info(
            f'Adding suffix "_{hapname}" to the sequence ID\'s of {fasta}'
        )
        fh = open(fasta)
        line = fh.readline()  # the first fasta header
        line = f'{line.rstrip().split()[0]}_{hapname}'
        seqout.write(f'{line}\n')
        lenout.write(f'{line[1:]}\t')
        seq_len = 0
        for line in fh:
            if line[0] == '>':
                line = f'{line.rstrip().split()[0]}_{hapname}\n'
                lenout.write(f'{seq_len}\n{line[1:].rstrip()}\t')
                seq_len = 0
            else:
                seq_len += len(line.rstrip())
            seqout.write(line)
        fh.close()
        lenout.write(f'{seq_len}\n')
    seqout.close()
    lenout.close()

    # Build bowtie index for the pooled transcriptome
    if build_bowtie_index:
        out_index = f'{outbase}.bowtie1'
        logger.info('Building bowtie1 index (could take some time)')
        status = subprocess.call(f'bowtie-build {output_file} {out_index}', shell=True)

    logger.info('Done')


def get_common_alignments(
    emase_files: list[str],
    output_file: str = None,
    comp_lib: str = 'zlib'
) -> None:
    if output_file is None:
        output_file = f'alignments.common.{os.path.basename(emase_files[0])}'

    for x in emase_files:
        logger.info(f'EMASE file: {x}')
    logger.info(f'Output File: {output_file}')
    logger.info(f'Compression Library: {comp_lib}')

    logger.info(f'Loading EMASE file: {emase_files[0]}')
    aln_mat = AlignmentPropertyMatrix(h5file=emase_files[0])
    logger.debug(f'Number Loci: {aln_mat.num_loci}')
    logger.debug(f'Number Haplotypes: {aln_mat.num_haplotypes}')
    logger.debug(f'Number Reads: {aln_mat.num_reads}')

    for f in emase_files[1:]:
        logger.info(f'Loading EMASE file: {f}')
        aln_mat_next = AlignmentPropertyMatrix(h5file=f)
        logger.debug(f'Number Loci: {aln_mat_next.num_loci}')
        logger.debug(f'Number Haplotypes: {aln_mat_next.num_haplotypes}')
        logger.debug(f'Number Reads: {aln_mat_next.num_reads}')

        if np.all(aln_mat.rname == aln_mat_next.rname):
            aln_mat = aln_mat * aln_mat_next
        else:
            logger.error('The read ID\'s are not compatible.')
            raise ValueError('The read ID\'s are not compatible.')

        logger.debug(f'Total Number Loci: {aln_mat.num_loci}')
        logger.debug(f'Total Number Haplotypes: {aln_mat.num_haplotypes}')
        logger.debug(f'Total Number Reads: {aln_mat.num_reads}')

    logger.info(f'Saving EMASE Formatted File: {output_file}')
    aln_mat.save(h5file=output_file, complib=comp_lib)
    logger.info('Done')


def pull_out_unique_reads(
    alignment_file: str,
    output_file: str,
    group_file: str = None,
    shallow: bool = False,
    ignore_alleles: bool = False
) -> None:
    logger.info(f'Alignment File: {alignment_file}')
    logger.info(f'Group File: {group_file}')
    logger.info(f'Output File: {output_file}')
    logger.info(f'Sahallow: {shallow}')
    logger.info(f'Ignore Alleles: {ignore_alleles}')

    logger.info(f'Loading EMASE file: {alignment_file}')
    aln_mat = AlignmentPropertyMatrix(h5file=alignment_file, grpfile=group_file)
    logger.debug(f'Number Loci: {aln_mat.num_loci}')
    logger.debug(f'Number Haplotypes: {aln_mat.num_haplotypes}')
    logger.debug(f'Number Reads: {aln_mat.num_reads}')

    logger.info('Getting unique reads')
    if group_file:
        logger.debug('Using group file')
        aln_mat_g = aln_mat.bundle(reset=True, shallow=shallow)
        aln_mat_g_uniq = aln_mat_g.get_unique_reads(
            ignore_haplotype=ignore_alleles, shallow=shallow
        )
        num_alns_per_read = aln_mat_g_uniq.sum(axis=AlignmentPropertyMatrix.Axis.LOCUS).sum(
            axis=AlignmentPropertyMatrix.Axis.HAPLOTYPE
        )
        aln_mat_uniq = aln_mat.pull_alignments_from(
            (num_alns_per_read > 0), shallow=shallow
        )
    else:
        logger.debug('Not using group file')
        aln_mat_uniq = aln_mat.get_unique_reads(
            ignore_haplotype=ignore_alleles, shallow=shallow
        )

    logger.info(f'Saving EMASE Formatted File: {output_file}')
    aln_mat_uniq.save(h5file=output_file, shallow=shallow)
    logger.info('Done')


def parse_gtf(gtf_fh):
    gdb = dict()
    tdb = dict()
    for line in dropwhile(utils.is_comment, gtf_fh):
        item = line.rstrip().split('\t')
        attribute = defaultdict(list)
        for e in item[8].split('; '):
            e1, e2 = e.replace(';', '').split('"')[:2]
            attribute[e1.strip()].append(e2.strip())
        attribute = dict(attribute)
        feature = item[2]
        s = int(item[3])
        e = int(item[4])
        # chro = item[0].replace('chr', '')
        chro = item[0]  # Leave the chromosome names as is in gtf file
        strand = item[6]
        if feature == 'gene':
            # Seqnature currently has some issue processing this entry
            gid = attribute.pop('gene_id')[0]
            if gid in gdb:
                print(f'[Error] Duplicate entry: {gid}')
            else:
                gdb[gid] = {
                    'chr': chro,
                    'strand': strand,
                    'start': s,
                    'end': e,
                    'isoform': set(),
                }
                for k, v in attribute.items():
                    if len(v) == 1:
                        v = v.pop()
                    gdb[gid][k] = v
        elif feature == 'transcript':
            gid = attribute['gene_id'][0]
            tid = attribute.pop('transcript_id')[0]
            if tid in tdb:
                print(f'[Error] Duplicate entry: {tid}')
            else:
                tdb[tid] = {
                    'chr': chro,
                    'strand': strand,
                    'start': s,
                    'end': e,
                    'eid': [],
                    'exon': [],
                    'exon_number': [],  # Do we need this?
                    'UTR': [],
                    'five_prime_utr': [],
                    'three_prime_utr': [],
                    'Selenocysteine': [],
                    'start_codon': [],
                    'stop_codon': [],
                    'start_codon_frame': [],
                    'stop_codon_frame': [],
                }
                for k, v in attribute.items():
                    if len(v) == 1:
                        v = v.pop()
                    tdb[tid][k] = v
                gid = tdb[tid]['gene_id']
                if gid in gdb:
                    gdb[gid]['isoform'].add(tid)
                else:  # This is a non-standard case where the input gtf does not have detailed gene-level annotation
                    gdb[gid] = dict()
                    gdb[gid]['chr'] = chro
                    gdb[gid]['isoform'] = set(tid)
        else:
            gid = attribute['gene_id'][0]
            tid = attribute['transcript_id'][0]
            if feature == 'exon':
                if tid in tdb:
                    try:
                        enu = int(attribute.pop('exon_number')[0])
                        tdb[tid]['exon_number'].append(enu)
                    except:
                        pass
                    try:
                        eid = attribute.pop('exon_id')[0]
                        tdb[tid]['eid'].append(eid)
                    except:
                        pass
                    tdb[tid]['exon'].append((s, e))
                else:  # This is a non-standard case where input gtf does not have detailed transcript-level annotation
                    tdb[tid] = dict()
                    tdb[tid]['chr'] = chro
                    tdb[tid]['strand'] = strand
                    tdb[tid]['exon'] = [(s, e)]
            elif feature in (
                'UTR',
                'five_prime_utr',
                'three_prime_utr',
                'Selenocysteine',
            ):
                tdb[tid][feature].append((s, e))
            elif feature in ('start_codon', 'stop_codon'):
                tdb[tid][feature].append((s, e))
                tdb[tid][feature + '_frame'].append(int(item[7]))
            elif feature == 'CDS':
                pass
            else:
                logger.error(f'[Unknown feature: {feature}]/n{line}')
    return gdb, tdb


#
# Get the regions of genes
#
def get_fragment(start, end, chro, strand, genome):
    fragment = genome[chro].seq[(start - 1): end]
    if strand == '-':
        fragment = fragment.reverse_complement()
    return fragment


def prepare(
    genome_files: list[str],
    haplotypes: list[str] = None,
    gtf_files: list[str] = None,
    out_dir: str = None,
    save_g2tmap: bool = False,
    save_dbs: bool = False,
    no_bowtie_index: bool = False
):
    """
    Prepare EMASE.

    Args:
        genome_files: list of genome files
        haplotypes: list of haplotypes
        gtf_files: list of gtf files
        out_dir: output directory
        save_g2tmap: whether to save g2tmap or not
        save_dbs: whether to save the db or not
        no_bowtie_index: True to not create index
    """
    for x in genome_files:
        logger.info(f'Genome File: {x}')
    logger.info(f'Haplotypes: {haplotypes}')
    for x in gtf_files:
        logger.info(f'GTF File: {x}')
    logger.info(f'Output Dir: {out_dir}')
    logger.info(f'Save Gene 2 Transcript Map: {save_g2tmap}')
    logger.info(f'Save DBs: {save_dbs}')
    logger.info(f'No bowtie Index: {no_bowtie_index}')

    num_haps = len(genome_files)

    if haplotypes is None:
        haplotypes = list(string.ascii_uppercase[:num_haps])
        if num_haps == 1:
            logger.info("Assuming single genome analysis. No suffix will be added to ID's")
        else:
            logger.info(f'Default haplotype names will be used: {", ".join(haplotypes)}')

    if gtf_files is None:
        gtf_files = []
        for genome_file in genome_files:
            gtf_files.append(f'{os.path.splitext(genome_file)[0]}.gtf')
        gtf_files_str = "\n".join(gtf_files)
        logger.info(f'Assuming there exist the following GTF files:\n{gtf_files_str}')

    if len(haplotypes) != num_haps or len(gtf_files) != num_haps:
        logger.warning('The number of gtf files or specified haplotypes is not matching to the number of genomes.')

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # Get pooled transcriptome
    if num_haps > 1:
        transcriptome_file = os.path.join(out_dir, 'emase.pooled.transcripts.fa')
        len_file = os.path.join(out_dir, 'emase.pooled.transcripts.info')
    else:
        transcriptome_file = os.path.join(out_dir, 'emase.transcripts.fa')
        len_file = os.path.join(out_dir, 'emase.transcripts.info')

    seq_out = open(transcriptome_file, 'w')
    len_out = open(len_file, 'w')

    for hid in range(num_haps):
        genome_file = genome_files[hid]
        genome_name = os.path.splitext(os.path.basename(genome_file))[0]
        hap_name = haplotypes[hid]
        gtf_file = gtf_files[hid]
        logger.info(f'Loading {genome_name} genome...')

        if os.path.splitext(genome_file)[1] == '.gz':
            genome_fh = gzip.open(genome_file, 'rb')
        else:
            genome_fh = open(genome_file)

        genome = SeqIO.to_dict(SeqIO.parse(genome_fh, 'fasta'))
        genome_fh.close()

        if os.path.splitext(gtf_file)[1] == '.gz':
            anno_fh = gzip.open(gtf_file, 'rb')
        else:
            anno_fh = open(gtf_file)

        logger.info(f'Parsing {os.path.basename(gtf_file)}...')
        gdb, tdb = parse_gtf(anno_fh)
        anno_fh.close()

        if num_haps == 1:
            print(
                f'Building {genome_name} transcriptome (Note: No suffix added to ID\'s)...'
            )
            if save_dbs:
                import _pickle as cPickle

                cPickle.dump(
                    gdb, open(os.path.join(out_dir, 'emase.gdb.pkl'), 'wb')
                )
                cPickle.dump(
                    tdb, open(os.path.join(out_dir, 'emase.tdb.pkl'), 'wb')
                )
        elif num_haps > 1:
            logger.info(
                f'Building {genome_name} transcriptome using suffix "_{hap_name}"...',
            )

        for tid in sorted(list(tdb.keys())):
            tinfo = tdb[tid]
            if tinfo['chr'] in genome:
                # Filter out transcripts from chromosome that the input genome does not contain
                if num_haps > 1:
                    # No need to add suffix if we deal with a single genome
                    tid = f'{tid}_{hap_name}'
                fragment = Seq('')
                for exon in tinfo['exon']:
                    fragment += get_fragment(
                        exon[0],
                        exon[1],
                        tinfo['chr'],
                        tinfo['strand'],
                        genome,
                    )
                if len(fragment) > 0:
                    SeqIO.write(
                        SeqRecord(fragment, tid, '', ''), seq_out, 'fasta'
                    )
                len_out.write(f'{tid}\t{len(fragment)}\n')
            else:
                print(
                    f'Skipping Transcript {tid} of Chromosome {tinfo["chr"]}...',
                )
    seq_out.close()
    len_out.close()

    if save_g2tmap:
        with open(
            os.path.join(out_dir, 'emase.gene2transcripts.tsv'), 'w'
        ) as fhout:
            logger.info("Recording mapping of gene id to transcript id's...",)
            for gid in sorted(list(gdb.keys())):
                if gdb[gid]['chr'] in genome:
                    item = [gid]
                    item = item + list(gdb[gid]['isoform'])
                    fhout.write('\t'.join(item) + '\n')

    #
    # Build bowtie index for the pooled transcriptome
    if not no_bowtie_index:
        out_index = os.path.join(
            os.path.dirname(transcriptome_file), 'bowtie.transcripts'
        )
        logger.info('Building bowtie index...')
        status = subprocess.call(
            f'bowtie-build {transcriptome_file} {out_index}', shell=True
        )

    logger.info('Done')


def run(
    alignment_file: str,
    group_file: str = None,
    length_file: str = None,
    outbase: str = 'emase',
    multiread_model: int = 4,
    read_length: int = 100,
    pseudocount: float = 0.0,
    max_iters: int = 999,
    tolerance: float = 0.0001,
    report_alignment_counts: bool = False,
    report_posterior: bool = False
):
    """
    Run EMASE.

    Args:
        alignment_file: The alignment file
        group_file: The group file
        length_file: The lengths file
        outbase: basename of all the generated output files
        multiread_model: which EMASE model to use
        read_length: the read length
        pseudocount: the prior read length
        max_iters: maximum iterations for EM iteration
        tolerance: tolerance for EM termination (default: 0.0001 in TPM)
        report_alignment_counts: whether to report alignment counts
        report_posterior: whether to report posterior probabilities
    """
    report_group_counts = (
        group_file is not None
    )

    logger.info(f'Alignment File: {alignment_file}')
    logger.info(f'Group File: {group_file}')
    logger.info(f'Read Length File: {length_file}')
    logger.info(f'Read Length: {read_length}')
    logger.info(f'Outbase: {outbase}')
    logger.info(f'Multiread Model: {multiread_model}')
    logger.info(f'Pseudocount: {pseudocount}')
    logger.info(f'Tolerance: {tolerance}')
    logger.info(f'Report Alignment Counts: {report_alignment_counts}')
    logger.info(f'Report Posterior: {report_posterior}')

    # load alignment incidence matrix ('alignment_file' is assumed to be in multiway transcriptome)
    logger.info(f'Loading EMASE file: {alignment_file}')
    aln_mat = AlignmentPropertyMatrix(h5file=alignment_file, grpfile=group_file)
    logger.debug(f'Number Loci: {aln_mat.num_loci}')
    logger.debug(f'Number Haplotypes: {aln_mat.num_haplotypes}')
    logger.debug(f'Number Reads: {aln_mat.num_reads}')

    logger.info('Running EMASE')
    em_factory = EMfactory(aln_mat)
    em_factory.prepare(pseudocount=pseudocount, lenfile=length_file, read_length=read_length)
    em_factory.run(
        model=multiread_model, tol=tolerance, max_iters=max_iters, verbose=True
    )

    logger.info(f'Generating isoform TPMs: {outbase}.isoforms.tpm')
    em_factory.report_depths(
        filename=f'{outbase}.isoforms.tpm', tpm=True
    )

    logger.info(f'Generating isoform Read Counts: {outbase}.isoforms.expected_read_counts')
    em_factory.report_read_counts(
        filename=f'{outbase}.isoforms.expected_read_counts'
    )

    if report_posterior:
        logger.info(f'Generating Posterior Probabilities: {outbase}.posterior.h5')
        em_factory.export_posterior_probability(
            filename=f'{outbase}.posterior.h5'
        )

    if report_group_counts:
        logger.info(f'Generating gene TPMs: {outbase}.genes.tpm')
        em_factory.report_depths(
            filename=f'{outbase}.genes.tpm',
            tpm=True,
            grp_wise=True
        )

        logger.info(f'Generating gene Read Counts: {outbase}.genes.expected_read_counts')
        em_factory.report_read_counts(
            filename=f'{outbase}.genes.expected_read_counts',
            grp_wise=True
        )

    if report_alignment_counts:
        aln_mat = AlignmentPropertyMatrix(h5file=alignment_file, grpfile=group_file)

        logger.info(f'Generating isoform Alignment Counts: {outbase}.isoforms.alignment_counts')
        aln_mat.report_alignment_counts(
            filename=f'{outbase}.isoforms.alignment_counts'
        )

        if report_group_counts:
            logger.info(f'Generating gene Alignment Counts: {outbase}.genes.alignment_counts')
            aln_mat._bundle_inline(reset=True)
            aln_mat.report_alignment_counts(
                filename=f'{outbase}.genes.alignment_counts'
            )
    logger.debug('Done')
