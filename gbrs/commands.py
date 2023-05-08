# standard library imports
from pathlib import Path
from typing import Annotated
import logging

# 3rd party library imports
from rich.logging import RichHandler
import typer

# local library imports
from gbrs import emase_utils
from gbrs import gbrs_utils


ch1 = RichHandler(level=logging.NOTSET, show_level=True, show_time=True, show_path=False, omit_repeated_times=False)
ch1.addFilter(gbrs_utils.NoDebugLogFilter())
ch2 = RichHandler(level=logging.NOTSET, show_level=True, show_time=True, show_path=True, omit_repeated_times=False)
ch2.addFilter(gbrs_utils.DebugLogFilter())

logging.basicConfig(
    level="NOTSET",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[ch1, ch2]
)

app = typer.Typer()


@app.command()
def hybridize(
    fasta_list: Annotated[list[Path], typer.Option('-F', '--target-files', exists=True, dir_okay=False, resolve_path=True)],
    hap_list: Annotated[str, typer.Option('-s', '--suffices')],
    out_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True)] = 'gbrs.hybridized.targets.fa',
    build_bowtie_index: bool = False,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help="specify multiple times for more verbose output")] = 0
) -> None:
    logger = gbrs_utils.configure_logging(verbose)
    logger.debug('hybridize')
    try:
        emase_utils.hybridize(
            fasta_list=fasta_list,
            hap_list=hap_list,
            out_file=str(out_file),
            build_bowtie_index=build_bowtie_index
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help="convert a BAM file to EMASE format")
def bam2emase(
    alignment_file: Annotated[Path, typer.Option('-i', '--alignment-file', exists=True, dir_okay=False, resolve_path=True, help="bam file to convert")],
    haplotypes: Annotated[list[str], typer.Option('-h', '--haplotype-char', help='haplotype, one per -h option, i.e. -h A -h B -h C, also a shortcut can be -h A,B,C')],
    locusid_file: Annotated[Path, typer.Option('-m', '--locus-ids', exists=True, dir_okay=False, resolve_path=True, help='filename for the locus (usually transcripts) info')] = None,
    out_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True, help="EMASE file (hdf5 format)")] = None,
    delim: Annotated[str, typer.Option('-d', '--delim', help='delimiter string between locus and haplotype in BAM file')] = '_',
    index_dtype: Annotated[str, typer.Option('--index-dtype', help='advanced option, see internal code')] = 'uint32',
    data_dtype: Annotated[str, typer.Option('--data-dtype', help='advanced_option, see internal code')] = 'uint8',
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help="specify multiple times for more verbose output")] = 0
) -> None:
    logger = gbrs_utils.configure_logging(verbose)
    logger.debug('bam2emase')
    try:
        # haplotype shortcut: the following command line options are all equal
        # -h A,B,C,D,E,F,G,H
        # -h A -h B -h C -h D -h E -h F -h G -h H
        # -h A,B,C,D -h E -h F -h G,H
        all_haplotypes: list[str] = []
        for x in haplotypes:
            all_haplotypes.extend(x.split(','))

        emase_utils.bam2emase(
            alignment_file=alignment_file,
            haplotypes=all_haplotypes,
            locusid_file=locusid_file,
            out_file=out_file,
            delim=delim,
            index_dtype=index_dtype,
            data_dtype=data_dtype
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help="compress EMASE format alignment incidence matrix")
def compress(
    emase_files: Annotated[list[Path], typer.Option('-i', '--emase-file', exists=True, dir_okay=False, resolve_path=True, help='EMASE file to compress, can seperate files by "," or have multiple -i')],
    out_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True, help='name of the compressed EMASE file')],
    comp_lib: Annotated[str, typer.Option('-c', '--comp-lib', help='compression library to use')] = 'zlib',
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help="specify multiple times for more verbose output")] = 0
) -> None:
    logger = gbrs_utils.configure_logging(verbose)
    logger.debug('compress')
    try:
        # file shortcut: the following command line options are all equal
        # -i abc.h5 -i def.h5
        # -i abc.h5,def.h5
        all_emase_files: list[str] = []
        for x in emase_files:
            all_emase_files.extend(str(x).split(','))

        emase_utils.compress(
            emase_files=all_emase_files,
            out_file=out_file,
            comp_lib=comp_lib
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help="quantify allele-specific expressions")
def quantify(
    alignment_file: Annotated[Path, typer.Option('-i', '--alignment-file', exists=True, dir_okay=False, resolve_path=True, help='EMASE alignment incidence file (in hdf5 format)')],
    group_file: Annotated[Path, typer.Option('-g', '--group-file', exists=True, dir_okay=False, resolve_path=True, help='tab delimited file of gene to transcript mapping')] = None,
    length_file: Annotated[Path, typer.Option('-L', '--length-file', exists=True, dir_okay=False, resolve_path=True, help='tab delimited file of locus(transcript) and length')] = None,
    genotype_file: Annotated[Path, typer.Option('-G', '--genotype', exists=True, dir_okay=False, resolve_path=True, help='tab delimited file of locus(transcipt) and diplotype')] = None,
    outbase: Annotated[str, typer.Option('-o', '--outbase', help='basename of all the generated output files')] = 'gbrs.quantified',
    multiread_model: Annotated[int, typer.Option('-M', '--multiread-model', help='emase model (default: 4)')] = 4,
    pseudocount: Annotated[float, typer.Option('-p', '--pseudocount', help='prior read count (default: 0.0)')] = 0.0,
    max_iters: Annotated[int, typer.Option('-m', '--max-iters', help='maximum iterations for EM iteration')] = 999,
    tolerance: Annotated[float, typer.Option('-t', '--tolerance', help='tolerance for EM termination (default: 0.0001 in TPM)')] = 0.0001,
    report_alignment_counts: Annotated[bool, typer.Option('-a', '--report-alignment-counts', help='whether to report alignment counts')] = False,
    report_posterior: Annotated[bool, typer.Option('-w', '--report-alignment-counts', help='whether to report posterior probabilities')] = False,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help="specify multiple times for more verbose output")] = 0
) -> None:
    logger = gbrs_utils.configure_logging(verbose)
    logger.debug('quantify')
    try:
        if not multiread_model in (1, 2, 3, 4):
            raise typer.Abort(f'-M, --multiread-model must be one of 1, 2, 3, or 4')

        emase_utils.quantify(
            alignment_file=alignment_file,
            group_file=group_file,
            length_file=length_file,
            genotype_file=genotype_file,
            outbase=outbase,
            multiread_model=multiread_model,
            pseudocount=pseudocount,
            max_iters=max_iters,
            tolerance=tolerance,
            report_alignment_counts=report_alignment_counts,
            report_posterior=report_posterior
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help="reconstruct the genome based upon gene-level TPM quantities")
def reconstruct(
    expression_file: Annotated[Path, typer.Option('-e', '--expr-file', exists=True, dir_okay=False, resolve_path=True, help='file containing gene-level TPM quantities')],
    tprob_file: Annotated[Path, typer.Option('-t', '--tprob-file', exists=True, dir_okay=False, resolve_path=True, help='transition probabilities file')],
    avec_file: Annotated[Path, typer.Option('-x', '--avec-file', exists=True, dir_okay=False, resolve_path=True)] = None,
    gpos_file: Annotated[Path, typer.Option('-g', '--gpos-file', exists=True, dir_okay=False, resolve_path=True)] = None,
    expr_threshold: Annotated[float, typer.Option('-c', '--expr-threshold')] = 1.5,
    sigma: Annotated[float, typer.Option('-s', '--sigma')] = 0.12,
    outbase: Annotated[str, typer.Option('-o', '--outbase', help='basename of all the generated output files')] = None,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help="specify multiple times for more verbose output")] = 0
) -> None:
    logger = gbrs_utils.configure_logging(verbose)
    logger.debug('reconstruct')
    try:
        gbrs_utils.reconstruct(
            expression_file=expression_file,
            tprob_file=tprob_file,
            avec_file=avec_file,
            gpos_file=gpos_file,
            expr_threshold=expr_threshold,
            sigma=sigma,
            outbase=outbase
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help="interpolate probability on a decently-spaced grid")
def interpolate(
    probability_file: Annotated[Path, typer.Option('-i', '--probability-file', exists=True, dir_okay=False, resolve_path=True)],
    grid_file: Annotated[Path, typer.Option('-g', '--grid-file', exists=True, dir_okay=False, resolve_path=True)] = None,
    gpos_file: Annotated[Path, typer.Option('-p', '--gpos-file', exists=True, dir_okay=False, resolve_path=True)] = None,
    out_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True)] = None,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help="specify multiple times for more verbose output")] = 0
) -> None:
    logger = gbrs_utils.configure_logging(verbose)
    logger.debug('interpolate')
    try:
        gbrs_utils.interpolate(
            probability_file=probability_file,
            grid_file=grid_file,
            gpos_file=gpos_file,
            out_file=out_file
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help="export to GBRS quant format")
def export(
    genoprob_file: Annotated[Path, typer.Option('-i', '--genoprob-file', exists=True, dir_okay=False, resolve_path=True)],
    strains: Annotated[str, typer.Option('-s', '--strains')],
    grid_file: Annotated[Path, typer.Option('-g', '--grid-file', exists=True, dir_okay=False, resolve_path=True)] = None,
    out_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True)] = None,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help="specify multiple times for more verbose output")] = 0
) -> None:
    logger = gbrs_utils.configure_logging(verbose)
    logger.debug('export')
    try:
        gbrs_utils.export(
            genoprob_file=genoprob_file,
            strains=strains,
            grid_file=grid_file,
            out_file=out_file
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help="plot a reconstructed genome")
def plot(
    genoprob_file: Annotated[Path, typer.Option('-i', '--genoprob-file', exists=True, dir_okay=False, resolve_path=True)],
    out_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True)] = None,
    sample_name: Annotated[str, typer.Option('-n', '--sample_name')] = None,
    grid_size: Annotated[int, typer.Option('-g', '--num-grids')] = 2,
    xt_max: Annotated[int, typer.Option('-m', '--xt-max')] = 5000,
    xt_size: Annotated[int, typer.Option('-s', '--xt_size')] = 500,
    grid_width: Annotated[float, typer.Option('-w', '--grid-width')] = 0.01,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help="specify multiple times for more verbose output")] = 0
):
    logger = gbrs_utils.configure_logging(verbose)
    logger.debug('plot')
    try:
        gbrs_utils.plot(
            genoprob_file=genoprob_file,
            out_file=out_file,
            sample_name=sample_name,
            grid_size=grid_size,
            xt_max=xt_max,
            xt_size=xt_size,
            grid_width=grid_width
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command()
def get_transition_prob(
    marker_file: Annotated[Path, typer.Option('-i', '--marker-file', exists=True, dir_okay=False, resolve_path=True)],
    haplotypes: Annotated[str, typer.Option('-s', '--haplotypes')] = 'A,B',
    mating_scheme: Annotated[str, typer.Option('-m', '--mating-scheme')] = 'RI',
    gamma_scale: Annotated[float, typer.Option('-g', '--gamma-scale')] = 0.01,
    epsilon: Annotated[float, typer.Option('-e', '--epsilon')] = 0.000001,
    out_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True)] = 'tranprob.npz',
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help="specify multiple times for more verbose output")] = 0
) -> None:
    logger = gbrs_utils.configure_logging(verbose)
    logger.debug('get-transition-prob')
    try:
        gbrs_utils.get_transition_prob(
            marker_file=marker_file,
            haplotypes=haplotypes,
            mating_scheme=mating_scheme,
            gamma_scale=gamma_scale,
            epsilon=epsilon,
            out_file=out_file
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command()
def get_alignment_spec(
    sample_file: Annotated[Path, typer.Option('-i', '--sample-file', exists=True, dir_okay=False, resolve_path=True)],
    haplotypes: Annotated[str, typer.Option('-s', '--parental-strains')],
    min_expr: Annotated[float, typer.Option('-m', '--min-expr')] = 2.0,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help="specify multiple times for more verbose output")] = 0
) -> None:
    logger = gbrs_utils.configure_logging(verbose)
    logger.debug('get_alignment_spec')
    try:
        gbrs_utils.get_alignment_spec(
            sample_file=sample_file,
            haplotypes=haplotypes,
            min_expr=min_expr
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help="apply genotype calls to multi-way alignment incidence matrix")
def stencil(
    alignment_file: Annotated[Path, typer.Option('-i', '--alignment-file', exists=True, dir_okay=False, resolve_path=True, help="alignment incidence file (h5)")],
    genotype_file: Annotated[Path, typer.Option('-G', '--genotype', exists=True, dir_okay=False, resolve_path=True, help="genotype calls by GBRS (tsv)")],
    group_file: Annotated[Path, typer.Option('-g', '--group-file', exists=True, dir_okay=False, resolve_path=True, help="gene ID to isoform ID mapping info (tsv)")] = None,
    out_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True, help="enotyped version of alignment incidence file (h5)")] = None,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help="specify multiple times for more verbose output")] = 0
) -> None:
    logger = gbrs_utils.configure_logging(verbose)
    logger.debug('stencil')
    try:
        emase_utils.stencil(
            alignment_file=alignment_file,
            genotype_file=genotype_file,
            group_file=group_file,
            out_file=out_file
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command()
def intersect(
    emase_files: Annotated[list[Path], typer.Option('-i', '--emase-file', exists=True, dir_okay=False, resolve_path=True)],
    out_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True)] = None,
    comp_lib: Annotated[str, typer.Option('-c', '--comp-lib')] = 'zlib',
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help="specify multiple times for more verbose output")] = 0
) -> None:
    logger = gbrs_utils.configure_logging(verbose)
    logger.debug('intersect')
    try:
        # file shortcut: the following command line options are all equal
        # -i abc.h5 -i def.h5
        # -i abc.h5,def.h5
        all_emase_files: list[str] = []
        for x in emase_files:
            all_emase_files.extend(str(x).split(','))

        emase_utils.intersect(
            emase_files=all_emase_files,
            out_file=out_file,
            comp_lib=comp_lib
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


if __name__ == '__main__':
    app()
