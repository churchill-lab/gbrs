# standard library imports
from pathlib import Path
from typing import Annotated, Optional
import importlib.metadata
import logging

# 3rd party library imports
import typer
from typer.core import TyperGroup
from rich.console import Console
from rich.panel import Panel

# local library imports
from gbrs import utils
from gbrs import h5_utils
from gbrs.emase import emase_utils

console = Console()

class SectionedGroup(TyperGroup):
    def format_help(self, ctx, formatter):
        console.print(Panel.fit(
            '[bold cyan]EMASE - Expression Modeling for Allele-Specific Expression[/bold cyan]'
        ))

        # Hardcoded groups
        main_commands = [
            'bam2emase', 'combine', 'count-alignments', 'count-shared-multireads-pairwise',
            'create-hybrid', 'get-common-alignments', 'get-common-alignments-optimized',
            'pull-out-unique-reads', 'prepare', 'run'
        ]
        utility_commands = ['h5-compare', 'h5-inspect']

        shown = set()

        console.print('\n[bold green]Main Commands:[/bold green]')
        for cmd in main_commands:
            if cmd in self.commands:
                help_text = self.commands[cmd].help or ''
                console.print(f'  • [cyan]{cmd}[/cyan] - {help_text}')
                shown.add(cmd)

        console.print('\n[bold yellow]Utility Commands:[/bold yellow]')
        for cmd in utility_commands:
            if cmd in self.commands:
                help_text = self.commands[cmd].help or ''
                console.print(f'  • [yellow]{cmd}[/yellow] - {help_text}')
                shown.add(cmd)

        # Show any commands not in your lists
        other_cmds = [cmd for cmd in self.commands if cmd not in shown]
        if other_cmds:
            console.print('\n[bold magenta]Other Commands:[/bold magenta]')
            for cmd in other_cmds:
                help_text = self.commands[cmd].help or ''
                console.print(f'  • [magenta]{cmd}[/magenta] - {help_text}')

        console.print('\n[dim]Use --help on each command for more details.[/dim]')


app = typer.Typer(cls=SectionedGroup, help='EMASE')


def version_callback(value: bool):
    if value:
        from gbrs import __logo_text__ as logo
        version = importlib.metadata.version('gbrs')
        typer.echo(f'{logo} {version}')
        raise typer.Exit()


@app.callback()
def common(
    ctx: typer.Context,
    version: bool = typer.Option(None, '--version', callback=version_callback),
):
    pass

@app.command(help='Convert BAM alignment files to EMASE format for allele-specific expression analysis')
def bam2emase(
    alignment_files: Annotated[list[Path], typer.Option('-i', '--alignment-files', exists=False, dir_okay=False, resolve_path=True, help='Input BAM file containing RNA-seq alignments, can separate files by "," or have multiple -i')],
    haplotypes: Annotated[list[str], typer.Option('-h', '--haplotype-char', help='Haplotype identifiers (e.g., A,B,C,D). Can specify multiple times or comma-separated')],
    locusid_file: Annotated[Path, typer.Option('-m', '--locus-ids', exists=True, dir_okay=False, resolve_path=True, help='Transcript/locus information file')],
    output_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True, help='Output EMASE file (HDF5 format). Auto-generated if not specified')] = None,
    delim: Annotated[str, typer.Option('-d', '--delim', help='Delimiter between transcript ID and haplotype in BAM file')] = '_',
    index_dtype: Annotated[str, typer.Option('--index-dtype', help='Data type for matrix indices (advanced users only)')] = 'uint32',
    data_dtype: Annotated[str, typer.Option('--data-dtype', help='Data type for matrix values (advanced users only)')] = 'uint8',
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('bam2emase_paired')
    try:
        all_alignment_files: list[str] = []
        for x in alignment_files:
            all_alignment_files.extend(str(x).split(','))

        for i, f in enumerate(all_alignment_files):
            all_alignment_files[i] = utils.check_file(f, 'r')

        all_haplotypes: list[str] = []
        for x in haplotypes:
            all_haplotypes.extend(x.split(','))

        locusid_file = str(locusid_file) if locusid_file else None
        output_file = str(output_file) if output_file else None

        emase_utils.bam2emase(
            alignment_files=all_alignment_files,
            haplotypes=all_haplotypes,
            locusid_file=locusid_file,
            output_file=output_file,
            delim=delim,
            index_dtype=index_dtype,
            data_dtype=data_dtype
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help='Merge multiple EMASE files into a single file')
def combine(
    emase_files: Annotated[list[Path], typer.Option('-i', '--emase-file', exists=False, dir_okay=False, resolve_path=True, help='EMASE files to merge. Can specify multiple files with comma separation or multiple -i flags')],
    output_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True, help='Output merged EMASE file')],
    comp_lib: Annotated[str, typer.Option('-c', '--comp-lib', help='Compression library for output file')] = 'zlib',
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('combine')
    try:
        # file shortcut: the following command line options are all equal
        # -i abc.h5 -i def.h5
        # -i abc.h5,def.h5
        all_emase_files: list[str] = []
        for x in emase_files:
            all_emase_files.extend(str(x).split(','))

        for i, f in enumerate(all_emase_files):
            all_emase_files[i] = utils.check_file(f, 'r')

        output_file = str(output_file) if output_file else None

        emase_utils.combine(
            emase_files=all_emase_files,
            output_file=output_file,
            comp_lib=comp_lib
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help='Generate alignment count statistics from EMASE files (isoform and gene levels)')
def count_alignments(
    alignment_file: Annotated[Path, typer.Option('-i', '--alignment-file', exists=True, dir_okay=False, resolve_path=True, help='Input EMASE file (HDF5 format)')],
    group_file: Annotated[Path, typer.Option('-g', '--group-file', exists=True, dir_okay=False, resolve_path=True, help='Gene-to-transcript mapping file (TSV format)')],
    outbase: Annotated[str, typer.Option('-o', '--outbase', help='Base name for output files (generates multiple files)')] = 'emase',
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('count_alignments')
    try:
        emase_utils.count_alignments(
            alignment_file=str(alignment_file),
            group_file=str(group_file),
            outbase=outbase
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help='Analyze pairwise shared multiread patterns between transcripts/genes')
def count_shared_multireads_pairwise(
    alignment_file: Annotated[Path, typer.Option('-i', '--alignment-file', exists=True, dir_okay=False, resolve_path=True, help='Input EMASE file (HDF5 format)')],
    group_file: Annotated[Path, typer.Option('-g', '--group-file', exists=True, dir_okay=False, resolve_path=True, help='Gene-to-transcript mapping file (TSV format)')],
    outbase: Annotated[str, typer.Option('-o', '--outbase', help='Base name for output files (generates multiple files)')] = 'emase',
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('count_shared_multireads_pairwise')
    try:
        emase_utils.count_shared_multireads_pairwise(
            alignment_file=str(alignment_file),
            group_file=str(group_file),
            outbase=outbase
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help='Create hybrid transcriptome from multiple haplotype-specific FASTA files')
def create_hybrid(
    fasta_list: Annotated[list[Path], typer.Option('-F', '--target-files', exists=False, dir_okay=False, resolve_path=True, help='Input FASTA files (one per haplotype). Can specify multiple files with comma separation or multiple -F flags')],
    haplotypes: Annotated[list[str], typer.Option('-s', '--suffices', help='Haplotype identifiers (e.g., A,B,C,D). Can specify multiple times or comma-separated')],
    output_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True)] = 'gbrs.hybridized.targets.fa',
    build_bowtie_index: bool = False,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('create_hybrid')
    try:
        # haplotype shortcut: the following command line options are all equal
        # -h A,B,C,D,E,F,G,H
        # -h A -h B -h C -h D -h E -h F -h G -h H
        # -h A,B,C,D -h E -h F -h G,H
        all_haplotypes: list[str] = []
        for x in haplotypes:
            all_haplotypes.extend(x.split(','))

        # file shortcut: the following command line options are all equal
        # -i abc.h5 -i def.h5
        # -i abc.h5,def.h5
        all_fasta_files: list[str] = []
        for x in fasta_list:
            all_fasta_files.extend(str(x).split(','))

        for i, f in enumerate(all_fasta_files):
            all_fasta_files[i] = utils.check_file(f, 'r')

        output_file = str(output_file) if output_file else None

        emase_utils.create_hybrid(
            fasta_list=all_fasta_files,
            haplotypes=all_haplotypes,
            output_file=str(output_file),
            build_bowtie_index=build_bowtie_index
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help='Find reads with identical alignment patterns across multiple EMASE files')
def get_common_alignments(
    emase_files: Annotated[list[Path], typer.Option('-i', '--emase-file', exists=False, dir_okay=False, resolve_path=True, help='EMASE files to compare. Can specify multiple files with comma separation or multiple -i flags')],
    output_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True, help='Output EMASE file containing only common alignments')] = None,
    comp_lib: Annotated[str, typer.Option('-c', '--comp-lib', help='Compression library for output file')] = 'zlib',
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('get_common_alignments')
    try:
        # file shortcut: the following command line options are all equal
        # -i abc.h5 -i def.h5
        # -i abc.h5,def.h5
        all_emase_files: list[str] = []
        for x in emase_files:
            all_emase_files.extend(str(x).split(','))

        for i, f in enumerate(all_emase_files):
            all_emase_files[i] = utils.check_file(f, 'r')

        output_file = str(output_file) if output_file else None

        emase_utils.get_common_alignments(
            emase_files=all_emase_files,
            output_file=output_file,
            comp_lib=comp_lib
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help='Optimized version of get-common-alignments with improved memory efficiency and validation options')
def get_common_alignments_optimized(
    emase_files: Annotated[list[Path], typer.Option('-i', '--emase-file', exists=False, dir_okay=False, resolve_path=True, help='EMASE files to compare. Can specify multiple files with comma separation or multiple -i flags')],
    output_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True, help='Output EMASE file containing only common alignments')] = None,
    comp_lib: Annotated[str, typer.Option('-c', '--comp-lib', help='Compression library for output file')] = 'zlib',
    skip_validate: Annotated[bool, typer.Option('-s', '--skip-validate', help='Skip validation of file compatibility (faster but less safe)')] = False,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('get_common_alignments_optimized')
    try:
        # file shortcut: the following command line options are all equal
        # -i abc.h5 -i def.h5
        # -i abc.h5,def.h5
        all_emase_files: list[str] = []
        for x in emase_files:
            all_emase_files.extend(str(x).split(','))

        for i, f in enumerate(all_emase_files):
            all_emase_files[i] = utils.check_file(f, 'r')

        output_file = str(output_file) if output_file else None

        emase_utils.get_common_alignments_optimized(
            emase_files=all_emase_files,
            output_file=output_file,
            comp_lib=comp_lib,
            validate=not skip_validate
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help='Compare two or more HDF5 files for differences in structure and data')
def h5_compare(
    h5_files: Annotated[list[Path], typer.Option('-i', '--h5-file', exists=False, dir_okay=False, resolve_path=True, help='HDF5 files to compare. Can specify multiple files with comma separation or multiple -i flags')],
    tolerant: Annotated[bool, typer.Option('-T', '--tolerant', help='Use tolerant comparison (allows small numerical differences)')] = False,
    tolerance: Annotated[float, typer.Option('-t', '--tolerance', help='Tolerance for floating point comparisons')] = 1e-6,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('h5_compare')
    try:
        # file shortcut: the following command line options are all equal
        # -i abc.h5 -i def.h5
        # -i abc.h5,def.h5
        all_h5_files: list[str] = []
        for x in h5_files:
            all_h5_files.extend(str(x).split(','))

        if len(h5_files) != 2:
            raise ValueError('Must provide exactly two h5 files to compare')

        for i, f in enumerate(all_h5_files):
            all_h5_files[i] = utils.check_file(f, 'r')

        h5_utils.compare_h5_files(all_h5_files[0], all_h5_files[1], tolerant, tolerance)
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help='Inspect bam2emase HDF5 file with debug information and matrix display')
def h5_inspect(
    h5_file: Annotated[Path, typer.Option('-i', '--h5-file', exists=False, dir_okay=False, resolve_path=True, help='HDF5 file to inspect')],
    haplotypes: Annotated[list[str], typer.Option('-h', '--haplotype', help='Haplotype names to display (e.g., A,B). Can specify multiple times or comma-separated. Default: all')] = None,
    show_matrix: Annotated[bool, typer.Option('-m', '--show-matrix', help='Show matrix for each haplotype')] = False,
    max_loci: Annotated[int, typer.Option('-l', '--max-loci', help='Maximum number of loci to display')] = 10,
    max_reads: Annotated[int, typer.Option('-r', '--max-reads', help='Maximum number of reads to display')] = 10,
    dense: Annotated[bool, typer.Option('-d', '--dense', help='Display matrices in dense format instead of sparse')] = False,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('h5_inspect')
    try:
        h5_file = utils.check_file(str(h5_file), 'r')
        
        # Process haplotype list
        all_haplotypes = None
        if haplotypes:
            all_haplotypes = []
            for x in haplotypes:
                all_haplotypes.extend(x.split(','))
        
        h5_utils.h5_inspect(h5_file, all_haplotypes, show_matrix, max_loci, max_reads, dense)
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help='Extract reads with unique alignment patterns from EMASE files')
def pull_out_unique_reads(
    alignment_file: Annotated[Path, typer.Option('-i', '--alignment-file', exists=True, dir_okay=False, resolve_path=True, help='Input EMASE file (HDF5 format)')],
    output_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True, help='Output EMASE file containing only unique reads')],
    group_file: Annotated[Path, typer.Option('-g', '--group-file', exists=True, dir_okay=False, resolve_path=True, help='Gene-to-transcript mapping file (TSV format)')] = None,
    shallow: Annotated[bool, typer.Option('-s', '--shallow', help='Create shallow EMASE file (faster but less complete)')] = False,
    ignore_alleles: Annotated[bool, typer.Option('-a', '--ignore-alleles', help='Ignore allele-level uniqueness requirements')] = False,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('pull_out_unique_reads')
    try:
        output_file = str(output_file) if output_file else None
        group_file = str(group_file) if group_file else None

        emase_utils.pull_out_unique_reads(
            alignment_file=str(alignment_file),
            group_file=group_file,
            output_file=output_file,
            shallow=shallow,
            ignore_alleles=ignore_alleles
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help='Prepare EMASE reference files from genome and annotation files')
def prepare(
    genome_file: Annotated[list[Path], typer.Option('-G', '--genome-file', exists=False, dir_okay=False, resolve_path=True, help='Genome FASTA files (one per haplotype). Can specify multiple files with comma separation or multiple -G flags')],
    haplotypes: Annotated[list[str], typer.Option('-s', '--haplotype-char', help='Haplotype identifiers (e.g., A,B,C,D). Can specify multiple times or comma-separated')] = None,
    gtf_file: Annotated[list[Path], typer.Option('-g', '--gtf-file', exists=False, dir_okay=False, resolve_path=True, help='Gene annotation files (GTF format). Can specify multiple files with comma separation or multiple -g flags')] = None,
    out_dir: Annotated[str, typer.Option('-o', '--out_dir', help='Output directory for reference files (default: current directory)')] = None,
    save_g2tmap: Annotated[bool, typer.Option('-m', '--save-g2tmap', help='Save gene-to-transcript mapping file')] = False,
    save_dbs: Annotated[bool, typer.Option('-d', '--save-dbs', help='Save database files for debugging')] = False,
    no_bowtie_index: Annotated[bool, typer.Option('-m', '--no-bowtie-index', help='Skip building Bowtie index (saves time)')] = False,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('run')
    try:
        # haplotype shortcut: the following command line options are all equal
        # -h A,B,C,D,E,F,G,H
        # -h A -h B -h C -h D -h E -h F -h G -h H
        # -h A,B,C,D -h E -h F -h G,H
        all_haplotypes: list[str] = []
        for x in haplotypes:
            all_haplotypes.extend(x.split(','))

        # file shortcut: the following command line options are all equal
        # -i abc.h5 -i def.h5
        # -i abc.h5,def.h5
        all_genome_files: list[str] = []
        for x in genome_file:
            all_genome_files.extend(str(x).split(','))

        for i, f in enumerate(all_genome_files):
            all_genome_files[i] = utils.check_file(f, 'r')

        all_gtf_files: list[str] = []
        for x in gtf_file:
            all_gtf_files.extend(str(x).split(','))

        for i, f in enumerate(all_gtf_files):
            all_gtf_files[i] = utils.check_file(f, 'r')

        emase_utils.prepare(
            genome_files=all_genome_files,
            haplotypes=all_haplotypes,
            gtf_files=all_gtf_files,
            out_dir=out_dir,
            save_g2tmap=save_g2tmap,
            save_dbs=save_dbs,
            no_bowtie_index=no_bowtie_index
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help='Run EMASE algorithm to quantify allele-specific expression from alignment data')
def run(
    alignment_file: Annotated[Path, typer.Option('-i', '--alignment-file', exists=True, dir_okay=False, resolve_path=True, help='Input EMASE file (HDF5 format)')],
    group_file: Annotated[Path, typer.Option('-g', '--group-file', exists=True, dir_okay=False, resolve_path=True, help='Gene-to-transcript mapping file (TSV format)')] = None,
    length_file: Annotated[Path, typer.Option('-L', '--length-file', exists=True, dir_okay=False, resolve_path=True, help='Transcript length file (TSV format)')] = None,
    outbase: Annotated[str, typer.Option('-o', '--outbase', help='Base name for output files (generates multiple files)')] = 'emase',
    multiread_model: Annotated[int, typer.Option('-M', '--multiread-model', help='EMASE multiread model (1-4, default: 4)')] = 4,
    pseudocount: Annotated[float, typer.Option('-p', '--pseudocount', help='Prior read count for allele specificity')] = 0.0,
    read_length: Annotated[int, typer.Option('-l', '--read-length', help='Read length for length normalization')] = 100,
    max_iters: Annotated[int, typer.Option('-m', '--max-iters', help='Maximum EM iterations')] = 999,
    tolerance: Annotated[float, typer.Option('-t', '--tolerance', help='Convergence tolerance in TPM')] = 0.0001,
    report_alignment_counts: Annotated[bool, typer.Option('-c', '--report-alignment-counts', help='Generate alignment count files')] = False,
    report_posterior: Annotated[bool, typer.Option('-w', '--report-posterior', help='Generate posterior probability files')] = False,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('run')
    try:
        if not multiread_model in (1, 2, 3, 4):
            raise typer.Abort(f'-M, --multiread-model must be one of 1, 2, 3, or 4')

        group_file = str(group_file) if group_file else None
        length_file = str(length_file) if length_file else None

        emase_utils.run(
            alignment_file=str(alignment_file),
            group_file=group_file,
            length_file=length_file,
            outbase=outbase,
            multiread_model=multiread_model,
            read_length=read_length,
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


if __name__ == '__main__':
    app()
