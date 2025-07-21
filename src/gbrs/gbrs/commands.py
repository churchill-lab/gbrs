# standard library imports
import importlib.metadata
import logging
from pathlib import Path
from typing import Annotated, Optional

# 3rd party library imports
import typer
from typer.core import TyperGroup
from rich.console import Console
from rich.panel import Panel

# local library imports
from gbrs.gbrs import bam_utils
from gbrs.gbrs import emase_utils
from gbrs.gbrs import gbrs_utils
from gbrs import utils

console = Console()

class SectionedGroup(TyperGroup):
    def format_help(self, ctx, formatter):
        console.print(Panel.fit('[bold cyan]GBRS - Genome Reconstruction and Quantification Suite[/bold cyan]'))
        
        # Hardcoded groups
        main_commands = [
            'compress', 'compress', 'quantify', 'reconstruct', 'interpolate',
            'export', 'plot',  'get-transition-prob', 'get-alignment-spec', 'stencil'
        ]
        utility_commands = ['debug-genoprob', 'generate_bam']

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


app = typer.Typer(cls=SectionedGroup, help='GBRS')


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


@app.command(help='Compress EMASE files to reduce storage size while maintaining data integrity')
def compress(
    emase_files: Annotated[list[Path], typer.Option('-i', '--emase-file', exists=False, dir_okay=False, resolve_path=True, help='EMASE files to compress. Can specify multiple files with comma separation or multiple -i flags')],
    output_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True, help='Output compressed EMASE file')],
    comp_lib: Annotated[str, typer.Option('-c', '--comp-lib', help='Compression library for output file')] = 'zlib',
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('compress')
    try:
        # file shortcut: the following command line options are all equal
        # -i abc.h5 -i def.h5
        # -i abc.h5,def.h5
        all_emase_files: list[str] = []
        for x in emase_files:
            all_emase_files.extend(str(x).split(','))

        for i, f in enumerate(all_emase_files):
            all_emase_files[i] = utils.check_file(f, 'r')

        emase_utils.compress(
            emase_files=all_emase_files,
            output_file=str(output_file),
            comp_lib=comp_lib
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help='Quantify allele-specific expression from EMASE alignment data using the GBRS model')
def quantify(
    alignment_file: Annotated[Path, typer.Option('-i', '--alignment-file', exists=True, dir_okay=False, resolve_path=True, help='Input EMASE alignment file (HDF5 format)')],
    group_file: Annotated[Path, typer.Option('-g', '--group-file', exists=True, dir_okay=False, resolve_path=True, help='Gene-to-transcript mapping file (TSV format)')] = None,
    length_file: Annotated[Path, typer.Option('-L', '--length-file', exists=True, dir_okay=False, resolve_path=True, help='Transcript length file (TSV format)')] = None,
    genotype_file: Annotated[Path, typer.Option('-G', '--genotype', exists=True, dir_okay=False, resolve_path=True, help='Genotype calls file (TSV format)')] = None,
    outbase: Annotated[str, typer.Option('-o', '--outbase', help='Base name for output files (generates multiple files)')] = 'gbrs.quantified',
    multiread_model: Annotated[int, typer.Option('-M', '--multiread-model', help='GBRS multiread model (1-4, default: 4)')] = 4,
    pseudocount: Annotated[float, typer.Option('-p', '--pseudocount', help='Prior read count for allele specificity')] = 0.0,
    max_iters: Annotated[int, typer.Option('-m', '--max-iters', help='Maximum EM iterations')] = 999,
    tolerance: Annotated[float, typer.Option('-t', '--tolerance', help='Convergence tolerance in TPM')] = 0.0001,
    report_alignment_counts: Annotated[bool, typer.Option('-a', '--report-alignment-counts', help='Generate alignment count files')] = False,
    report_posterior: Annotated[bool, typer.Option('-w', '--report-posterior', help='Generate posterior probability files')] = False,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('quantify')
    try:
        if multiread_model not in (1, 2, 3, 4):
            raise typer.Abort('-M, --multiread-model must be one of 1, 2, 3, or 4')

        group_file = str(group_file) if group_file else None
        length_file = str(length_file) if length_file else None
        genotype_file = str(genotype_file) if genotype_file else None

        emase_utils.quantify(
            alignment_file=str(alignment_file),
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


@app.command(help='Reconstruct the genome based on gene-level TPM quantities using the GBRS model')
def reconstruct(
    expression_file: Annotated[Path, typer.Option('-e', '--expr-file', exists=True, dir_okay=False, resolve_path=True, help='File containing gene-level TPM quantities')],
    tprob_file: Annotated[Path, typer.Option('-t', '--tprob-file', exists=True, dir_okay=False, resolve_path=True, help='Transition probabilities file (NPZ format)')],
    avec_file: Annotated[Path, typer.Option('-x', '--avec-file', exists=True, dir_okay=False, resolve_path=True, help='Alignment specificity file (optional)')] = None,
    gpos_file: Annotated[Path, typer.Option('-g', '--gpos-file', exists=True, dir_okay=False, resolve_path=True, help='Gene meta information file (chromosome, id, location)')] = None,
    expr_threshold: Annotated[float, typer.Option('-c', '--expr-threshold', help='Minimum TPM threshold for expression')] = 1.5,
    sigma: Annotated[float, typer.Option('-s', '--sigma', help='Standard deviation for HMM emission probabilities')] = 0.12,
    outbase: Annotated[str, typer.Option('-o', '--outbase', help='Base name for output files (generates multiple files)')] = None,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('reconstruct')
    try:
        avec_file = str(avec_file) if avec_file else None
        gpos_file = str(gpos_file) if gpos_file else None

        gbrs_utils.reconstruct(
            expression_file=str(expression_file),
            tprob_file=str(tprob_file),
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


@app.command(help='Interpolate genotype probabilities onto a regular genomic grid for downstream analysis')
def interpolate(
    genoprob_file: Annotated[Path, typer.Option('-i', '--genoprob-file', exists=True, dir_okay=False, resolve_path=True, help='Input genotype probability file (NPZ or TSV format)')],
    grid_file: Annotated[Path, typer.Option('-g', '--grid-file', exists=True, dir_okay=False, resolve_path=True, help='Grid file specifying genomic positions (e.g., ref.genome_grid.64k.txt)')] = None,
    gpos_file: Annotated[Path, typer.Option('-p', '--gpos-file', exists=True, dir_okay=False, resolve_path=True, help='Gene meta information file (chromosome, id, location)')] = None,
    output_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True, help='Output file in GBRS quant format')] = None,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('interpolate')
    try:
        grid_file = str(grid_file) if grid_file else None
        gpos_file = str(gpos_file) if gpos_file else None
        output_file = str(output_file) if output_file else None

        gbrs_utils.interpolate(
            genoprob_file=str(genoprob_file),
            grid_file=grid_file,
            gpos_file=gpos_file,
            output_file=output_file
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help='Export genotype probabilities to GBRS quant format for downstream analysis')
def export(
    genoprob_file: Annotated[Path, typer.Option('-i', '--genoprob-file', exists=True, dir_okay=False, resolve_path=True, help='Input genotype probability file (NPZ or TSV format)')],
    strains: Annotated[list[str], typer.Option('-s', '--strains', help='Strain identifiers (e.g., A,B,C,D). Can specify multiple times or comma-separated')],
    grid_file: Annotated[Path, typer.Option('-g', '--grid-file', exists=True, dir_okay=False, resolve_path=True, help='Grid file specifying genomic positions (e.g., ref.genome_grid.64k.txt)')] = None,
    output_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True, help='Output file in GBRS quant format')] = None,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('export')
    try:
        # strains shortcut: the following command line options are all equal
        # -s A,B,C,D,E,F,G,H
        # -s A,B,C,D -s E -s F -s G,H
        all_strains: list[str] = []
        for x in strains:
            all_strains.extend(x.split(','))

        grid_file = str(grid_file) if grid_file else None
        output_file = str(output_file) if output_file else None

        gbrs_utils.export(
            genoprob_file=str(genoprob_file),
            strains=all_strains,
            grid_file=grid_file,
            output_file=output_file
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help='Plot a reconstructed genome from genotype probabilities')
def plot(
    genoprob_file: Annotated[Path, typer.Option('-i', '--genoprob-file', exists=True, dir_okay=False, resolve_path=True, help='Input genotype probabilities (NPZ format)')],
    output_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True, help='Output plot file (e.g., PDF)')] = None,
    output_format: Annotated[str, typer.Option('-f', '--format', help='Output file format (e.g., pdf, png)')] = 'pdf',
    sample_name: Annotated[str, typer.Option('-n', '--sample_name', help='Sample name for labeling the plot')] = None,
    grid_size: Annotated[int, typer.Option('-g', '--num-grids', hidden=True)] = 2,
    xt_max: Annotated[int, typer.Option('-m', '--xt-max', hidden=True)] = 5000,
    xt_size: Annotated[int, typer.Option('-s', '--xt_size', hidden=True)] = 500,
    grid_width: Annotated[float, typer.Option('-w', '--grid-width', hidden=True)] = 0.01,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('plot')
    try:
        output_file = str(output_file) if output_file else None

        gbrs_utils.plot(
            genoprob_file=str(genoprob_file),
            output_file=output_file,
            output_format=output_format,
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

@app.command(help='Plot genotype probabilities across the genome with customizable visualization options')
def plot_alt(
    genoprobs_file: Annotated[Path, typer.Option('-i', '--genoprobs-file', exists=True, dir_okay=False, resolve_path=True, help='Input genotype probabilities file (NPZ format)')],
    output_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True, help='Output plot file (e.g., PDF)')] = None,
    output_format: Annotated[str, typer.Option('-f', '--format', help='Output file format (e.g., pdf, png)')] = 'pdf',
    sample_name: Annotated[str, typer.Option('-n', '--sample_name', help='Sample name for labeling the plot')] = None,
    founder_colors_file: Annotated[Path, typer.Option('-C', '--colors', exists=True, dir_okay=False, resolve_path=True, help='Tab-delimited file of founder and hex color')] = None,
    genome_pos_file: Annotated[Path, typer.Option('-P', '--positions', exists=True, dir_okay=False, resolve_path=True, help='Gene position .npz file for cM positions')] = None,
    dpi: Annotated[int, typer.Option('-D', '--dpi', help='DPI for output')] = 300,
    bar_height: Annotated[float, typer.Option('--bar-height', help='Height of each haplotype bar')] = 0.5,
    haplotype_gap: Annotated[float, typer.Option('--founder-gap', help='Gap between haplotype bars')] = 0.1,
    chrom_spacing: Annotated[float, typer.Option('--chrom-spacing', help='Vertical spacing between chromosomes')] = 2.0,
    fig_width: Annotated[int, typer.Option('--fig-width', help='Figure width in inches')] = 18,
    min_fig_height: Annotated[int, typer.Option('--min-fig-height', help='Minimum figure height in inches')] = 8,
    height_per_chrom: Annotated[float, typer.Option('--height-per-chrom', help='Height per chromosome in inches')] = 0.8,
    font_size_title: Annotated[int, typer.Option('--font-size-title', help='Title font size')] = 22,
    font_size_axis: Annotated[int, typer.Option('--font-size-axis', help='Axis label font size')] = 18,
    font_size_tick: Annotated[int, typer.Option('--font-size-tick', help='Tick label font size')] = 16,
    font_size_legend: Annotated[int, typer.Option('--font-size-legend', help='Legend font size')] = 14,
    legend_line_width: Annotated[int, typer.Option('--legend-line-width', help='Legend line width')] = 8,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('plot_alt')
    try:
        output_file = str(output_file) if output_file else None
        founder_colors_file = str(founder_colors_file) if founder_colors_file else None
        genome_pos_file = str(genome_pos_file) if genome_pos_file else None

        gbrs_utils.plot_genoprobs(
            genoprobs_file=str(genoprobs_file),
            output_file=output_file,
            output_format=output_format,
            sample_name=sample_name,
            founder_colors_file=founder_colors_file,
            genome_pos_file=genome_pos_file,
            dpi=dpi,
            bar_height=bar_height,
            haplotype_gap=haplotype_gap,
            chrom_spacing=chrom_spacing,
            fig_width=fig_width,
            min_fig_height=min_fig_height,
            height_per_chrom=height_per_chrom,
            font_size_title=font_size_title,
            font_size_axis=font_size_axis,
            font_size_tick=font_size_tick,
            font_size_legend=font_size_legend,
            legend_line_width=legend_line_width
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)

@app.command(help='Calculate transition probabilities for a given marker set and mating scheme')
def get_transition_prob(
    marker_file: Annotated[Path, typer.Option('-i', '--marker-file', exists=True, dir_okay=False, resolve_path=True, help='Input marker file (TSV Format)')],
    haplotypes: Annotated[str, typer.Option('-s', '--haplotypes', help='Comma-separated list of haplotype identifiers (e.g., A,B)')] = 'A,B',
    mating_scheme: Annotated[str, typer.Option('-m', '--mating-scheme', help='Mating scheme (e.g., RI, F2, HS)')] = 'RI',
    gamma_scale: Annotated[float, typer.Option('-g', '--gamma-scale', help='Scaling factor for transition probabilities')] = 0.01,
    epsilon: Annotated[float, typer.Option('-e', '--epsilon', help='Minimum transition probability')] = 0.000001,
    output_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True, help='Output transition probability file (NPZ format)')] = 'tranprob.npz',
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('get-transition-prob')
    try:
        output_file = str(output_file) if output_file else None

        gbrs_utils.get_transition_prob(
            marker_file=str(marker_file),
            haplotypes=haplotypes,
            mating_scheme=mating_scheme,
            gamma_scale=gamma_scale,
            epsilon=epsilon,
            output_file=output_file
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help='Generate alignment specificity metrics for a sample and set of parental strains')
def get_alignment_spec(
    sample_file: Annotated[Path, typer.Option('-i', '--sample-file', exists=True, dir_okay=False, resolve_path=True, help='Input sample file')],
    haplotypes: Annotated[list[str], typer.Option('-s', '--parental-strains', help='List of parental strain identifiers (e.g., A,B,C,D)')],
    min_expr: Annotated[float, typer.Option('-m', '--min-expr', help='Minimum expression threshold')] = 2.0,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('get_alignment_spec')
    try:
        # haplotype shortcut: the following command line options are all equal
        # -h A,B,C,D,E,F,G,H
        # -h A -h B -h C -h D -h E -h F -h G -h H
        # -h A,B,C,D -h E -h F -h G,H
        all_haplotypes: list[str] = []
        for x in haplotypes:
            all_haplotypes.extend(x.split(','))

        gbrs_utils.get_alignment_spec(
            sample_file=str(sample_file),
            haplotypes=all_haplotypes,
            min_expr=min_expr
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help='Apply genotype calls to a multi-way alignment incidence matrix to produce a genotyped alignment file')
def stencil(
    alignment_file: Annotated[Path, typer.Option('-i', '--alignment-file', exists=True, dir_okay=False, resolve_path=True, help='Input alignment incidence file (HDF5 format)')],
    genotype_file: Annotated[Path, typer.Option('-G', '--genotype', exists=True, dir_okay=False, resolve_path=True, help='Genotype calls by GBRS (TSV format)')],
    group_file: Annotated[Path, typer.Option('-g', '--group-file', exists=True, dir_okay=False, resolve_path=True, help='Gene-to-isoform mapping file (TSV format)')] = None,
    output_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True, help='Output genotyped alignment incidence file (HDF5 format)')] = None,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('stencil')
    try:
        group_file = str(group_file) if group_file else None
        output_file = str(output_file) if output_file else None

        emase_utils.stencil(
            alignment_file=str(alignment_file),
            genotype_file=str(genotype_file),
            group_file=group_file,
            output_file=output_file
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help='Debug or inspect genotype probability file and export')
def debug_genoprob(
    genoprob_file: Annotated[Path, typer.Option('-i', '--genoprob-file', exists=True, dir_okay=False, resolve_path=True, help='Input genotype probability file (NPZ format)')],
    strains: Annotated[list[str], typer.Option('-s', '--strains', help='Strain identifiers (e.g., A,B,C,D). Can specify multiple times or comma-separated')],
    output_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True, help='Output file')] = None,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('debug_genoprobs')
    try:
        # strains shortcut: the following command line options are all equal
        # -s A,B,C,D,E,F,G,H
        # -s A,B,C,D -s E -s F -s G,H
        all_strains: list[str] = []
        for x in strains:
            all_strains.extend(x.split(','))

        output_file = str(output_file) if output_file else None

        gbrs_utils.debug_genoprob(
            genoprob_file=str(genoprob_file),
            strains=all_strains,
            output_file=output_file
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help='Generate synthetic BAM files for testing with realistic read patterns')
def generate_bam(
    output_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True, help='Output BAM file path')],
    num_loci: Annotated[int, typer.Option('-l', '--loci', help='Number of loci/transcripts to generate')] = 10,
    num_haplotypes: Annotated[int, typer.Option('-h', '--haplotypes', help='Number of haplotypes (A, B, C, D, etc.)')] = 4,
    num_read_names: Annotated[int, typer.Option('-r', '--read-names', help='Number of unique read names to generate')] = 100,
    alignments_r1: Annotated[int, typer.Option('-a', '--alignments', help='Number of alignments in first BAM file (R1)')] = 100,
    alignments_r2: Annotated[int, typer.Option('--a2', '--alignments-r2', help='Number of alignments in second BAM file (R2) for paired-end')] = None,
    read_length: Annotated[int, typer.Option('--read-length', help='Length of each read')] = 100,
    coverage_distribution: Annotated[str, typer.Option('--coverage', help='Coverage distribution: uniform, exponential, or realistic')] = 'realistic',
    prefix_locus: Annotated[str, typer.Option('--prefix-locus', help='Prefix for locus names')] = 'ENSMUST',
    prefix_read: Annotated[str, typer.Option('--prefix-read', help='Prefix for read names')] = 'READ',
    custom_haplotypes: Annotated[str, typer.Option('--custom-haplotypes', help='Custom haplotype string (e.g., A,B,C)')] = None,
    randomize_read_names: Annotated[bool, typer.Option('--randomize-reads', help='Randomize read name order')] = False,
    paired_end: Annotated[bool, typer.Option('--paired-end', help='Generate paired-end reads')] = False,
    second_file: Annotated[Path, typer.Option('--second-file', exists=False, dir_okay=False, writable=True, resolve_path=True, help='Second BAM file (R2) for paired-end')] = None,
    locus_file: Annotated[Path, typer.Option('--locus-file', exists=False, dir_okay=False, writable=True, resolve_path=True, help='Output locus file for bam2emase')] = None,
    random_seed: Annotated[int, typer.Option('--seed', help='Random seed for reproducible results')] = None,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help='Increase verbosity (use multiple times for more detail)')] = 0
) -> None:
    """
    Generate synthetic BAM files for testing GBRS pipeline.
    
    This command creates BAM files with configurable parameters to simulate real RNA-seq data.
    The generated files are compatible with the bam2emase pipeline and can be used for testing
    and debugging purposes.
    
    PARAMETER DETAILS:
    
    Basic Parameters:
    - --loci: Number of loci/transcripts to generate (default: 10)
    - --haplotypes: Number of haplotypes A, B, C, D, etc. (default: 4)
    - --read-names: Number of unique read names to generate (default: 100)
    - --alignments: Number of alignments in first BAM file (R1) (default: 100)
    - --a2: Number of alignments in second BAM file (R2) for paired-end (default: same as R1)
    - --read-length: Length of each read in base pairs (default: 100)
    
    Coverage Distribution (--coverage):
    - uniform: Equal probability for all combinations
    - exponential: Weight by locus index (first loci more likely)
    - realistic: Log-normal distribution mimicking real gene expression (default)
    
    Prefix Options:
    - --prefix-locus: Prefix for locus names (default: ENSMUST)
    - --prefix-read: Prefix for read names (default: READ)
    - --custom-haplotypes: Custom haplotype string (e.g., A,B,C) overrides -h
    
    Output Options:
    - --paired-end: Generate paired-end reads (creates R1 and R2 files)
    - --locus-file: Generate locus file for bam2emase compatibility
    
    EXAMPLES:
    
    Basic example:
    gbrs generate-bam -o test.bam -l 3 -h 2 -r 5 -a 10 --coverage realistic
    
    Custom haplotypes:
    gbrs generate-bam -o test.bam -l 3 --custom-haplotypes A,B,C -r 5 -a 15 --coverage uniform
    
    Paired-end with different alignment counts:
    gbrs generate-bam -o test.bam --paired-end -l 10 -h 4 -r 20 -a 50 --a2 45 --prefix-locus GENE --prefix-read SEQ
    """
    logger = utils.configure_logging('gbrs', verbose)
    logger.debug('generate_bam')
    
    try:
        # Validate coverage distribution
        if coverage_distribution not in ['uniform', 'exponential', 'realistic']:
            raise ValueError(f"Invalid coverage distribution: {coverage_distribution}. Must be one of: uniform, exponential, realistic")
        
        # Validate paired-end parameters
        if paired_end:
            if second_file is None:
                # Auto-generate second file name
                output_path = Path(output_file)
                second_file = output_path.parent / f"{output_path.stem}_R2{output_path.suffix}"
            
            if alignments_r2 is None:
                # Use same number of alignments for R2 if not specified
                alignments_r2 = alignments_r1
        
        # Auto-generate locus file if not specified
        if locus_file is None:
            output_path = Path(output_file)
            locus_file = output_path.parent / "loci.txt"
        
        # Generate BAM file using new implementation
        stats = bam_utils.generate_test_bam(
            output_file=str(output_file),
            num_loci=num_loci,
            num_haplotypes=num_haplotypes,
            num_read_names=num_read_names,
            alignments_r1=alignments_r1,
            alignments_r2=alignments_r2,
            read_length=read_length,
            coverage_distribution=coverage_distribution,
            prefix_locus=prefix_locus,
            prefix_read=prefix_read,
            haplotypes=custom_haplotypes,
            randomize_read_names=randomize_read_names,
            paired_end=paired_end,
            second_file=str(second_file) if second_file else None,
            locus_file=str(locus_file),
            random_seed=random_seed
        )
        
        # Print statistics
        logger.info(f"Generated BAM file: {output_file}")
        if paired_end and second_file:
            logger.info(f"Generated paired-end BAM file: {second_file}")
        logger.info(f"Generated locus file: {locus_file}")
        
        logger.info(f"Statistics:")
        logger.info(f"  - Loci: {stats['num_loci']}")
        logger.info(f"  - Haplotypes: {stats['num_haplotypes']} ({', '.join(stats['haplotypes'])})")
        logger.info(f"  - Read names: {stats['num_read_names']}")
        logger.info(f"  - Alignments R1: {stats['alignments_r1']}")
        if paired_end:
            logger.info(f"  - Alignments R2: {stats['alignments_r2']}")
        logger.info(f"  - Read length: {stats['read_length']}")
        logger.info(f"  - Coverage distribution: {stats['coverage_distribution']}")
        logger.info(f"  - Prefix locus: {stats['prefix_locus']}")
        logger.info(f"  - Prefix read: {stats['prefix_read']}")
        logger.info(f"  - Selected combinations: {stats['selected_combinations']}")
        logger.info(f"  - Total references: {stats['total_references']}")
        
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)
        raise typer.Exit(1)



if __name__ == '__main__':
    app()
