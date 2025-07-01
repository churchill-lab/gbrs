# GBRS: Genome Reconstruction from RNA-Seq

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.1101/2020.10.11.335323.svg)](https://doi.org/10.1101/2020.10.11.335323)

**GBRS** (Genome Reconstruction from RNA-Seq) is a comprehensive suite of tools for reconstructing genomes using RNA-Seq data from multiparent populations and quantifying allele-specific expression. GBRS employs Hidden Markov Models (HMMs) to infer underlying genetic structure from gene expression patterns, enabling high-resolution genome reconstruction without requiring DNA sequencing.

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Workflow](#pipeline-workflow)
- [Output Files](#output-files)
- [Parameters and Configuration](#parameters-and-configuration)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Support](#support)

## Overview

GBRS reconstructs genomes by leveraging allele-specific expression patterns in multiparent populations. The pipeline consists of two main phases:

1. **Expression Quantification**: Quantify allele-specific expression using the EMASE algorithm
2. **Genome Reconstruction**: Infer the most likely diplotype sequence using HMM algorithms

### Key Features

- **High-resolution genome reconstruction** from RNA-Seq data alone
- **Support for multiparent populations** (tested with Diversity Outbred and Collaborative Cross mice)
- **Dual algorithm approach** providing both most likely sequences and uncertainty quantification
- **Comprehensive output formats** for downstream analysis
- **Efficient memory usage** with compressed data formats

### Supported Populations

While tested primarily with mouse models, GBRS is designed to work with any multiparent population. Required data files for [Diversity Outbred](https://www.jax.org/strain/009376) and [Collaborative Cross](https://www.jax.org/mouse-search/?straingroup=Collaborative%20Cross) mice are available [here](https://zenodo.org/records/8289936).

## Installation

### Prerequisites

- Python 3.8 or higher
- [Bowtie](http://bowtie-bio.sourceforge.net/) for read alignment
- [SAMtools](http://samtools.sourceforge.net/) for BAM file manipulation

### Recommended Setup

We strongly recommend using a Python virtual environment to avoid dependency conflicts:

```bash
# Create virtual environment
python -m venv gbrs_env

# Activate virtual environment
# On macOS/Linux:
source gbrs_env/bin/activate
# On Windows:
gbrs_env\Scripts\activate

# Install GBRS
pip install git+https://github.com/churchill-lab/gbrs
```

### Verification

After installation, verify that GBRS is working correctly:

```bash
gbrs --help
```

## Quick Start

For a complete example workflow, see the [Pipeline Workflow](#pipeline-workflow) section below. Here's a minimal example:

```bash
# 1. Align reads
bowtie -q -a --best --strata --sam -v 3 ${GBRS_DATA}/bowtie.transcriptome sample.fastq | samtools view -bS - > sample.bam

# 2. Convert to EMASE format
gbrs bam2emase -i sample.bam -m ${GBRS_DATA}/transcripts.info -h A,B,C,D,E,F,G,H -o sample.emase

# 3. Compress EMASE file
gbrs compress -i sample.emase -o sample.compressed.emase

# 4. Quantify expression
gbrs quantify -i sample.compressed.emase -g ${GBRS_DATA}/ref.gene2transcripts.tsv -L ${GBRS_DATA}/gbrs.hybridized.targets.info -M 4 --report-alignment-counts

# 5. Reconstruct genome
gbrs reconstruct -e gbrs.quantified.multiway.genes.tpm -t ${GBRS_DATA}/tranprob.DO.G20.F.npz -x ${GBRS_DATA}/avecs.npz -g ${GBRS_DATA}/ref.gene_pos.ordered.npz
```

## Pipeline Workflow

The GBRS pipeline consists of 9 main steps, each building upon the previous step's output.

### Step 1: Read Alignment

Align RNA-Seq reads against the pooled transcriptome of all founder strains using Bowtie.

```bash
bowtie \
    -q -a --best --strata --sam \
    -v 3 ${GBRS_DATA}/bowtie.transcriptome ${FASTQ} \
| samtools view -bS - > ${BAM_FILE}
```

**Parameters:**
- `${FASTQ}`: Input FASTQ file
- `${BAM_FILE}`: Output BAM file
- `${GBRS_DATA}`: Directory containing GBRS reference files

**Note:** For paired-end data, align R1 and R2 reads separately.

### Step 2: BAM to EMASE Conversion

Convert BAM files to EMASE format for allele-specific expression analysis.

```bash
gbrs bam2emase \
    -i ${BAM_FILE} \
    -m ${GBRS_DATA}/transcripts.info \
    -h ${COMMA_SEPARATED_HAPLOTYPES} \
    -o ${EMASE_FILE}
```

**Parameters:**
- `${BAM_FILE}`: Input BAM file from Step 1
- `${COMMA_SEPARATED_HAPLOTYPES}`: Haplotype codes (e.g., A,B,C,D,E,F,G,H)
- `${EMASE_FILE}`: Output EMASE file

### Step 3: EMASE Compression

Compress EMASE files to reduce storage requirements and enable merging of replicates.

```bash
# Single file compression
gbrs compress \
    -i ${EMASE_FILE} \
    -o ${COMPRESSED_EMASE_FILE}

# Merge multiple files (e.g., technical replicates)
gbrs compress \
    -i ${COMPRESSED_EMASE_FILE1},${COMPRESSED_EMASE_FILE2},... \
    -o ${MERGED_COMPRESSED_EMASE_FILE}
```

**Storage Optimization:** After compression, you can safely delete the original BAM and EMASE files.

### Step 4: Multiway Expression Quantification

Quantify allele-specific expression using the EMASE algorithm.

```bash
gbrs quantify \
    -i ${MERGED_COMPRESSED_EMASE_FILE} \
    -g ${GBRS_DATA}/ref.gene2transcripts.tsv \
    -L ${GBRS_DATA}/gbrs.hybridized.targets.info \
    -M 4 \
    --report-alignment-counts
```

**Output Files:**
- `gbrs.quantified.multiway.genes.tpm`: Gene-level TPM values
- `gbrs.quantified.multiway.genes.expected_read_counts`: Expected read counts
- `gbrs.quantified.multiway.genes.alignment_counts`: Raw alignment counts

### Step 5: Genome Reconstruction

Reconstruct the genome using HMM algorithms based on gene-level TPM quantities.

```bash
gbrs reconstruct \
    -e gbrs.quantified.multiway.genes.tpm \
    -t ${GBRS_DATA}/tranprob.DO.G20.F.npz \
    -x ${GBRS_DATA}/avecs.npz \
    -g ${GBRS_DATA}/ref.gene_pos.ordered.npz
```

**Output Files:**
- `*.genotypes.tsv`: Tab-separated genotype calls
- `*.genotypes.npz`: Viterbi algorithm output (recommended for recombination analysis)
- `*.genoprobs.npz`: Forward-backward algorithm output (uncertainty quantification)

**Important:** See [Output Files](#output-files) section for detailed explanation of differences between these files.

### Step 6: Diploid Expression Quantification

Quantify allele-specific expression on the reconstructed diploid transcriptome.

```bash
gbrs quantify \
    -i ${MERGED_COMPRESSED_EMASE_FILE} \
    -G gbrs.reconstructed.genotypes.tsv \
    -g ${GBRS_DATA}/ref.gene2transcripts.tsv \
    -L ${GBRS_DATA}/gbrs.hybridized.targets.info \
    -M 4 \
    --report-alignment-counts
```

### Step 7: Genotype Interpolation

Interpolate genotype probabilities to a standardized genomic grid for cross-sample comparison.

```bash
gbrs interpolate \
    -i gbrs.reconstructed.genoprobs.npz \
    -g ${GBRS_DATA}/ref.genome_grid.69k.txt \
    -p ${GBRS_DATA}/ref.gene_pos.ordered.npz \
    -o gbrs.interpolated.genoprobs.npz
```

### Step 8: Genome Visualization

Generate publication-quality plots of the reconstructed genome.

```bash
gbrs plot \
    -i gbrs.interpolated.genoprobs.npz \
    -o gbrs.plotted.genome.pdf \
    -n ${SAMPLE_ID}
```

### Step 9: Data Export

Export genotype probabilities in standard formats for downstream analysis.

```bash
gbrs export \
    -i ${interpolated_genoprobs} \
    -s ${gbrs_strain_list} \
    -g ${genotype_grid} \
    -o ${sampleID}.gbrs.interpolated.genoprobs.tsv
```

## Output Files

### Genotype Reconstruction Outputs

The `gbrs reconstruct` command generates three distinct output files, each serving different analytical purposes:

#### 1. `*.genotypes.tsv` - Tab-Separated Genotype Calls
- **Format**: Standard tab-separated values (TSV)
- **Structure**: Two columns (`Gene_ID`, `Diplotype`)
- **Content**: Most likely diplotype for each gene
- **Use Case**: Easy-to-read format for downstream analysis, input for diploid quantification
- **Example**:
  ```
  Gene_ID    Diplotype
  ENSMUSG00000000001    AB
  ENSMUSG00000000002    CH
  ENSMUSG00000000003    DE
  ```

#### 2. `*.genotypes.npz` - Viterbi Algorithm Output
- **Format**: Compressed NumPy array
- **Structure**: `{chromosome: array_of_diplotypes}`
- **Content**: Most likely diplotype sequence per chromosome
- **Algorithm**: Viterbi algorithm (global optimization)
- **Use Case**: **Recommended for recombination counting and sequence analysis**
- **Advantage**: Provides coherent, biologically plausible sequences

#### 3. `*.genoprobs.npz` - Forward-Backward Algorithm Output
- **Format**: Compressed NumPy array
- **Structure**: `{chromosome: array_of_shape(36, num_genes)}`
- **Content**: Posterior probabilities for all 36 possible diplotypes at each gene position
- **Algorithm**: Forward-backward algorithm (local probability estimation)
- **Use Case**: **Recommended for uncertainty assessment and probability-based analysis**
- **Advantage**: Shows confidence levels and uncertainty in each call

### Understanding Algorithm Differences

**Important**: The genotypes and genoprobs files may show different calls for the same position. This is **expected behavior** and reflects fundamental differences between the two HMM algorithms:

#### Viterbi Algorithm (genotypes.npz)
- **Objective**: Find the single most likely path through all states across the entire sequence
- **Optimization**: Global best path considering transition probabilities between adjacent positions
- **Output**: Coherent sequence respecting recombination patterns
- **When to Use**: Recombination counting, sequence analysis, applications requiring single best estimate

#### Forward-Backward Algorithm (genoprobs.npz)
- **Objective**: Compute posterior probabilities for each possible diplotype at each gene position
- **Optimization**: Local probability estimation at each position independently
- **Output**: Probability distribution showing uncertainty in each call
- **When to Use**: Uncertainty assessment, probability-based analysis, confidence evaluation

#### Example of Algorithm Differences
```
Position 2140:
- Viterbi (genotypes.npz): 'AB' 
- Forward-Backward (genoprobs.npz): 'AC' (prob: 0.3581 vs 0.3441 for AB)
```

The Viterbi algorithm may choose a less probable state at one position to maintain consistency with neighboring positions and achieve a more probable overall sequence.

### Expression Quantification Outputs

- **`*.genes.tpm`**: Gene-level TPM values for each founder strain
- **`*.genes.expected_read_counts`**: Expected read counts for each founder strain
- **`*.genes.alignment_counts`**: Raw alignment counts for each founder strain
- **`*.isoforms.tpm`**: Transcript-level TPM values (transcript-level analysis)

### Interpolated and Visualization Outputs

- **`*.interpolated.genoprobs.npz`**: Genotype probabilities interpolated to standardized genomic grid
- **`*.plotted.genome.pdf`**: Publication-quality visualization of reconstructed genome

## Parameters and Configuration

### Reconstruction Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `expr_threshold` | 1.5 | Minimum expression threshold for genes to be included in reconstruction |
| `sigma` | 0.12 | Controls emission probability precision in the HMM |
| `haplotypes` | A,B,C,D,E,F,G,H | List of founder strain codes |

### Quantification Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `multiread_model` | 4 | EMASE model for handling multi-mapping reads |
| `pseudocount` | 0.0 | Prior read count for allele specificity estimation |
| `tolerance` | 0.0001 | Convergence tolerance for EM algorithm |
| `max_iters` | 999 | Maximum iterations for EM algorithm |

### Performance Optimization

- **Memory Usage**: Large datasets may require significant RAM. Consider reducing gene set or expression threshold
- **Storage**: Use compressed EMASE files to minimize storage requirements
- **Testing**: Run reconstruction on gene subsets for initial testing and parameter optimization

## Troubleshooting

### Common Issues

#### 1. Memory Errors
**Symptoms**: OutOfMemoryError or system crashes during reconstruction
**Solutions**:
- Reduce `expr_threshold` to include fewer genes
- Use subset of chromosomes for testing
- Increase system RAM or use compute cluster

#### 2. File Format Errors
**Symptoms**: "File not found" or "Invalid format" errors
**Solutions**:
- Verify all input files exist and are properly formatted
- Check file permissions
- Ensure correct file paths in `${GBRS_DATA}`

#### 3. Missing Dependencies
**Symptoms**: "Command not found" errors
**Solutions**:
- Install Bowtie: `conda install bowtie` or download from source
- Install SAMtools: `conda install samtools` or download from source
- Verify Python virtual environment activation

#### 4. Low-Quality Reconstructions
**Symptoms**: Poor concordance with known genotypes or excessive recombination
**Solutions**:
- Increase `expr_threshold` to include only high-expression genes
- Adjust `sigma` parameter (lower values = higher precision)
- Verify transition probability files match your population and generation

### Performance Tips

1. **Use compressed EMASE files** to reduce storage requirements
2. **Test on chromosome subsets** before full genome reconstruction
3. **Monitor memory usage** during reconstruction step
4. **Use appropriate transition probability files** for your population and generation
5. **Consider parallel processing** for multiple samples

### Getting Help

If you encounter issues not covered here:

1. Check the [GBRS documentation](https://github.com/churchill-lab/gbrs)
2. Review error messages carefully for specific guidance
3. Contact the development team with detailed error information
4. Include system specifications and parameter settings in bug reports

## Citation

If you use GBRS in your research, please cite:

```
[GBRS publication reference - please add the specific citation from the publication you mentioned]
```

For the data files used in this study, please cite:

```
[Zenodo reference for data files]
```

## Support

### Documentation
- [GBRS GitHub Repository](https://github.com/churchill-lab/gbrs)
- [User Guide](https://github.com/churchill-lab/gbrs/wiki)
- [API Documentation](https://github.com/churchill-lab/gbrs/tree/main/docs)

### Contact
- **Issues**: [GitHub Issues](https://github.com/churchill-lab/gbrs/issues)
- **Questions**: [GitHub Discussions](https://github.com/churchill-lab/gbrs/discussions)
- **Email**: [Development Team Contact]

### Contributing
We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

---

**License**: MIT License - see [LICENSE](LICENSE) file for details.

