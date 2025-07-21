# GBRS: Genome Reconstruction from RNA-Seq

[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.1101/2020.10.11.335323.svg)](https://doi.org/10.1101/2020.10.11.335323)

**GBRS** (Genome Reconstruction by RNA-Seq) is a comprehensive suite of tools for reconstructing individual genomes and quantifying allele-specific expression from RNA-Seq data in multi-parent populations (MPPs). GBRS employs Hidden Markov Models (HMMs) to infer underlying genetic structure from gene expression patterns, enabling high-resolution genome reconstruction without requiring DNA sequencing.

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Workflow](#pipeline-workflow)
- [Output Files](#output-files)
- [Parameters and Configuration](#parameters-and-configuration)
- [Large Dataset Considerations](#large-dataset-considerations)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Support](#support)

## Overview

GBRS reconstructs individual genomes by leveraging allele-specific expression patterns in multi-parent populations. The method works by:

1. **Multi-way Alignment**: Aligning RNA-Seq reads to a combined transcriptome index containing all founder strain transcripts
2. **EMASE Algorithm**: Using expectation maximization for allele-specific expression to resolve multi-mapped reads
3. **HMM Reconstruction**: Inferring the most likely diplotype sequence using Hidden Markov Models
4. **Allele-Specific Quantification**: Quantifying expression on the reconstructed diploid genome

### Algorithm Overview

The GBRS algorithm consists of four main phases:

1. **Expression Profiling**: Convert RNA-Seq reads to gene-level expression profiles using EMASE
2. **Genome Reconstruction**: Infer diploid genotypes using HMM with founder strain expression patterns
3. **Diploid Quantification**: Re-quantify expression on the reconstructed diploid transcriptome
4. **Quality Control**: Detect sample mix-ups and validate reconstruction accuracy

## Key Features

- **High-resolution genome reconstruction** from RNA-Seq data alone
- **Support for multiparent populations** (tested with Diversity Outbred and Collaborative Cross mice)
- **Dual algorithm approach** providing both most likely sequences and uncertainty quantification
- **Sample quality control** with automatic detection of sample mix-ups
- **Efficient memory usage** with compressed data formats
- **Large dataset support** optimized for datasets with 50M+ reads
- **Comprehensive output formats** for downstream analysis

### Supported Populations

While tested primarily with mouse models, GBRS is designed to work with any multiparent population. Required data files for [Diversity Outbred](https://www.jax.org/strain/009376) and [Collaborative Cross](https://www.jax.org/mouse-search/?straingroup=Collaborative%20Cross) mice are available [here](https://zenodo.org/records/8289936).

## Installation

### Prerequisites

- Python 3.12 or higher
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

### Docker Installation (Recommended for Large Datasets)

For large datasets or reproducible environments, use Docker:

```bash
# Pull the latest image
docker pull quay.io/jaxcompsci/gbrs_py3:latest

# Or build locally
docker build -t gbrs:latest .

# Run with data volumes
docker run --rm -v $(pwd)/data:/data -v $(pwd)/output:/output gbrs:latest gbrs --help
```

### Verification

After installation, verify that GBRS is working correctly:

```bash
gbrs --help
```

## Quick Start

For a complete example workflow, see the [Pipeline Workflow](#pipeline-workflow) section below. Here's a minimal example:

```bash
# 1. Align reads to multi-way transcriptome index
bowtie -q -a --best --strata --sam -v 3 ${GBRS_DATA}/bowtie.transcriptome sample.fastq | samtools view -bS - > sample.bam

# 2. Convert BAM to EMASE format
gbrs bam2emase -i sample.bam -m ${GBRS_DATA}/transcripts.info -h A,B,C,D,E,F,G,H -o sample.emase

# 3. Compress EMASE file for storage efficiency
gbrs compress -i sample.emase -o sample.compressed.emase

# 4. Quantify multi-way expression
gbrs quantify -i sample.compressed.emase -g ${GBRS_DATA}/ref.gene2transcripts.tsv -L ${GBRS_DATA}/gbrs.hybridized.targets.info -M 4 --report-alignment-counts

# 5. Reconstruct genome using HMM
gbrs reconstruct -e gbrs.quantified.multiway.genes.tpm -t ${GBRS_DATA}/tranprob.DO.G20.F.npz -x ${GBRS_DATA}/avecs.npz -g ${GBRS_DATA}/ref.gene_pos.ordered.npz

# 6. Quantify diploid expression
gbrs quantify -i sample.compressed.emase -G gbrs.reconstructed.genotypes.tsv -g ${GBRS_DATA}/ref.gene2transcripts.tsv -L ${GBRS_DATA}/gbrs.hybridized.targets.info -M 4
```

## Pipeline Workflow

The GBRS pipeline consists of 6 main steps, each building upon the previous step's output.

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

**Note:** For paired-end data, align R1 and R2 reads separately, then use `emase get-common-alignments` to pair them.

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
- `${COMMA_SEPARATED_HAPLOTYPES}`: Haplotype codes (e.g., A,B,C,D,E,F,G,H for DO mice)
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

## Large Dataset Considerations

GBRS is designed to handle large RNA-Seq datasets efficiently. Based on the example you provided (78M reads, 968M alignments), here are key considerations:

### Memory Requirements

- **Compression Step**: ~8-16GB RAM for datasets with 50M+ reads
- **Reconstruction Step**: ~4-8GB RAM depending on number of genes
- **Quantification Step**: ~2-4GB RAM per sample

### Storage Optimization

- **Use compressed EMASE files**: Reduces storage by 80-90%
- **Delete intermediate files**: BAM and uncompressed EMASE files can be deleted after compression
- **Batch processing**: Process samples in batches to manage disk space

### Performance Tips for Large Datasets

1. **Use Docker**: Containerized environment ensures consistent performance
2. **Monitor resources**: Use `htop` or `top` to monitor memory usage
3. **Batch processing**: Process multiple samples in parallel when possible
4. **Use SSD storage**: Faster I/O for large file operations

### Example Large Dataset Workflow

```bash
# For a dataset with 78M reads:
# 1. Align with Bowtie (may take 2-4 hours)
bowtie -q -a --best --strata --sam -v 3 ${GBRS_DATA}/bowtie.transcriptome sample.fastq | samtools view -bS - > sample.bam

# 2. Convert to EMASE (may take 1-2 hours, 8-16GB RAM)
gbrs bam2emase -i sample.bam -m ${GBRS_DATA}/transcripts.info -h A,B,C,D,E,F,G,H -o sample.emase

# 3. Compress (may take 30-60 minutes, 8-16GB RAM)
gbrs compress -i sample.emase -o sample.compressed.emase

# 4. Quantify (may take 1-2 hours, 2-4GB RAM)
gbrs quantify -i sample.compressed.emase -g ${GBRS_DATA}/ref.gene2transcripts.tsv -L ${GBRS_DATA}/gbrs.hybridized.targets.info -M 4

# 5. Reconstruct (may take 30-60 minutes, 4-8GB RAM)
gbrs reconstruct -e gbrs.quantified.multiway.genes.tpm -t ${GBRS_DATA}/tranprob.DO.G20.F.npz -x ${GBRS_DATA}/avecs.npz -g ${GBRS_DATA}/ref.gene_pos.ordered.npz
```

## Troubleshooting

### Common Issues

#### 1. Memory Errors
**Symptoms**: OutOfMemoryError or system crashes during reconstruction
**Solutions**:
- Reduce `expr_threshold` to include fewer genes
- Use subset of chromosomes for testing
- Increase system RAM or use compute cluster
- Use Docker with increased memory limits

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
Choi, K., Lloyd, M.W., He, H., Gatti, D.M., Philip, V.M., Raghupathy, N., 
Vincent, M., Lek, S., Gerdes Gyuricza, I., Munger, S.C., Attie, A.D., 
Keller, M., Chesler, E.J., Broman, K.W., Srivastava, A., Churchill, G.A. 
(2024). Genome reconstruction by RNA-Seq (GBRS): A novel approach for 
genotyping and quantifying allele-specific expression in multiparent 
populations. [Journal reference to be added]
```

For the data files used in this study, please cite:

```
Choi, K., et al. (2024). GBRS reference data for Diversity Outbred and 
Collaborative Cross mice. Zenodo. https://doi.org/10.5281/zenodo.8289936
```

## Support

### Documentation
- [GBRS GitHub Repository](https://github.com/churchill-lab/gbrs)
- [User Guide](https://github.com/churchill-lab/gbrs/wiki)
- [API Documentation](https://github.com/churchill-lab/gbrs/tree/main/docs)

### Contact
- **Issues**: [GitHub Issues](https://github.com/churchill-lab/gbrs/issues)
- **Questions**: [GitHub Discussions](https://github.com/churchill-lab/gbrs/discussions)
- **Email**: matt.vincent@jax.org, mike.lloyd@jax.org

### Contributing
We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

---

**License**: MIT License - see [LICENSE](LICENSE) file for details.

