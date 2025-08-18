# GBRS: Genome Reconstruction by RNA-Seq

[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/) [![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE) [![DOI](https://zenodo.org/badge/DOI/10.1101/2020.10.11.335323.svg)](https://doi.org/10.1101/2020.10.11.335323)

**GBRS** (Genome Reconstruction by RNA-Seq) reconstructs individual genomes and quantifies allele-specific expression directly from RNA-Seq data in multi-parent populations.  For theory and benchmarks see the [GBRS paper](https://www.biorxiv.org/content/10.1101/2020.10.11.335323v3).

---

## Prerequisites

* Python ≥ **3.12**
* Bowtie ≥ **1.3.1**
* SAMtools ≥ **1.17**

## Installation

```bash
# Latest GBRS from GitHub main
pip install git+https://github.com/churchill-lab/gbrs

# – or – Reproducible Docker image (no local deps)
docker pull quay.io/jaxcompsci/gbrs_py3:latest
```

---

## Quick Start (paired-end example)

Below is a minimal, end-to-end run using **your own RNA-Seq FASTQs** together with the publicly available supporting-files bundle (download from [Zenodo 10.5281/zenodo.8289936](https://zenodo.org/records/8289936)).

```bash
# ---- paths --------------------------------------------------------------
FASTQ_R1=mySample_R1.fastq.gz          # change to your paired-end files
FASTQ_R2=mySample_R2.fastq.gz
THREADS=8
HAPS=A,B,C,D,E,F,G,H                  # founder order

# Directory created after unpacking the Zenodo archive
export GBRS_DATA=/path/to/gbrs_supporting_files
# -------------------------------------------------------------------------

# 1) Align reads to the pooled transcriptome (R1 / R2 separately)
zcat ${FASTQ_R1} | bowtie -p ${THREADS} -q -a --best --strata --sam -v 3 \
      ${GBRS_DATA}/bowtie.transcriptome - \
  2> mySample.R1.log | samtools view -bS - > mySample.R1.bam

zcat ${FASTQ_R2} | bowtie -p ${THREADS} -q -a --best --strata --sam -v 3 \
      ${GBRS_DATA}/bowtie.transcriptome - \
  2> mySample.R2.log | samtools view -bS - > mySample.R2.bam

# 2) Convert BAM → EMASE
emase bam2emase -i mySample.R1.bam -m ${GBRS_DATA}/emase.fullTranscripts.info \
                -h ${HAPS} -o mySample.R1.h5
emase bam2emase -i mySample.R2.bam -m ${GBRS_DATA}/emase.fullTranscripts.info \
                -h ${HAPS} -o mySample.R2.h5

# 3) Intersect paired-end alignments & compress
emase get-common-alignments -i mySample.R1.h5 -i mySample.R2.h5 \
                            -o mySample.R1R2.h5
gbrs compress -i mySample.R1R2.h5 -o mySample.R1R2.compressed.h5

# 4) Quantify multi-way expression
gbrs quantify -i mySample.R1R2.compressed.h5 \
              -g ${GBRS_DATA}/emase.gene2transcripts.tsv \
              -L ${GBRS_DATA}/emase.pooled.fullTranscripts.info \
              -M 4 -a -o mySample

# 5) Reconstruct genotypes
gbrs reconstruct -e mySample.multiway.genes.tpm \
                 -t ${GBRS_DATA}/transition_probabilities/tranprob.DO.G20.F.npz \
                 -x ${GBRS_DATA}/gbrs_emissions_all_tissues.avecs.npz \
                 -g ${GBRS_DATA}/ref.gene_pos.ordered_ensBuild_105.npz \
                 -o mySample

# 6) Quantify on reconstructed diploid genome
gbrs quantify -i mySample.R1R2.compressed.h5 \
              -g ${GBRS_DATA}/emase.gene2transcripts.tsv \
              -L ${GBRS_DATA}/emase.pooled.fullTranscripts.info \
              -G mySample.genotypes.tsv -M 4 -a -o mySample

# 7) (optional) Interpolate / plot / export  → see docs/users.md
```

*Single-end data?*  Run `emase bam2emase` once, skip the `get-common-alignments` step, and continue from compression onward.

---

## Need more detail?

This README is intentionally brief — **see `docs/users.md` for the complete user guide, reference-data specs, command reference, file-format docs, troubleshooting, and more.**

---

MIT License.  Please cite the GBRS paper when publishing research that uses this software.
