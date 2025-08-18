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

## Quick Start (paired-end, ships with repository)

A miniature dataset lives in `examples/example3/`.  The snippet below reproduces the full pipeline; adjust paths if you move the dataset.

```bash
# ------- variables -----------------------------------------------------
S=examples/example3/example          # output prefix
THREADS=8
HAPS=A,B,C,D,E,F,G,H                 # founder list
INDEX=examples/example3/supporting_files/emase_gbrs/rel_2112_v8/bowtie/bowtie.transcripts
META=examples/example3/supporting_files/emase_gbrs/rel_2112_v8/emase.fullTranscripts.info
G2T=examples/example3/supporting_files/emase_gbrs/rel_2112_v8/emase.gene2transcripts.tsv
LEN=examples/example3/supporting_files/emase_gbrs/rel_2112_v8/emase.pooled.fullTranscripts.info
EMISS=examples/example3/supporting_files/emase_gbrs/rel_2112_v8/gbrs_emissions_all_tissues.avecs.npz
TP=examples/example3/supporting_files/emase_gbrs/rel_2112_v8/transition_probabilities/tranprob.DO.G17.M.npz
GPOS=examples/example3/supporting_files/emase_gbrs/rel_2112_v8/ref.gene_pos.ordered_ensBuild_105.npz
GRID=examples/example3/supporting_files/emase_gbrs/rel_2112_v8/ref.genome_grid.GRCm39.tsv
# ----------------------------------------------------------------------

# 1) Align paired-end reads (BAMs already provided; shown for completeness)
# bowtie (see docs/users.md) → ${S}.R1.bam / ${S}.R2.bam

# 2) Convert BAM → EMASE
emase bam2emase -i ${S}.R1.bam -m ${META} -h ${HAPS} -o ${S}.R1.h5
emase bam2emase -i ${S}.R2.bam -m ${META} -h ${HAPS} -o ${S}.R2.h5

# 3) Intersect alignments (pair consistency)
emase get-common-alignments -i ${S}.R1.h5 -i ${S}.R2.h5 -o ${S}.R1R2.h5

# 4) Compress
gbrs compress -i ${S}.R1R2.h5 -o ${S}.R1R2.compressed.h5

# 5) Quantify multi-way expression
gbrs quantify -i ${S}.R1R2.compressed.h5 -g ${G2T} -L ${LEN} -M 4 -a -o ${S}

# 6) Reconstruct genotype
gbrs reconstruct -e ${S}.multiway.genes.tpm -t ${TP} -x ${EMISS} -g ${GPOS} -o ${S}

# 7) Quantify on reconstructed diploid
gbrs quantify -i ${S}.R1R2.compressed.h5 -g ${G2T} -L ${LEN} -G ${S}.genotypes.tsv -M 4 -a -o ${S}

# 8) (optional) Interpolate / plot / export → see docs/users.md
```

*Single-end data?*  Skip step 3 and run `emase bam2emase` once.

---

## Need more detail?

This README is intentionally brief — **see `docs/users.md` for the complete user guide, reference-data specs, command reference, file-format docs, troubleshooting, and more.**

---

MIT License.  Please cite the GBRS paper when publishing research that uses this software.
