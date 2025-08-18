## GBRS User Guide (Final)

This guide is self-contained: if you read only this file, you can run the full GBRS pipeline. Every
command lists what it does, all parameters, the exact inputs and outputs, and a runnable example. 

---

## 0) Prerequisites and installation
- Python 3.12+
- Bowtie ≥ 1.3.1 (multi-way alignment)
- SAMtools ≥ 1.17 (SAM→BAM)

Virtual environment (recommended):
```bash
python -m venv gbrs_env
source gbrs_env/bin/activate
pip install git+https://github.com/churchill-lab/gbrs.git@v1.1.0
```
Docker (alternative):
```bash
docker pull churchilllab/gbrs:1.1.0
```

---

## 1) Required reference data
You need these files before running the pipeline:

### Core alignment files
- **Bowtie index**: `bowtie.transcripts` (multi-way pooled transcriptome)
  - **What it contains**: Indexed transcriptome sequences from all founder strains, with transcript IDs formatted as `TRANSCRIPTID_HAPLOTYPE` (e.g., `ENSMUST00000000001_A`)
  - **Format**: Bowtie index files (`.1.ebwt`, `.2.ebwt`, `.3.ebwt`, `.4.ebwt`, `.rev.1.ebwt`, `.rev.2.ebwt`)
  - **Usage**: Used by Bowtie for RNA-Seq read alignment to the multi-way transcriptome
  - **How to get**: Build from pooled transcriptome using `bowtie-build` (see Appendix B)

### EMASE metadata files
- **Transcript info**: `emase.fullTranscripts.info`
  - **What it contains**: List of transcript IDs that define the loci/transcripts in the multi-way transcriptome
  - **Format**: TSV with columns: `transcript_id`, `metadata` (second column is typically 0.0 and ignored)
  - **Example**: 
    ```
    ENSMUST00000000001      0.0
    ENSMUST00000000003      0.0
    ENSMUST00000000010      0.0
    ```
  - **Usage**: Used by `emase bam2emase` to define which transcripts exist in the multi-way transcriptome and create the locus mapping
  - **Note**: The second column values are ignored - only the transcript IDs matter

- **Gene↔transcript map**: `emase.gene2transcripts.tsv`
  - **What it contains**: Mapping between gene IDs and their transcript IDs for aggregating transcript-level results to gene-level
  - **Format**: TSV with gene ID in first column, followed by all transcript IDs for that gene
  - **Example**:
    ```
    ENSMUSG00000000001      ENSMUST00000000001
    ENSMUSG00000000003      ENSMUST00000114041      ENSMUST00000000003
    ENSMUSG00000000028      ENSMUST00000231819      ENSMUST00000096990      ENSMUST00000000028
    ```
  - **Usage**: Used by `gbrs quantify` to aggregate transcript-level expression to gene-level and enable gene-level analysis

- **Transcript lengths**: `emase.pooled.fullTranscripts.info`
  - **What it contains**: Transcript ID to length mapping for the pooled transcriptome, with haplotype suffixes
  - **Format**: TSV with columns: `transcript_id_haplotype`, `length`
  - **Example**:
    ```
    ENSMUST00000000001_A    3262
    ENSMUST00000000003_A    902
    ENSMUST00000000010_A    2576
    ```
  - **Usage**: Used by `gbrs quantify` for length bias correction in TPM calculation (normalizes by effective transcript length)
  - **Note**: The EMASE algorithm uses these lengths to adjust for transcript length bias in expression estimation

### GBRS model files
- **Emission profiles**: `gbrs_emissions_all_tissues.avecs.npz`
  - **What it contains**: Pre-trained emission model parameters for each gene, representing founder-specific expression patterns
  - **Format**: NPZ file with one key per gene ID, each containing a 2D array of founder expression profiles
  - **Usage**: Used by `gbrs reconstruct` to compute emission probabilities for the HMM based on founder expression patterns
  
- **Transition probabilities**: `tranprob.*.npz`
  - **What it contains**: HMM transition matrices for each chromosome, modeling recombination probabilities between adjacent genes
  - **Format**: NPZ file with one key per chromosome, each containing transition probability matrices
  - **Naming convention**: `tranprob.{POPULATION}.{GENERATION}.{SEX}.npz` (e.g., `tranprob.DO.G17.M.npz`)
  - **Usage**: Used by `gbrs reconstruct` for the forward-backward algorithm to model genetic linkage between genes

- **Gene positions**: `ref.gene_pos.ordered_ensBuild_105.npz`
  - **What it contains**: Chromosome and genetic position (cM) for each gene, ordered by genomic position
  - **Format**: NPZ file with one key per chromosome, each containing arrays of `(gene_id, position_cM)` pairs
  - **Usage**: Used by `gbrs reconstruct` to order genes by position and by `gbrs interpolate` for grid interpolation

### Genome reference files
- **FASTA index**: `ref.fa.fai`
  - **What it contains**: Chromosome names and lengths from the reference genome assembly
  - **Format**: FASTA index format with columns: `chromosome`, `length`, `offset`, `line_length`, `line_width`
  - **Example**:
    ```
    1       195154279       56      60      61
    10      130530862       198406965       60      61
    11      121973369       331113400       60      61
    ```
  - **Usage**: Used throughout the pipeline by `get_chromosome_info()` for:
    - Chromosome ordering and validation in `gbrs reconstruct`
    - Chromosome filtering in `gbrs interpolate` 
    - Chromosome sorting and plotting in `gbrs plot`
  
### Optional files
- **Grid file**: `ref.genome_grid.GRCm39.tsv` (optional, for interpolation)
  - **What it contains**: Uniform grid of positions across the genome for interpolation
  - **Format**: TSV with columns: `marker`, `chr`, `pos`, `cM`, `bp`
  - **Example**:
    ```
    marker  chr     pos     cM      bp
    1_3000000       1       3000000 0.00001000000000001     3000000
    1_3039563       1       3039563 0.02001 3039563
    ```
  - **Usage**: Used by `gbrs interpolate` to create uniform grid positions for plotting/QTL analysis

- **Founder colors**: `founder.hexcolor.info`
  - **What it contains**: Color assignments for each founder strain for visualization
  - **Format**: TSV with columns: `founder_id`, `hex_color`
  - **Example**:
    ```
    A       #F0F000
    B       #808080
    C       #F08080
    D       #1010F0
    ```
  - **Usage**: Used by `gbrs plot` for consistent color coding of founder strains

### Notes
- Mouse DO/CC bundles are published (see paper/Zenodo).
- You can build your own using `emase prepare` and related commands (see Appendix B).
- **Critical**: All files must use the same founder strain order (e.g., `A,B,C,D,E,F,G,H`) throughout the pipeline.

---

## 2) End-to-end quick start (paired-end example)
```bash
S=example
THREADS=8
HAPS=A,B,C,D,E,F,G,H
INDEX=./supporting_files/.../bowtie/bowtie.transcripts
META=./supporting_files/.../emase.fullTranscripts.info
G2T=./supporting_files/.../emase.gene2transcripts.tsv
LEN=./supporting_files/.../emase.pooled.fullTranscripts.info
EMISS=./supporting_files/.../gbrs_emissions_all_tissues.avecs.npz
TP=./supporting_files/.../transition_probabilities/tranprob.DO.G17.M.npz
GPOS=./supporting_files/.../ref.gene_pos.ordered_ensBuild_105.npz
GRID=./supporting_files/.../ref.genome_grid.GRCm39.tsv

# 1) Align (example3 ships BAMs already)
# zcat ${S}_1.fastq.gz | bowtie -p ${THREADS} -q -a --best --strata --sam -v 3 -x ${INDEX} - 2> ${S}.R1.log | samtools view -bS - > ${S}_R1.bam
# zcat ${S}_2.fastq.gz | bowtie -p ${THREADS} -q -a --best --strata --sam -v 3 -x ${INDEX} - 2> ${S}.R2.log | samtools view -bS - > ${S}_R2.bam

# 2) BAM → EMASE
emase bam2emase -i ${S}.R1.bam -m ${META} -h ${HAPS} -o ${S}.R1.h5
emase bam2emase -i ${S}.R2.bam -m ${META} -h ${HAPS} -o ${S}.R2.h5

# 3) Pair consistency
emase get-common-alignments -i ${S}.R1.h5 -i ${S}.R2.h5 -o ${S}.R1R2.h5

# 4) Compress
gbrs compress -i ${S}.R1R2.h5 -o ${S}.R1R2.compressed.h5

# 5) Quantify (multi-way)
gbrs quantify -i ${S}.R1R2.compressed.h5 -g ${G2T} -L ${LEN} -M 4 -a -o ${S}

# 6) Reconstruct
gbrs reconstruct -e ${S}.multiway.genes.tpm -t ${TP} -x ${EMISS} -g ${GPOS} -o ${S}

# 7) Quantify (diploid)
gbrs quantify -i ${S}.R1R2.compressed.h5 -g ${G2T} -L ${LEN} -G ${S}.genotypes.tsv -M 4 -a -o ${S}

# 8) Interpolate (optional)
gbrs interpolate -i ${S}.genoprobs.npz -g ${GRID} -p ${GPOS} -o ${S}.gbrs.interpolated.genoprobs.npz

# 9) Plot
gbrs plot -i ${S}.gbrs.interpolated.genoprobs.npz -o ${S}.gbrs.plotted.genome.pdf -n ${S}

# 10) Export TSV
gbrs export -i ${S}.gbrs.interpolated.genoprobs.npz -s ${HAPS} -g ${GRID} -o ${S}.gbrs.interpolated.genoprobs.tsv
```
Outputs (matching `examples/example3/`)
- Multi-way expression: `example.multiway.*`
- Reconstruction: `example.genoprobs.npz`, `example.genotypes.tsv`, `example.genotypes.npz`
- Diploid expression: `example.diploid.*`
- Interpolated: `example.gbrs.interpolated.genoprobs.npz|tsv`
- Plot: `example.gbrs.plotted.genome.pdf`

---

## 3) Command reference (purpose, parameters, inputs/outputs, examples)

Below, each command is documented with: What it does, Parameters (table), Inputs, Outputs, and an Example.

### emase bam2emase — Convert BAM → EMASE HDF5
**What it does**: Builds a sparse 3‑D incidence matrix (read × locus × allele) from BAM alignments to the multi‑way transcriptome.

Parameters

| Flag | Long | Type | Default | Required | Description |
|------|------|------|---------|----------|-------------|
| -i | --alignment-files | path[,path...] or repeated | — | Yes | Input BAM(s) aligned to pooled transcriptome; BAM refs must be TRANSCRIPTID_<HAP>. |
| -m | --locus-ids | path | — | Yes | Transcript info file `emase.fullTranscripts.info`. |
| -h | --haplotype-char | list[str] | — | Yes | Founder list/order, e.g. `A,B,C,D,E,F,G,H`. |
| -o | --output | path | auto | No | Output EMASE file `.h5`. |
| -d | --delim | str | _ | No | Delimiter between transcript ID and haplotype in BAM refs. |
| — | --index-dtype | str | uint32 | No | Index width: `uint16`/`uint32`/`uint64`. |
| — | --data-dtype | str | uint8 | No | Value dtype; keep `uint8` for incidence. |

Inputs
- BAM alignments (multi‑way index), `emase.fullTranscripts.info`, founder list.

Outputs
- `<out>.h5` EMASE (HDF5) with `/lname`, `/rname`, `/h0..hN/indices,indptr`.

Example
```bash
emase bam2emase -i example.R1.bam -m ./supporting_files/.../emase.fullTranscripts.info -h A,B,C,D,E,F,G,H -o example.R1.h5
```

---

### emase get-common-alignments — Pair-end intersection
**What it does**: Keeps reads that align identically in all provided EMASE files (R1/R2 pairing).

Parameters

| Flag | Long | Type | Default | Required | Description |
|------|------|------|---------|----------|-------------|
| -i | --emase-file | path | — | Yes (2×) | EMASE files to intersect (R1 + R2). |
| -o | --output | path | auto | No | Output EMASE file `.h5`. |

Inputs
- Two (or more) EMASE `.h5` files with identical `/rname`, `/lname`, `hname`.

Outputs
- `merged.h5` containing only reads present with identical locus/allele patterns across inputs.

Example
```bash
emase get-common-alignments -i example.R1.h5 -i example.R2.h5 -o example.R1R2.h5
```

---

### gbrs compress — Compress EMASE HDF5
**What it does**: Reduces size by compressing EMASE alignments without changing logical content.

Parameters

| Flag | Long | Type | Default | Required | Description |
|------|------|------|---------|----------|-------------|
| -i | --emase-file | path | — | Yes | Input EMASE `.h5`. |
| -o | --output | path | — | Yes | Output compressed `.h5`. |
| -c | --comp-lib | str | zlib | No | HDF5 compression library. |

Inputs
- `example.R1R2.h5`.

Outputs
- `example.R1R2.compressed.h5`.

Example
```bash
gbrs compress -i example.R1R2.h5 -o example.R1R2.compressed.h5
```

---

### gbrs quantify — Expression with EMASE (multi‑way or diploid)
**What it does**: Estimates TPM and expected counts. With `-G`, constrains allocation to the reconstructed diploid genome.

Parameters

| Flag | Long | Type | Default | Required | Description |
|------|------|------|---------|----------|-------------|
| -i | --alignment-file | path | — | Yes | EMASE `.h5` (compressed or per‑read). |
| -g | --group-file | path | — | Yes | Gene↔transcript map (`emase.gene2transcripts.tsv`). First col transcript, second col gene. |
| -L | --length-file | path | — | Yes | Transcript lengths (`emase.pooled.fullTranscripts.info`). |
| -G | --genotype | path | — | No | Viterbi calls (`<prefix>.genotypes.tsv`) for diploid mode. |
| -o | --outbase | str | gbrs.quantified | No | Output prefix. |
| -M | --multiread-model | int | 4 | No | EMASE multiread model (1–4). |
| -p | --pseudocount | float | 0.0 | No | Prior read count. |
| -m | --max-iters | int | 999 | No | EM max iterations. |
| -t | --tolerance | float | 0.0001 | No | EM convergence tolerance. |
| -a | --report-alignment-counts | flag | off | No | Write `*.alignment_counts`. |
| -w | --report-posterior | flag | off | No | Write posterior HDF5. |
| -v | --verbose | count | 0 | No | Increase log verbosity (repeatable). |

Inputs
- Multi‑way mode: EMASE `.h5`, gene↔transcript map, lengths file.
- Diploid mode: add `-G <genotypes.tsv>` from `gbrs reconstruct`.

Outputs
- Multi‑way: `<prefix>.multiway.{genes,isoforms}.{tpm,expected_read_counts}` (+ optional `*.alignment_counts`, `*.posterior.h5`).
- Diploid: `<prefix>.diploid.{genes,isoforms}.{tpm,expected_read_counts}` (+ optional reports).

Examples
```bash
# Multi‑way
gbrs quantify -i example.R1R2.compressed.h5 -g ./.../emase.gene2transcripts.tsv -L ./.../emase.pooled.fullTranscripts.info -M 4 -a -o example
# Diploid
gbrs quantify -i example.R1R2.compressed.h5 -g ./.../emase.gene2transcripts.tsv -L ./.../emase.pooled.fullTranscripts.info -G example.genotypes.tsv -M 4 -a -o example
```

---

### gbrs reconstruct — Genotype HMM (forward‑backward + Viterbi)
**What it does**: Infers genotype probabilities and Viterbi diplotype path from multi‑way gene TPM.

Parameters

| Flag | Long | Type | Default | Required | Description |
|------|------|------|---------|----------|-------------|
| -e | --expr-file | path | — | Yes | `<prefix>.multiway.genes.tpm` from quantify. |
| -t | --tprob-file | path | — | Yes | Transition matrices `tranprob.*.npz`. |
| -x | --avec-file | path | — | Yes | Emission profiles `gbrs_emissions*.avecs.npz`. |
| -g | --gpos-file | path | — | Yes | Gene positions NPZ (chrom, gene, cM). |
| -o | --outbase | str | auto | No | Output prefix. |
| -c | --expr-threshold | float | 1.5 | No | Minimum TPM for gene to contribute informative emission. |
| -s | --sigma | float | 0.12 | No | Emission kernel width (smaller = stricter). |
| -v | --verbose | count | 0 | No | Increase verbosity. |

Inputs
- Multi‑way gene TPM, emissions NPZ, transitions NPZ, gene positions NPZ.

Outputs
- `<prefix>.genoprobs.npz` (36×N posteriors per chromosome)
- `<prefix>.genotypes.tsv` (Viterbi diplotype per gene)
- `<prefix>.genotypes.npz` (ordered Viterbi by chromosome)

Example
```bash
gbrs reconstruct -e example.multiway.genes.tpm -t ./.../tranprob.DO.G17.M.npz -x ./.../gbrs_emissions_all_tissues.avecs.npz -g ./.../ref.gene_pos.ordered_ensBuild_105.npz -o example -c 1.5 -s 0.12
```

---

### gbrs interpolate — Gene → grid (optional)
**What it does**: Interpolates gene‑position probabilities onto a uniform grid for plotting/QTL.

Parameters

| Flag | Long | Type | Default | Required | Description |
|------|------|------|---------|----------|-------------|
| -i | --genoprob-file | path | — | Yes | Gene‑based `*.genoprobs.npz`. |
| -g | --grid-file | path | $GBRS_DATA/ref.genome_grid.64k.txt | No | Grid TSV: `grid_id, chr, bp, cM`. |
| -p | --gpos-file | path | $GBRS_DATA/ref.gene_pos.ordered.npz | No | Gene positions NPZ. |
| -o | --output | path | auto | No | Output NPZ filename. |

Inputs
- `example.genoprobs.npz`, grid TSV, gene positions NPZ.

Outputs
- `example.gbrs.interpolated.genoprobs.npz` (36×M per chromosome at grid positions).

Example
```bash
gbrs interpolate -i example.genoprobs.npz -g ./.../ref.genome_grid.GRCm39.tsv -p ./.../ref.gene_pos.ordered_ensBuild_105.npz -o example.gbrs.interpolated.genoprobs.npz
```

---

### gbrs plot — Genome mosaic figure
**What it does**: Draws founder mosaic across chromosomes from genotype probabilities.

Parameters

| Flag | Long | Type | Default | Required | Description |
|------|------|------|---------|----------|-------------|
| -i | --genoprob-file | path | — | Yes | Interpolated or gene‑based NPZ. |
| -o | --output | path | auto | No | Output figure path. |
| -f | --format | str | pdf | No | pdf|png|svg|... |
| -n | --sample_name | str | '' | No | Label for title. |

Inputs
- `example.gbrs.interpolated.genoprobs.npz`.

Outputs
- `example.gbrs.plotted.genome.pdf`.

Example
```bash
gbrs plot -i example.gbrs.interpolated.genoprobs.npz -o example.gbrs.plotted.genome.pdf -n example
```

---

### gbrs export — NPZ → TSV for QTL
**What it does**: Converts diplotype probabilities to founder dosages at each position.

Parameters

| Flag | Long | Type | Default | Required | Description |
|------|------|------|---------|----------|-------------|
| -i | --genoprob-file | path | — | Yes | Interpolated (or gene‑based) NPZ. |
| -s | --strains | list[str] | — | Yes | Founder list; controls column order. |
| -g | --grid-file | path | $GBRS_DATA/ref.genome_grid.64k.txt | No | Grid TSV (for position structure). |
| -o | --output | path | auto | No | Output TSV.

Inputs
- Genoprob NPZ, strains list, grid file.

Outputs
- `example.gbrs.interpolated.genoprobs.tsv` with columns for founder dosages.

Example
```bash
gbrs export -i example.gbrs.interpolated.genoprobs.npz -s A,B,C,D,E,F,G,H -g ./.../ref.genome_grid.GRCm39.tsv -o example.gbrs.interpolated.genoprobs.tsv
```

---

## 4) File formats (what's inside)

### EMASE HDF5 files (`.h5`)
**Structure**: Sparse 3D incidence matrices stored in HDF5 format

**Root attributes**:
- `@shape`: Tuple `(num_loci, num_haplotypes, num_reads)` defining matrix dimensions
- `@hname`: List of haplotype identifiers (e.g., `['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']`)
- `@mtype`: Matrix type (typically `'csc_matrix'`)
- `@incidence_only`: Boolean indicating if file contains only presence/absence data

**Arrays**:
- `/lname`: Array of locus/transcript names (e.g., `['ENSMUST00000000001', 'ENSMUST00000000002', ...]`)
- `/rname`: Array of read names (e.g., `['SRR123456.1', 'SRR123456.2', ...]`) - only in per-read files

**Per-haplotype sparse matrices** (one per founder strain):
- `/h0/`, `/h1/`, `/h2/`, etc.: Each contains:
  - `indices`: CSC matrix indices array (read positions)
  - `indptr`: CSC matrix indptr array (locus boundaries)
  - `data`: CSC matrix data array (if not incidence-only)

**Example structure**:
```
/ (root)
├── @shape: (50000, 8, 1000000)  # 50k transcripts, 8 haplotypes, 1M reads
├── @hname: ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
├── @mtype: 'csc_matrix'
├── @incidence_only: True
├── /lname: ['ENSMUST00000000001', 'ENSMUST00000000002', ...]
├── /rname: ['SRR123456.1', 'SRR123456.2', ...]  # per-read files only
├── /h0:  # Haplotype A
│   ├── indices: [0, 5, 12, 23, ...]  # read positions
│   └── indptr: [0, 2, 5, 8, ...]    # locus boundaries
├── /h1:  # Haplotype B
│   ├── indices: [1, 6, 13, 24, ...]
│   └── indptr: [0, 3, 6, 9, ...]
└── ... (h2-h7 for other haplotypes)
```

**File types**:
- **Per-read files**: Contain `/rname` array, larger size, used for detailed analysis
- **Compressed files**: Remove `/rname`, group identical alignment patterns, much smaller

---

### Expression files (`.tpm`, `.expected_read_counts`, `.alignment_counts`)

#### TPM files (`*.genes.tpm`, `*.isoforms.tpm`)
**What TPM means**: Transcripts Per Million - a normalized expression measure where the sum of all TPM values across all transcripts/genes equals 1,000,000. This normalization accounts for transcript length and sequencing depth, making expression values comparable across samples and transcripts of different lengths.

**Format**: Tab-separated values with transcript/gene IDs, founder-specific TPM values, and a total column

**Structure**:
- Header: `locus` followed by founder haplotypes (A, B, C, D, E, F, G, H) and `total`
- Data: One row per transcript/gene with TPM values per founder and total sum
- **Key property**: The sum of all TPM values across the entire file equals 1,000,000

**Example** (`example.multiway.genes.tpm`):
```
locus	A	B	C	D	E	F	G	H	total
ENSMUSG00000000001	10.725671960377081	10.725671960377081	10.725671960377081	10.725671960377081	0.09539518961440861	48.504229214838546	0.4850186772875907	10.725671960377081	102.71300288362595
ENSMUSG00000000003	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
ENSMUSG00000000028	0.9907943734569702	0.9907943734569702	0.9907943734569702	0.9907943734569702	0.9907943734569702	5.820297633794942	0.9907942138369417	0.13670760068211651	11.9017713155988
```

**Example** (`example.multiway.isoforms.tpm`):
```
locus	A	B	C	D	E	F	G	H	total
ENSMUST00000000001	10.725671960377081	10.725671960377081	10.725671960377081	10.725671960377081	0.09539518961440861	48.504229214838546	0.4850186772875907	10.725671960377081	102.71300288362595
ENSMUST00000000003	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
ENSMUST00000000010	0.059598394399775954	0.021156672138640406	0.021156672138640406	0.021156672138640406	1.3119900993063451	0.03204314277283058	1.3569891849390443	0.021156672138640406	2.845247509972557
```

**TPM interpretation**:
- Each value represents the proportion of transcripts (per million) from that founder strain
- Values are length-normalized and sequencing-depth normalized
- The `total` column shows the sum across all founders for that transcript/gene
- **Global normalization**: Sum of all TPM values in the file = 1,000,000

#### Expected read count files (`*.expected_read_counts`)
**Format**: Same structure as TPM files but with expected read counts instead of TPM

**What it contains**: Expected read counts from the EMASE algorithm after accounting for multi-mapping reads and length bias

**Example** (`example.multiway.genes.expected_read_counts`):
```
locus	A	B	C	D	E	F	G	H	total
ENSMUSG00000000001	1171.7059455000963	1171.7059455000963	1171.7059455000963	1171.7059455000963	10.414679342516747	5302.104139110948	52.95145404605297	1171.7059455000963	11224.0
ENSMUSG00000000003	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
ENSMUSG00000000028	69.88292662201697	69.88292662201697	69.88292662201697	69.88292662201697	69.88292662201697	410.241524254529	69.88292258104079	3.460920054345358	833.0
```

**Interpretation**:
- Values represent expected read counts after EMASE processing
- The `total` column shows the total expected reads for that transcript/gene
- These are the raw counts before TPM normalization

#### Alignment count files (`*.alignment_counts`)
**Format**: Raw alignment counts with detailed breakdown of unique vs. multi-mapping reads

**Structure**:
- Header: `locus`, `aln_A` through `aln_H` (total alignments per founder), `uniq_A` through `uniq_H` (unique alignments per founder), `locus_uniq` (total unique alignments for the locus)

**What it contains**: Raw alignment statistics before EMASE processing, showing both total alignments and unique alignments per founder strain

**Example** (`example.multiway.genes.alignment_counts`):
```
locus	aln_A	aln_B	aln_C	aln_D	aln_E	aln_F	aln_G	aln_H	uniq_A	uniq_B	uniq_C	uniq_D	uniq_E	uniq_F	uniq_G	uniq_H	locus_uniq
ENSMUSG00000000001	9014.0	9014.0	9014.0	9014.0	2584.0	8794.0	3022.0	9014.0	0.0	0.0	0.0	2.0	1726.0	8.0	0.0	11224.0
ENSMUSG00000000003	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
ENSMUSG00000000028	752.0	752.0	752.0	752.0	752.0	751.0	751.0	463.0	0.0	0.0	0.0	0.0	80.0	0.0	1.0	833.0
```

**Interpretation**:
- `aln_X`: Total number of alignments to founder X (includes multi-mapping reads)
- `uniq_X`: Number of uniquely mapping reads to founder X
- `locus_uniq`: Total number of uniquely mapping reads for this locus across all founders
- **Key insight**: Multi-mapping reads appear in multiple founder columns, while unique reads appear only in their true founder

---

### Genotype probability files (`.genoprobs.npz`)

**Format**: Compressed NumPy archive with one key per chromosome

**Structure**: Each chromosome contains a `float64[36, N]` array where:
- **36 rows**: All possible diplotype combinations for 8 founders (8×8 + 8 = 36)
- **N columns**: Number of genes/positions on that chromosome
- **Values**: Log probabilities of each diplotype at each position

**Diplotype ordering**: 
- Rows 0-7: Homozygous (AA, BB, CC, DD, EE, FF, GG, HH)
- Rows 8-35: Heterozygous (AB, AC, AD, AE, AF, AG, AH, BC, BD, BE, BF, BG, BH, CD, CE, CF, CG, CH, DE, DF, DG, DH, EF, EG, EH, FG, FH, GH)

**Example structure**:
```
example.genoprobs.npz
├── 1: float64[36, 2341]    # Chromosome 1: 36 diplotypes × 2341 genes
├── 2: float64[36, 2189]    # Chromosome 2: 36 diplotypes × 2189 genes
├── 3: float64[36, 1956]    # Chromosome 3: 36 diplotypes × 1956 genes
├── X: float64[36, 1567]    # X chromosome: 36 diplotypes × 1567 genes
└── Y: float64[36, 234]     # Y chromosome: 36 diplotypes × 234 genes
```

**Data interpretation**:
- `genoprobs['1'][0, 100]`: Log probability of diplotype AA at gene 100 on chromosome 1
- `genoprobs['1'][8, 100]`: Log probability of diplotype AB at gene 100 on chromosome 1
- Higher (less negative) values = higher probability

---

### Genotype call files (`.genotypes.tsv`)

**Format**: Tab-separated values with gene IDs and called diplotypes

**Structure**:
- Header: `#Gene_ID	Diplotype`
- Data: One row per gene with the most likely diplotype (Viterbi path)

**Example** (`example.genotypes.tsv`):
```
#Gene_ID	Diplotype
ENSMUSG00000000001	DF
ENSMUSG00000000003	CC
ENSMUSG00000000028	CF
ENSMUSG00000000031	GH
```

**Diplotype interpretation**:
- `DF`: Heterozygous with haplotypes D and F
- `CC`: Homozygous for haplotype C
- `CF`: Heterozygous with haplotypes C and F
- `GH`: Heterozygous with haplotypes G and H

---

### Genotype NPZ files (`.genotypes.npz`)

**Format**: Compressed NumPy archive with Viterbi diplotype path per chromosome

**Structure**: Each chromosome contains an array of diplotype indices corresponding to the most likely genotype at each position

**Example structure**:
```
example.genotypes.npz
├── 1: int64[2341]    # Chromosome 1: diplotype indices for 2341 genes
├── 2: int64[2189]    # Chromosome 2: diplotype indices for 2189 genes
├── 3: int64[1956]    # Chromosome 3: diplotype indices for 1956 genes
├── X: int64[1567]    # X chromosome: diplotype indices for 1567 genes
└── Y: int64[234]     # Y chromosome: diplotype indices for 234 genes
```

**Data interpretation**:
- `genotypes['1'][100]`: Diplotype index for gene 100 on chromosome 1
- Index 0 = diplotype AA, index 8 = diplotype AB, etc.
- Maps directly to the row indices in the corresponding `genoprobs.npz` file

---

### Interpolated genotype files (`.gbrs.interpolated.genoprobs.npz`)

**Format**: Same structure as `.genoprobs.npz` but interpolated to uniform grid positions

**Structure**: Each chromosome contains a `float64[36, M]` array where:
- **36 rows**: Same diplotype combinations as gene-based files
- **M columns**: Number of grid positions (typically 64k positions per chromosome)
- **Values**: Interpolated log probabilities at each grid position

**Example structure**:
```
example.gbrs.interpolated.genoprobs.npz
├── 1: float64[36, 64000]    # Chromosome 1: 36 diplotypes × 64k grid positions
├── 2: float64[36, 64000]    # Chromosome 2: 36 diplotypes × 64k grid positions
├── 3: float64[36, 64000]    # Chromosome 3: 36 diplotypes × 64k grid positions
└── X: float64[36, 64000]    # X chromosome: 36 diplotypes × 64k grid positions
```

**Usage**: Used for plotting and QTL analysis where uniform position spacing is required

---

### Export TSV files (`.gbrs.interpolated.genoprobs.tsv`)

**Format**: Tab-separated values with founder dosages at each grid position

**Structure**:
- Header: Founder strain columns (A, B, C, D, E, F, G, H)
- Data: One row per grid position with founder dosages (0.0, 1.0, or 2.0)
- **Note**: This file does NOT contain genomic coordinates - it's just the dosage matrix

**Example** (`example.gbrs.interpolated.genoprobs.tsv`):
```
# A	B	C	D	E	F	G	H
0.500000	0.000000	0.000000	0.000000	0.000000	0.500000	0.000000	0.000000
0.500000	0.000000	0.000000	0.000000	0.000000	0.500000	0.000000	0.000000
0.500000	0.000000	0.000000	0.000000	0.000000	0.500000	0.000000	0.000000
0.500000	0.000000	0.000000	0.000000	0.000000	0.500000	0.000000	0.000000
```

**Dosage interpretation**:
- `0.0`: No copies of this founder's haplotype
- `1.0`: One copy of this founder's haplotype  
- `2.0`: Two copies of this founder's haplotype (homozygous)
- **Key insight**: Each row sums to 2.0 (diploid genome), representing the founder composition at that position

**Usage**: Direct input for QTL analysis software (e.g., R/qtl2) where founder dosages are required. The file contains only the dosage matrix without genomic coordinates, so you need to combine it with position information from the grid file for complete analysis.

---

## 5) Common pitfalls and fixes
- “No alignments” → confirm Bowtie index path and that BAM refs use transcriptID_haplotype.
- “Haplotype mismatch” → use the same founder list order in all steps.
- Fragmented Viterbi path → wrong transition file for your population/generation.
- Empty emission for many genes → increase `-c/--expr-threshold` or check group/length files.

---

## Appendix A: Example3 mapping
From `examples/example3/run.sh`, each step and outputs:
- bam2emase → `example.R1.h5`, `example.R2.h5` (or directly `example.R1R2.h5`)
- get-common-alignments → `example.R1R2.h5`
- compress → `example.R1R2.compressed.h5`
- quantify (multi-way) → `example.multiway.*`
- reconstruct → `example.genoprobs.npz`, `example.genotypes.tsv`, `example.genotypes.npz`
- quantify (diploid) → `example.diploid.*`
- interpolate → `example.gbrs.interpolated.genoprobs.npz`
- plot → `example.gbrs.plotted.genome.pdf`
- export → `example.gbrs.interpolated.genoprobs.tsv`

---

## Appendix B: Building references (brief)
- Build pooled transcriptome and index: `emase prepare` or `create-hybrid` (see `src/gbrs/emase/emase_utils.py`).
- Train emissions: see `get_alignment_spec` in `gbrs_utils` (requires founder TPMs).
- Generate transitions: `gbrs get-transition-prob -i markers.tsv -s A,B,... -m RI|F2|CC|DO -o tranprob.npz`.

For full details, see the paper (`GBRS.md`).
