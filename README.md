# GBRS

GBRS is a suite of tools for reconstructing genomes using RNA-Seq data from multiparent population and quantifying allele specific expression.  Although we tested it with mouse models only, GBRS should work for any multiparent populations. For the [Diversity Outbred](https://www.jax.org/strain/009376) and [Collaborative Cross](https://www.jax.org/mouse-search/?straingroup=Collaborative%20Cross) mice, the required data files are available [here](https://zenodo.org/records/8289936).

## Installation

*Note: To avoid conflicts among dependencies, we highly recommend using a [Python virtual environment](https://realpython.com/python-virtual-environments-a-primer/).*

GBRS requires Python 3+ to run.  Install GBRS and all its dependencies from the command line:

```
pip install git+https://github.com/churchill-lab/gbrs
```

The **gbrs** script should now be installed and you should be able to run GBRS from the command line. 

## Usage

*Note: In the steps below,* `{GBRS_DATA}` *refers to a local GBRS directory that contains numerous supporting files to run GBRS.  You can create your own or download them [here](https://zenodo.org/records/8289936).*

##### Step 1: Map reads with BOWTIE and convert SAM to BAM  
#
*Note: R1 and R2 are mapped separately for paired-end data.*

The first step is to align our RNA-Seq reads against the pooled transcriptome of all founder strains:

```
bowtie \
        -q -a --best --strata --sam \
        -v 3 ${GBRS_DATA}/bowtie.transcriptome ${FASTQ} \
    | samtools view -bS - > ${BAM_FILE}
```

where:
`${FASTQ}` is the FASTQ file you are aligning
`${BAM_FILE}` is the resulting BAM file

##### Step 2: Convert BAM to EMASE
#
*Note: R1 and R2 are converted separately for paired-end data.*  

Before quantifying multiway allele specificity, bam file should be converted into emase format:

```
gbrs bam2emase \
        -i ${BAM_FILE} \
        -m ${GBRS_DATA}\transcripts.info \
        -h ${COMMA_SEPARATED_LIST_OF_HAPLOTYPE_CODES} \
        -o ${EMASE_FILE}
```

where:
`${BAM_FILE}` is the output BAM file from STEP 1
`${COMMA_SEPARATED_LIST_OF_HAPLOTYPE_CODES}` is a list of haplotypes (i.e, A,B,C,D,E,F,G,H)
`${EMASE_FILE}` is the resulting EMASE output file

##### Step 3: Single-end data: Compress EMASE file. 
##### Step 3: Paired-end data: Find common alignments and compress EMASE file.  
#
We can compress EMASE format alignment incidence matrix:

```
gbrs compress \
        -i ${EMASE_FILE} \
        -o ${COMPRESSED_EMASE_FILE}
```

where:
`${EMASE_FILE}` is the output EMASE file from STEP 2
`${COMPRESSED_EMASE_FILE}` is the resulting **compressed** EMASE output file

If storage space is tight, you may want to delete `${BAM_FILE}` or `${EMASE_FILE}` at this point since `${COMPRESSED_EMASE_FILE}` has all the information the following steps would need. If you want to merge emase format files in order to, for example, pool technical replicates, you run ‘compress’ once more listing files you want to merge with commas:

```
gbrs compress \
        -i ${COMPRESSED_EMASE_FILE1},${COMPRESSED_EMASE_FILE2},... \
        -o ${MERGED_COMPRESSED_EMASE_FILE}
```

and use `${MERGED_COMPRESSED_EMASE_FILE}` in the following steps. Now we are ready to quantify multiway allele specificity.

##### Step 4: Quantify multiway expression  
#
```
gbrs quantify \
        -i ${COMPRESSED_EMASE_FILE} \
        -g ${GBRS_DATA}/ref.gene2transcripts.tsv \
        -L ${GBRS_DATA}/gbrs.hybridized.targets.info \
        -M 4 \
        --report-alignment-counts
```

##### Step 5: Genotype reconstruction
#

Then, we reconstruct the genome based upon gene-level TPM quantities (assuming the sample is a female from the 20th generation Diversity Outbred mice population)

```
gbrs reconstruct \
        -e gbrs.quantified.multiway.genes.tpm \
        -t ${GBRS_DATA}/tranprob.DO.G20.F.npz \
        -x ${GBRS_DATA}/avecs.npz \
        -g ${GBRS_DATA}/ref.gene_pos.ordered.npz
```

##### Step 6: Quantify diploid expression with GBRS   
#

We can now quantify allele-specific expressions on diploid transcriptome:

```
gbrs quantify \
        -i ${COMPRESSED_EMASE_FILE} \
        -G gbrs.reconstructed.genotypes.tsv \
        -g ${GBRS_DATA}/ref.gene2transcripts.tsv \
        -L ${GBRS_DATA}/gbrs.hybridized.targets.info \
        -M 4 \
        --report-alignment-counts
```

##### Step 7: Interpolate genotypes and genotype probabilities  
#

Genotype probabilities are on a grid of genes. For eQTL mapping or plotting genome reconstruction, we may want to interpolate probability on a decently-spaced grid of the reference genome:

```
gbrs interpolate \
        -i gbrs.reconstructed.genoprobs.npz \
        -g ${GBRS_DATA}/ref.genome_grid.69k.txt \
        -p ${GBRS_DATA}/ref.gene_pos.ordered.npz \
        -o gbrs.interpolated.genoprobs.npz
```

##### Step 8: Plot inferred genotypes  
#
To plot a reconstructed genome:

```
gbrs plot \
        -i gbrs.interpolated.genoprobs.npz \
        -o gbrs.plotted.genome.pdf \
        -n ${SAMPLE_ID}
```

##### Step 9: Export genotype probabilities
#

```
gbrs export \
        -i ${interpolated_genoprobs} \
        -s ${params.gbrs_strain_list} \
        -g ${params.genotype_grid} \
        -o ${sampleID}.gbrs.interpolated.genoprobs.tsv
```


