=====
Usage
=====

To use GBRS in command line
~~~~~~~~~~~~~~~~~~~~~~~~~~~

First of all, we need to align our RNA-Seq reads against pooled transcriptome of all founder strains::

    bowtie -q -a --best --strata --sam -v 3 ${GBRS_DIR}/bowtie.transcriptome ${FASTQ} \
        | samtools view -bS - > ${BAM_FILE}

Before quantifying multiway allele specificity, bam file should be converted into emase format::

    gbrs bam2emase -i ${BAM_FILE} \
                   -m ${GBRS_DATA}/ref.transcripts.info \
                   -s ${COMMA_SEPARATED_LIST_OF_HAPLOTYPE_CODES} \
                   -o ${EMASE_FILE}

We can compress EMASE format alignment incidence matrix::

    gbrs compress -i ${EMASE_FILE} \
                  -o ${COMPRESSED_EMASE_FILE}

If storage space is tight, you may want to delete ${BAM_FILE} or ${EMASE_FILE} at this point since ${COMPRESSED_EMASE_FILE} has all the information the following steps need. Now we are ready to quantify multiway allele specificity::

    gbrs quantify -i ${COMPRESSED_EMASE_FILE} \
                  -g ${GBRS_DATA}/ref.gene2transcripts.tsv \
                  -L ${GBRS_DATA}/gbrs.hybridized.targets.info \
                  -M 4 --report-counts

Then, we reconstruct the genome based upon gene-level TPM quantities (assuming the sample is a female from the 20th generation Diversity Outbred mice population) ::

    gbrs reconstruct -e gbrs.quantified.multiway.genes.tpm \
                     -t ${GBRS_DATA}/tranprob.DO.G20.F.npz \
                     -x ${GBRS_DATA}/avecs.npz \
                     -g ${GBRS_DATA}/ref.gene_ids.ordered.npz

We can now quantify allele-specific expressions on diploid transcriptome::

    gbrs quantify -i ${COMPRESSED_EMASE_FILE} \
                  -G gbrs.reconstructed.genotypes.tsv \
                  -g ${GBRS_DATA}/ref.gene2transcripts.tsv \
                  -L ${GBRS_DATA}/gbrs.hybridized.targets.info \
                  -M 4 --report-counts

Genotype probabilities are on a grid of genes. For eQTL mapping or plotting genome reconstruction, we may want to interpolate probability on a grid at the genome scale.::

    gbrs interpolate -i gbrs.reconstructed.genoprobs.npz \
                     -g ${GBRS_DATA}/ref.genome_grid.64k.txt \
                     -p ${GBRS_DATA}/ref.gene_pos.ordered.npz \
                     -o gbrs.interpolated.genoprobs.npz

To plot a reconstructed genome::

    gbrs plot -i gbrs.interpolated.genoprobs.npz \
              -o gbrs.plotted.genome.pdf \
              -n ${SAMPLE_ID}


To use GBRS in a project
~~~~~~~~~~~~~~~~~~~~~~~~

All the features are available as a python module. You can simply::

    import gbrs

