=====
Usage
=====

To use GBRS in command line
~~~~~~~~~~~~~~~~~~~~~~~~~~~

First of all, we need to align our RNA-Seq reads against pooled transcriptome of all founder strains::

    bowtie -q -a --best --strata --sam -v 3 ${GBRS_DIR}/bowtie.pooled.transcriptome ${FASTQ} \
        | samtools view -bS - > ${BAM_FILE}

Before quantifying multiway allele specificity, bam file should be converted to emase format::

    gbrs bam2emase -i ${BAM_FILE} -o ${EMASE_FILE}

We can compress EMASE format alignment incidence matrix::

    gbrs compress -i ${EMASE_FILE} -o ${COMPRESSED_EMASE_FILE}

Now we are ready to quantify multiway allele specificity::

    gbrs quantify -i ${COMPRESSED_EMASE_FILE} \
                  -g ${GBRS_DATA}/REF/emase.gene2transcripts.tsv \
                  -L ${GBRS_DATA}/8-way/emase.pooled.transcripts.info \
                  -M ${MULTIREAD_MODEL} --report-counts

Then, we run HMM to reconstruct the genome::

    gbrs reconstruct -e emase.genes.tpm \
                     -x ${GBRS_DATA}/avecs.npz \
                     -g ${GBRS_DATA}/gene_ids.ordered.npz \
                     -t ${GBRS_DATA}/tranprob.DO.G1.F.npz \
                     -o ${OUTDIR}

We can now quantify allele-specific expressions on diploid transcriptome::

    gbrs quantify -i bowtie.8-way.transcriptome.h5 \
                  -g ${GBRS_DATA}/REF/emase.gene2transcripts.tsv \
                  -G gbrs.reconstructed.genotypes.tsv \
                  -L ${GBRS_DATA}/8-way/emase.pooled.transcripts.info \
                  -M ${MULTIREAD_MODEL} --report-counts


Genotype probabilities are on a grid of genes. For eQTL mapping or plotting genome reconstruction, we may want to interpolate probability on a grid at the genome scale.::

    gbrs interpolate -g marker_grid_64K.wYwMT.txt \
                     -i gbrs.reconstructed.genoprobs.npz \
                     -o gbrs.interpolated.genoprobs.npz


To plot a reconstructed genome::

    gbrs plot -i gbrs.gamma.on_grid.npz -o gbrs.genome.pdf -n SAMPLE_ID


To use GBRS in a project
~~~~~~~~~~~~~~~~~~~~~~~~

All the features are available as a python module. You can simply::

    import gbrs
