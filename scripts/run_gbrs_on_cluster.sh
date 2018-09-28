#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=32gb,walltime=23:59:00

module load bowtie
module load samtools
module load Anaconda
source activate gbrs

cd $PBS_O_WORKDIR

if [ ! -d ${OUTDIR} ]; then
    mkdir -p ${OUTDIR}
fi

mkfifo ${OUTDIR}/bowtie.fifo
zcat $(ls -t ${INDIR}/*.fastq.gz) > ${OUTDIR}/bowtie.fifo &
bowtie -q -p 8 -a --best --strata --sam -v 3 -m 100 ${INDEXBASE} ${OUTDIR}/bowtie.fifo | samtools view -bS - > ${OUTDIR}/bowtie.bam
rm ${OUTDIR}/bowtie.fifo

gbrs bam2emase -i ${OUTDIR}/bowtie.bam \
               -m ${GBRS_DATA}/ref.transcripts.info \
               -s A,B,C,D,E,F,G,H \
               -o ${OUTDIR}/bowtie.h5

gbrs compress  -i ${OUTDIR}/bowtie.h5 -o ${OUTDIR}/bowtie.compressed.h5
rm ${OUTDIR}/bowtie.bam
rm ${OUTDIR}/bowtie.h5
mv ${OUTDIR}/bowtie.compressed.h5 ${OUTDIR}/bowtie.8-way.transcripts.h5

gbrs quantify -i ${OUTDIR}/bowtie.8-way.transcripts.h5 \
              -g ${GBRS_DATA}/ref.gene2transcripts.tsv \
              -L ${GBRS_DATA}/gbrs.hybridized.targets.info \
              -M 4 --report-alignment-counts \
              -o ${OUTDIR}/gbrs.quantified

gbrs reconstruct -e ${OUTDIR}/gbrs.quantified.multiway.genes.tpm \
                 -x ${GBRS_DATA}/avecs.npz \
                 -t ${GBRS_DATA}/tranprob.DO.G${GEN}.${SEX}.npz \
                 -g ${GBRS_DATA}/ref.gene_pos.ordered.npz \
                 -o ${OUTDIR}/gbrs.reconstructed \
                 -s ${SIGMA}

gbrs quantify -i ${OUTDIR}/bowtie.8-way.transcripts.h5 \
              -G ${OUTDIR}/gbrs.reconstructed.genotypes.tsv \
              -g ${GBRS_DATA}/ref.gene2transcripts.tsv \
              -L ${GBRS_DATA}/gbrs.hybridized.targets.info \
              -M 4 --report-alignment-counts \
              -o ${OUTDIR}/gbrs.quantified

gbrs interpolate -i ${OUTDIR}/gbrs.reconstructed.genoprobs.npz \
                 -g ${GBRS_DATA}/ref.genome_grid.txt \
                 -p ${GBRS_DATA}/ref.gene_pos.ordered.npz \
                 -o ${OUTDIR}/gbrs.interpolated.genoprobs.npz

export-genoprob-file -i ${OUTDIR}/gbrs.interpolated.genoprobs.npz \
                      -s A,B,C,D,E,F,G,H \
                      -g ${GBRS_DATA}/ref.genome_grid.txt

gbrs plot -i ${OUTDIR}/gbrs.interpolated.genoprobs.npz \
          -o ${OUTDIR}/gbrs.plotted.genome.pdf \
          -n ${SAMPLE_ID}

source deactivate
