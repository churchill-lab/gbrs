=====
Usage
=====

To use GBRS in command line
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once multi-way EMASE finish, we run HMM to reconstruct its genome.::

    gbrs reconstruct -e emase.genes.tpm -x ${GBRS_DATA}/avecs.npz \
                     -g ${GBRS_DATA}/ENSMUSG.ids.wYwMT.npz -t ${GBRS_DATA}/tranprob.DO.G1.F.npz \
                     -c 1.5 -o ${OUTDIR}

Genotype probabilities are on a grid of genes. For eQTL mapping or plotting genome reconstruction, we may want to interpolate probability on a grid in genome scale.::

    gbrs interpolate -g marker_grid_64K.wYwMT.txt -i gbrs.gamma.npz -o gbrs.gamma.on_grid.npz


To plot a reconstructed genome::

    gbrs plot -i gbrs.gamma.on_grid.npz -o gbrs.genome.pdf -n sample_id


To use GBRS in a project
~~~~~~~~~~~~~~~~~~~~~~~~

All the features are available as a python module. You can simply::

    import gbrs
