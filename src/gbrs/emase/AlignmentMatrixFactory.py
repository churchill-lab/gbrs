# standard library imports
import os
import struct

# 3rd party library imports
from scipy.sparse import coo_matrix
import numpy as np
import pysam
import tables

# local library imports
from gbrs import utils


logger = utils.get_logger('gbrs')


class AlignmentMatrixFactory:
    def __init__(self, alnfile):
        self.alnfile = alnfile
        self.hname = None
        self.lname = None
        self.rname = None
        self.tmpfiles = None

    def prepare(self, haplotypes: list[str], loci, delim: str = '_', outdir=None):
        if len(haplotypes) > 0:
            # Suffices given
            self.hname = haplotypes
        else:
            # Suffix not given
            self.hname = ['h0']
        self.lname = loci
        self.rname = set()

        save = pysam.set_verbosity(0)
        fh = pysam.AlignmentFile(self.alnfile, 'rb')
        pysam.set_verbosity(save)

        for aln in fh.fetch(until_eof=True):
            self.rname.add(aln.query_name)
        # loop over the file to get all read names.
        # The query_name attribute of the pysam.AlignedSegment object is the read name.
        # The read names are stored in a set to remove duplicates.

        self.rname = np.array(sorted(list(self.rname)))
        # sorts the elements of the NumPy array self.rname and updates the array with the sorted values.

        num_loci = len(self.lname)
        num_reads = len(self.rname)
        lid = dict(zip(self.lname, np.arange(num_loci)))
        rid = dict(zip(self.rname, np.arange(num_reads)))
        self.tmpfiles = dict.fromkeys(self.hname)

        if outdir is None:
            outdir = os.path.dirname(self.alnfile)

        fhout = dict.fromkeys(self.hname)
        for hap in self.hname:
            outfile = os.path.join(outdir, f'{hap}_{os.getpid()}.bin')
            self.tmpfiles[hap] = outfile
            fhout[hap] = open(outfile, 'wb')

        save = pysam.set_verbosity(0)
        fh = pysam.AlignmentFile(self.alnfile, 'rb')
        pysam.set_verbosity(save)

        if len(haplotypes) > 0:
            # Suffices given
            for aln in fh.fetch(until_eof=True):
                if aln.flag != 4 and aln.flag != 8:
                    locus, hap = fh.get_reference_name(aln.tid).split(delim)
                    fhout[hap].write(struct.pack('>I', rid[aln.qname]))
                    fhout[hap].write(struct.pack('>I', lid[locus]))
        else:
            # Suffix not given
            hap = self.hname[0]
            for aln in fh.fetch(until_eof=True):
                if aln.flag != 4 and aln.flag != 8:
                    locus = fh.get_reference_name(aln.tid)
                    fhout[hap].write(struct.pack('>I', rid[aln.qname]))
                    fhout[hap].write(struct.pack('>I', lid[locus]))
        for hap in self.hname:
            fhout[hap].close()
        # make temp files for each haplotype, the files contain the read name and locus name of each read in the alignment file.
        # the read name and locus name are written to the file in binary format.

    def produce(
        self,
        h5file: str,
        title='Alignments',
        index_dtype='uint32',
        data_dtype=float,
        complib='zlib',
        incidence_only=True,
    ):
        h5fh = tables.open_file(h5file, 'w', title=title)
        fil = tables.Filters(complevel=1, complib=complib)
        h5fh.set_node_attr(h5fh.root, 'incidence_only', incidence_only)
        h5fh.set_node_attr(h5fh.root, 'mtype', 'csc_matrix')
        h5fh.set_node_attr(
            h5fh.root,
            'shape',
            (len(self.lname), len(self.hname), len(self.rname)),
        )
        h5fh.set_node_attr(h5fh.root, 'hname', self.hname)
        h5fh.create_carray(
            h5fh.root,
            'lname',
            obj=self.lname,
            title='Locus Names',
            filters=fil,
        )
        h5fh.create_carray(
            h5fh.root, 'rname', obj=self.rname, title='Read Names', filters=fil
        )
        for hid in range(len(self.hname)):  # loop hap and read? Make the temp files above for each BAM independently and then loop through them here?
            hap = self.hname[hid]
            infile = self.tmpfiles[hap]
            dmat = np.fromfile(open(infile, 'rb'), dtype='>I')
            dmat = dmat.reshape((int(len(dmat) / 2), 2)).T
            if dmat.shape[0] > 2:
                dvec = dmat[2]
            else:
                dvec = np.ones(dmat.shape[1])
            spmat = coo_matrix(
                (dvec, dmat[:2]), shape=(len(self.rname), len(self.lname))
            )
            # spmat contains the read name and locus name of each read in the alignment file, in the form: ((read name, locus name), read count))
            # spmat is a sparse matrix with the read names as the row indices and the locus names as the column indices.
            # the data in the matrix is the number of reads that align to a specific locus.

            spmat = spmat.tocsc()
            hgroup = h5fh.create_group(
                h5fh.root,
                f'h{hid}',
                f'Sparse matrix components for Haplotype {hid}',
            )
            # add header
            h5fh.create_carray(
                hgroup,
                'indptr',
                obj=spmat.indptr.astype(index_dtype),
                filters=fil,
            )
            h5fh.create_carray(
                hgroup,
                'indices',
                obj=spmat.indices.astype(index_dtype),
                filters=fil,
            )
            if not incidence_only:
                h5fh.create_carray(
                    hgroup,
                    'data',
                    obj=spmat.data.astype(data_dtype),
                    filters=fil,
                )
            # apply sparse matrix indexing and indptr.
        h5fh.flush()
        h5fh.close()

    def cleanup(self):
        for tmpfile in self.tmpfiles.items():
            os.remove(tmpfile[1])
