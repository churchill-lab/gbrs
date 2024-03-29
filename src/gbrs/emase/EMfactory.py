# standard library imports
import time

# 3rd party library imports
from scipy.sparse import eye, lil_matrix
import numpy as np

# local library imports
from gbrs import utils
from gbrs.emase.AlignmentPropertyMatrix import AlignmentPropertyMatrix as APM

logger = utils.get_logger('gbrs')


class EMfactory:
    """
    A class that coordinate Expectation-Maximization
    """

    def __init__(self, alignments):
        self.probability = alignments
        self.allelic_expression = None
        self.grp_conv_mat = None
        self.t2t_mat = None
        self.target_lengths = None

    def prepare(
        self,
        pseudocount: float = 0.0,
        lenfile: str = None,
        read_length: int = 100
    ) -> None:
        """
        Initializes the probability of read origin according to the alignment profile

        Args:
            pseudocount: Uniform prior for allele specificity estimation
            lenfile: lengths file
            read_length: default read length
        """
        if self.probability.num_groups > 0:
            self.grp_conv_mat = lil_matrix(
                (self.probability.num_loci, self.probability.num_groups)
            )
            for i in range(self.probability.num_groups):
                self.grp_conv_mat[self.probability.groups[i], i] = 1.0
            self.grp_conv_mat = self.grp_conv_mat.tocsc()
            self.t2t_mat = eye(
                self.probability.num_loci, self.probability.num_loci
            )
            self.t2t_mat = self.t2t_mat.tolil()
            for tid_list in self.probability.groups:
                for ii in range(len(tid_list)):
                    for jj in range(ii):
                        i = tid_list[ii]
                        j = tid_list[jj]
                        self.t2t_mat[i, j] = 1
                        self.t2t_mat[j, i] = 1
            self.t2t_mat = self.t2t_mat.tocsc()
        if lenfile is not None:
            hid = dict(
                zip(
                    self.probability.hname,
                    np.arange(len(self.probability.hname)),
                )
            )
            self.target_lengths = np.zeros(
                (self.probability.num_loci, self.probability.num_haplotypes)
            )
            if self.probability.num_haplotypes > 1:
                with open(lenfile) as fh:
                    for curline in fh:
                        item = curline.rstrip().split('\t')
                        locus, hap = item[0].split('_')
                        self.target_lengths[
                            self.probability.lid[locus], hid[hap]
                        ] = max(float(item[1]) - read_length + 1.0, 1.0)
            elif self.probability.num_haplotypes > 0:
                with open(lenfile) as fh:
                    for curline in fh:
                        item = curline.rstrip().split('\t')
                        self.target_lengths[
                            self.probability.lid[item[0]], 0
                        ] = max(float(item[1]) - read_length + 1.0, 1.0)
            else:
                raise RuntimeError(
                    'There is something wrong with your emase-format alignment file.'
                )
            self.target_lengths = self.target_lengths.transpose()
            # self.target_lengths = self.target_lengths.transpose() / read_length  # lengths in terms of read counts
            if not np.all(self.target_lengths > 0.0):
                raise RuntimeError(
                    'There exist transcripts missing length information.'
                )
        self.probability.normalize_reads(
            axis=APM.Axis.READ
        )  # Initialize alignment probability matrix
        self.allelic_expression = self.probability.sum(axis=APM.Axis.READ)
        if (
            self.target_lengths is not None
        ):  # allelic_expression will be at depth-level
            self.allelic_expression = np.divide(
                self.allelic_expression, self.target_lengths
            )
        if pseudocount > 0.0:  # pseudocount is at depth-level
            orig_allelic_expression_sum = self.allelic_expression.sum()
            nzloci = np.nonzero(self.allelic_expression)[1]
            self.allelic_expression[:, nzloci] += pseudocount
            self.allelic_expression *= (
                orig_allelic_expression_sum / self.allelic_expression.sum()
            )  # original depth scale

    def reset(self, pseudocount: float = 0.0) -> None:
        """
        Initializes the probability of read origin according to the alignment
        profile.

        Args:
            pseudocount: Uniform prior for allele specificity estimation
        """
        self.probability.reset()
        self.probability.normalize_reads(
            axis=APM.Axis.READ
        )  # Initialize alignment probability matrix
        self.allelic_expression = self.probability.sum(axis=APM.Axis.READ)
        if (
            self.target_lengths is not None
        ):  # allelic_expression will be at depth-level
            self.allelic_expression = np.divide(
                self.allelic_expression, self.target_lengths
            )
        if pseudocount > 0.0:  # pseudocount is at depth-level
            orig_allelic_expression_sum = self.allelic_expression.sum()
            nzloci = np.nonzero(self.allelic_expression)[1]
            self.allelic_expression[:, nzloci] += pseudocount
            self.allelic_expression *= (
                orig_allelic_expression_sum / self.allelic_expression.sum()
            )  # original depth scale

    def get_allelic_expression(self, at_group_level: bool = False):
        if at_group_level:
            return self.allelic_expression * self.grp_conv_mat
        else:
            return self.allelic_expression.copy()

    def update_probability_at_read_level(self, model: int = 3) -> None:
        """
        Updates the probability of read origin at read level

        Normalization model:
            1: Gene->Allele->Isoform
            2: Gene->Isoform->Allele
            3: Gene->Isoform*Allele
            4: Gene*Isoform*Allele

        Args:
            model: Normalization model
        """
        self.probability.reset()  # reset to alignment incidence matrix
        if model == 1:
            self.probability.multiply(
                self.allelic_expression, axis=APM.Axis.READ
            )
            self.probability.normalize_reads(
                axis=APM.Axis.HAPLOGROUP, grouping_mat=self.t2t_mat
            )
            haplogroup_sum_mat = self.allelic_expression * self.t2t_mat
            self.probability.multiply(haplogroup_sum_mat, axis=APM.Axis.READ)
            self.probability.normalize_reads(
                axis=APM.Axis.GROUP, grouping_mat=self.t2t_mat
            )
            self.probability.multiply(
                haplogroup_sum_mat.sum(axis=0), axis=APM.Axis.HAPLOTYPE
            )
            self.probability.normalize_reads(axis=APM.Axis.READ)
        elif model == 2:
            self.probability.multiply(
                self.allelic_expression, axis=APM.Axis.READ
            )
            self.probability.normalize_reads(axis=APM.Axis.LOCUS)
            self.probability.multiply(
                self.allelic_expression.sum(axis=0), axis=APM.Axis.HAPLOTYPE
            )
            self.probability.normalize_reads(
                axis=APM.Axis.GROUP, grouping_mat=self.t2t_mat
            )
            self.probability.multiply(
                (self.allelic_expression * self.t2t_mat).sum(axis=0),
                axis=APM.Axis.HAPLOTYPE,
            )
            self.probability.normalize_reads(axis=APM.Axis.READ)
        elif model == 3:
            self.probability.multiply(
                self.allelic_expression, axis=APM.Axis.READ
            )
            self.probability.normalize_reads(
                axis=APM.Axis.GROUP, grouping_mat=self.t2t_mat
            )
            self.probability.multiply(
                (self.allelic_expression * self.t2t_mat).sum(axis=0),
                axis=APM.Axis.HAPLOTYPE,
            )
            self.probability.normalize_reads(axis=APM.Axis.READ)
        elif model == 4:
            self.probability.multiply(
                self.allelic_expression, axis=APM.Axis.READ
            )
            self.probability.normalize_reads(axis=APM.Axis.READ)
        else:
            raise RuntimeError(
                'The read normalization model should be 1, 2, 3, or 4.'
            )

    def update_allelic_expression(self, model: int = 3) -> None:
        """
        A single EM step: Update probability at read level and then re-estimate allelic specific expression

        Normalization model:
            1: Gene->Allele->Isoform
            2: Gene->Isoform->Allele
            3: Gene->Isoform*Allele
            4: Gene*Isoform*Allele

        Args:
            model: Normalization model
        """
        self.update_probability_at_read_level(model)
        self.allelic_expression = self.probability.sum(axis=APM.Axis.READ)
        if self.target_lengths is not None:
            self.allelic_expression = np.divide(
                self.allelic_expression, self.target_lengths
            )

    def run(
        self,
        model: int,
        tol: float = 0.001,
        max_iters: int = 999,
        verbose: bool = True
    ) -> None:
        """
        Runs EM iterations

        Normalization model:
            1: Gene->Allele->Isoform
            2: Gene->Isoform->Allele
            3: Gene->Isoform*Allele
            4: Gene*Isoform*Allele

        Args:
            model: Normalization model
            tol: Tolerance for termination
            max_iters: Maximum number of iterations until termination
            verbose: Display information on how EM is running
        """
        orig_err_states = np.seterr(all='raise')
        np.seterr(under='ignore')
        if verbose:
            print('')
            print('Iter No  Time (hh:mm:ss)    Total change (TPM)  ')
            print('-------  ---------------  ----------------------')

        num_iters = 0
        err_sum = 1000000.0
        time0 = time.time()
        target_err = 1000000.0 * tol
        while err_sum > target_err and num_iters < max_iters:
            prev_isoform_expression = self.get_allelic_expression().sum(axis=0)
            prev_isoform_expression *= (
                1000000.0 / prev_isoform_expression.sum()
            )
            self.update_allelic_expression(model=model)
            curr_isoform_expression = self.get_allelic_expression().sum(axis=0)
            curr_isoform_expression *= (
                1000000.0 / curr_isoform_expression.sum()
            )
            err = np.abs(curr_isoform_expression - prev_isoform_expression)
            err_sum = err.sum()
            num_iters += 1
            if verbose:
                time1 = time.time()
                delmin, s = divmod(int(time1 - time0), 60)
                h, m = divmod(delmin, 60)
                print(
                    ' %5d      %4d:%02d:%02d     %9.1f / 1000000'
                    % (num_iters, h, m, s, err_sum)
                )

    def report_read_counts(
        self, filename, grp_wise=False, reorder='as-is', notes=None
    ):
        """
        Export read counts

        Args:
            filename: file name for output
            grp_wise: whether the report is at isoform level or gene level
            reorder: whether the report should be either 'decreasing' or
                'increasing' order or just 'as-is'
            notes: notes for the group
        """
        expected_read_counts = self.probability.sum(axis=APM.Axis.READ)
        if grp_wise:
            lname = self.probability.gname
            expected_read_counts = expected_read_counts * self.grp_conv_mat
        else:
            lname = self.probability.lname
        total_read_counts = expected_read_counts.sum(axis=0)
        if reorder == 'decreasing':
            report_order = np.argsort(total_read_counts.flatten())
            report_order = report_order[::-1]
        elif reorder == 'increasing':
            report_order = np.argsort(total_read_counts.flatten())
        elif reorder == 'as-is':
            # report in the original locus order
            report_order = np.arange(len(lname))
        cntdata = np.vstack((expected_read_counts, total_read_counts))
        fhout = open(filename, 'w')
        fhout.write('locus\t' + '\t'.join(self.probability.hname) + '\ttotal')
        if notes is not None:
            fhout.write('\tnotes')
        fhout.write('\n')
        for locus_id in report_order:
            lname_cur = lname[locus_id]
            lout = [lname_cur]
            lout.extend(list(map(str, cntdata[:, locus_id].ravel())))
            fhout.write('\t'.join(lout))
            if notes is not None:
                fhout.write('\t%s' % notes[lname_cur])
            fhout.write('\n')
        fhout.close()

    def report_depths(
        self, filename, tpm=True, grp_wise=False, reorder='as-is', notes=None
    ) -> None:
        """
        Exports expected depths

        Args:
            filename: file name for output
            tpm: True for tpms
            grp_wise: whether the report is at isoform level or gene level
            reorder: whether the report should be either 'decreasing' or
                'increasing' order or just 'as-is'
            notes: notes for the group
        """
        if grp_wise:
            lname = self.probability.gname
            depths = self.allelic_expression * self.grp_conv_mat
        else:
            lname = self.probability.lname
            depths = self.allelic_expression
        if tpm:
            depths *= 1000000.0 / depths.sum()
        total_depths = depths.sum(axis=0)
        if reorder == 'decreasing':
            report_order = np.argsort(total_depths.flatten())
            report_order = report_order[::-1]
        elif reorder == 'increasing':
            report_order = np.argsort(total_depths.flatten())
        elif reorder == 'as-is':
            # report in the original locus order
            report_order = np.arange(len(lname))
        cntdata = np.vstack((depths, total_depths))
        fhout = open(filename, 'w')
        fhout.write('locus\t' + '\t'.join(self.probability.hname) + '\ttotal')
        if notes is not None:
            fhout.write('\tnotes')
        fhout.write('\n')
        for locus_id in report_order:
            lname_cur = lname[locus_id]
            fhout.write(
                '\t'.join(
                    [lname_cur] + list(map(str, cntdata[:, locus_id].ravel()))
                )
            )
            if notes is not None:
                fhout.write(f'\t{notes[lname_cur]}')
            fhout.write('\n')
        fhout.close()

    def export_posterior_probability(
        self, filename: str, title: str = 'Posterior Probability'
    ) -> None:
        """
        Writes the posterior probability of read origin.

        Args:
            filename: File name for output
            title: the title of the posterior probability matrix
        """
        self.probability.save(h5file=filename, title=title)
