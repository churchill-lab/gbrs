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
    Expectation-Maximization (EM) algorithm implementation for allele-specific expression analysis.
    
    This class implements the EMASE (Expectation-Maximization for Allele-Specific Expression)
    algorithm, which estimates allele-specific expression levels from RNA-seq alignment data.
    The algorithm iteratively updates alignment probabilities and expression estimates until
    convergence.
    
    The EM algorithm consists of two main steps that are repeated:
    1. **Expectation (E-step)**: Update alignment probabilities given current expression estimates
    2. **Maximization (M-step)**: Update expression estimates given current alignment probabilities
    
    Key Features:
    - Supports multiple normalization models for handling multi-mapping reads
    - Handles transcript length normalization for proper expression quantification
    - Provides convergence monitoring and iteration control
    - Supports pseudocount addition for regularization
    - Generates comprehensive output reports
    
    This implementation is based on the methodology described in:
    - Munger et al. (2014) Genetics 198:59-73
    - Choi et al. (2020) bioRxiv 10.1101/2020.10.11.335323
    
    Attributes:
        probability (AlignmentPropertyMatrix): Alignment probability matrix
        allelic_expression (numpy.ndarray): Current allele-specific expression estimates
        grp_conv_mat (scipy.sparse.csc_matrix): Group conversion matrix for transcript isoforms
        t2t_mat (scipy.sparse.csc_matrix): Transcript-to-transcript incidence matrix
        target_lengths (numpy.ndarray): Effective lengths for transcript length normalization
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
        Initialize the EM algorithm with alignment data and optional parameters.
        
        This method sets up the EM algorithm by:
        1. Creating group conversion matrices for transcript isoforms
        2. Loading transcript length information for normalization
        3. Initializing allele-specific expression estimates
        4. Applying pseudocount regularization if specified
        
        Args:
            pseudocount (float, optional): Uniform prior added to all expression estimates
                for regularization. Defaults to 0.0 (no regularization).
            lenfile (str, optional): Path to transcript length file. The file should contain
                tab-separated values with transcript name and effective length. If provided,
                expression estimates are normalized by transcript length.
            read_length (int, optional): Default read length for length normalization.
                Used to calculate effective transcript length. Defaults to 100.
                
        Raises:
            RuntimeError: If transcript length file is invalid or missing transcripts
            
        Note:
            The lenfile format should be:
            ```
            Transcript1    1000
            Transcript2    1500
            Transcript3    800
            ```
            
            Effective length = transcript_length - read_length + 1
            
            Pseudocount regularization helps prevent zero expression estimates
            and improves algorithm stability, especially for low-coverage transcripts.
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
        Update alignment probabilities at the read level (E-step of EM algorithm).
        
        This method implements the Expectation step of the EM algorithm, updating
        the probability that each read originated from each possible transcript-haplotype
        combination given the current expression estimates.
        
        The method supports four different normalization models for handling
        multi-mapping reads and transcript isoforms:
        
        **Model 1: Gene->Allele->Isoform**
        - Normalize first across haplotypes, then across isoforms within each haplotype
        - Assumes hierarchical structure: gene → allele → isoform
        
        **Model 2: Gene->Isoform->Allele**  
        - Normalize first across isoforms, then across haplotypes within each isoform
        - Assumes hierarchical structure: gene → isoform → allele
        
        **Model 3: Gene->Isoform*Allele** (default)
        - Normalize across both isoforms and haplotypes simultaneously
        - Assumes independent effects of isoform and allele
        
        **Model 4: Gene*Isoform*Allele**
        - No hierarchical normalization, treats all combinations equally
        - Assumes complete independence between gene, isoform, and allele
        
        Args:
            model (int): Normalization model to use (1, 2, 3, or 4). Defaults to 3.
                
        Note:
            This method modifies the probability matrix in-place. The choice of
            normalization model can significantly affect results, especially for
            genes with multiple isoforms or complex allele-specific patterns.
            
            Model 3 is recommended for most applications as it provides a good
            balance between biological realism and computational efficiency.
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
        Run the complete EM algorithm until convergence.
        
        This method executes the full EM algorithm, iterating between the
        Expectation and Maximization steps until convergence or maximum
        iterations is reached.
        
        The algorithm monitors convergence by tracking changes in TPM
        (Transcripts Per Million) values between iterations. Convergence
        is achieved when the total absolute change in TPM values across
        all transcripts falls below the specified tolerance.
        
        Args:
            model (int): Normalization model to use (1, 2, 3, or 4):
                - 1: Gene->Allele->Isoform
                - 2: Gene->Isoform->Allele  
                - 3: Gene->Isoform*Allele (recommended)
                - 4: Gene*Isoform*Allele
            tol (float, optional): Convergence tolerance. Algorithm stops when
                total TPM change < tol * 1,000,000. Defaults to 0.001.
            max_iters (int, optional): Maximum number of EM iterations.
                Defaults to 999.
            verbose (bool, optional): If True, print progress information.
                Defaults to True.
                
        Note:
            The convergence criterion is based on TPM values rather than
            raw expression values to ensure scale-invariant convergence.
            
            Typical convergence requires 10-50 iterations depending on
            data complexity and tolerance settings.
            
            If the algorithm doesn't converge within max_iters, the final
            estimates may not be optimal. Consider increasing max_iters
            or adjusting the tolerance.
            
        Examples:
            # Run with default settings
            em_factory.run(model=3)
            
            # Run with strict convergence
            em_factory.run(model=3, tol=0.0001, max_iters=2000)
            
            # Run silently
            em_factory.run(model=3, verbose=False)
        """
        np.seterr(all='raise')
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
        self,
        filename,
        grp_wise=False,
        reorder='as-is',
        notes=None
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
        self,
        filename,
        tpm=True,
        grp_wise=False,
        reorder='as-is',
        notes=None
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
        self,
        filename: str,
        title: str = 'Posterior Probability'
    ) -> None:
        """
        Writes the posterior probability of read origin.

        Args:
            filename: File name for output
            title: the title of the posterior probability matrix
        """
        self.probability.save(h5_file=filename, title=title)
