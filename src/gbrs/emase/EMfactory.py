# standard library imports
import time

# 3rd party library imports
from scipy.sparse import eye, lil_matrix, csc_matrix
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
    convergence, providing robust quantification of transcript abundances in multiparent populations.
    
    The EM algorithm consists of two main steps that are repeated:
    **Expectation (E-step)**: Update alignment probabilities given current expression estimates
    **Maximization (M-step)**: Update expression estimates given current alignment probabilities
    
    Key Features:
    - Supports multiple normalization models for handling multi-mapping reads
    - Handles transcript length normalization for proper expression quantification
    - Provides convergence monitoring and iteration control
    - Supports pseudocount addition for regularization
    - Generates comprehensive output reports
    
    This implementation is based on the methodology described in:
    - Munger et al. (2014) Genetics 198:59-73
    - Choi et al. (2020) bioRxiv 10.1101/2020.10.11.335323
    
    The algorithm is designed for multiparent populations (e.g., Diversity Outbred mice)
    where RNA-seq reads are aligned to a hybrid transcriptome containing transcripts
    from all founder strains. It estimates the contribution of each founder haplotype
    to the observed expression patterns.
    
    Attributes:
        probability (AlignmentPropertyMatrix): Alignment probability matrix containing
            read-to-transcript-haplotype alignment probabilities
        allelic_expression (numpy.ndarray): Current allele-specific expression estimates with shape
            (num_haplotypes, num_loci)
        grp_conv_mat (scipy.sparse.csc_matrix): Group conversion matrix for aggregating
            transcript-level estimates to gene-level estimates
        t2t_mat (scipy.sparse.csc_matrix): Transcript-to-transcript incidence matrix for handling
            transcript isoforms within genes
        target_lengths (numpy.ndarray): Effective lengths for transcript length normalization with
            shape (num_haplotypes, num_loci)
    """

    def __init__(self, alignments: APM) -> None:
        """
        Initialize the EMfactory with alignment data.
        
        Args:
            alignments: AlignmentPropertyMatrix containing RNA-seq alignment data with
                read-to-transcript-haplotype relationships
        """
        #self.probability: APM = alignments
        #self.allelic_expression: np.ndarray | None = None
        #self.grp_conv_mat: csc_matrix | None = None
        #self.t2t_mat: csc_matrix | None = None
        #self.target_lengths: np.ndarray | None = None
        self.probability = alignments
        self.allelic_expression = None
        self.grp_conv_mat = None
        self.t2t_mat = None
        self.target_lengths = None


    def prepare(
        self,
        pseudocount: float = 0.0,
        length_file: str | None = None,
        read_length: int = 100
    ) -> None:
        """
        Initialize the EM algorithm with alignment data and optional parameters.
        
        This method sets up the EM algorithm by:
        - Creating group conversion matrices for transcript isoforms
        - Loading transcript length information for normalization
        - Initializing allele-specific expression estimates
        - Applying pseudocount regularization if specified
        
        The method handles the complex structure of multiparent populations by creating
        matrices that capture relationships between transcripts, haplotypes, and genes.
        This is essential for proper allele-specific expression quantification.
        
        Args:
            pseudocount: Uniform prior added to all expression estimates for regularization.
                Defaults to 0.0 (no regularization). Values > 0 help prevent zero expression
                estimates and improve algorithm stability, especially for low-coverage transcripts.
            length_file: Path to transcript length file. The file should contain tab-separated
                values with transcript name and effective length. If provided, expression estimates
                are normalized by transcript length to account for length bias in RNA-seq.
            read_length: Default read length for length normalization. Used to calculate effective
                transcript length as transcript_length - read_length + 1.
                
        Raises:
            RuntimeError: If transcript length file is invalid or missing transcripts
            
        Note:
            The length file format should be:
            ```
            Transcript1_A    1000
            Transcript1_B    1000
            Transcript2_A    1500
            Transcript2_B    1500
            ```
            
            For multiparent populations, transcript names should include haplotype suffixes
            (e.g., "_A", "_B", "_C", etc.) to distinguish between founder strain variants.
            
            Effective length = transcript_length - read_length + 1
            
            Pseudocount regularization helps prevent zero expression estimates and improves
            algorithm stability, especially for low-coverage transcripts. The pseudocount is
            added at the depth level and then scaled to maintain the original total expression.
        """
        if self.probability.num_groups > 0:
            # create group conversion matrix for transcript-to-gene aggregation
            self.grp_conv_mat = lil_matrix((self.probability.num_loci, self.probability.num_groups))

            for i in range(self.probability.num_groups):
                self.grp_conv_mat[self.probability.groups[i], i] = 1.0

            self.grp_conv_mat = self.grp_conv_mat.tocsc()
            
            # create transcript-to-transcript incidence matrix for isoform handling
            self.t2t_mat = eye(self.probability.num_loci, self.probability.num_loci)
            self.t2t_mat = self.t2t_mat.tolil()

            for tid_list in self.probability.groups:
                for ii in range(len(tid_list)):
                    for jj in range(ii):
                        i = tid_list[ii]
                        j = tid_list[jj]
                        self.t2t_mat[i, j] = 1
                        self.t2t_mat[j, i] = 1

            self.t2t_mat = self.t2t_mat.tocsc()
            
        if length_file is not None:
            # create mapping from haplotype names to indices
            hid = dict(zip(self.probability.hname, np.arange(len(self.probability.hname)),))
            self.target_lengths = np.zeros(
                (self.probability.num_loci, self.probability.num_haplotypes)
            )

            if self.probability.num_haplotypes > 1:
                # handle multiparent populations with haplotype-specific transcripts
                with open(length_file) as fh:
                    for line in fh:
                        item = line.rstrip().split('\t')
                        locus, hap = item[0].split('_')

                        self.target_lengths[
                            self.probability.lid[locus], hid[hap]
                        ] = max(float(item[1]) - read_length + 1.0, 1.0)
            elif self.probability.num_haplotypes > 0:
                # handle single haplotype case
                with open(length_file) as fh:
                    for line in fh:
                        item = line.rstrip().split('\t')
                        self.target_lengths[
                            self.probability.lid[item[0]], 0
                        ] = max(float(item[1]) - read_length + 1.0, 1.0)
            else:
                raise RuntimeError('There is an issue with your emase-format alignment file.')
            self.target_lengths = self.target_lengths.transpose()

            if not np.all(self.target_lengths > 0.0):
                raise RuntimeError('There exist transcripts missing length information.')

        # initialize alignment probability matrix
        self.probability.normalize_reads(axis=APM.Axis.READ)
        self.allelic_expression = self.probability.sum(axis=APM.Axis.READ)

        if self.target_lengths is not None:
            # normalize by transcript length to get depth-level expression
            self.allelic_expression = np.divide(self.allelic_expression, self.target_lengths)

        if pseudocount > 0.0:
            # apply pseudocount regularization at depth level
            orig_allelic_expression_sum = self.allelic_expression.sum()
            nzloci = np.nonzero(self.allelic_expression)[1]
            self.allelic_expression[:, nzloci] += pseudocount
            # scale back to original depth scale
            self.allelic_expression *= (orig_allelic_expression_sum / self.allelic_expression.sum())


    def reset(self, pseudocount: float = 0.0) -> None:
        """
        Reset the EM algorithm to initial state with current alignment profile. This re-initializes
        the probability of read origin according to the current alignment profile and optionally
        applies pseudocount regularization. It's useful for restarting the EM algorithm with
        different parameters or after modifying the alignment data.
        
        Args:
            pseudocount: Uniform prior for allele specificity estimation. Adds a small prior to
                prevent zero probabilities in the EM algorithm.
        """
        self.probability.reset()

        # initialize alignment probability matrix
        self.probability.normalize_reads(axis=APM.Axis.READ)
        self.allelic_expression = self.probability.sum(axis=APM.Axis.READ)

        if self.target_lengths is not None:
            # normalize by transcript length to get depth-level expression
            self.allelic_expression = np.divide(self.allelic_expression, self.target_lengths)

        if pseudocount > 0.0:
            # apply pseudocount regularization at depth level
            orig_allelic_expression_sum = self.allelic_expression.sum()
            nzloci = np.nonzero(self.allelic_expression)[1]
            self.allelic_expression[:, nzloci] += pseudocount
            # Scale back to original depth scale
            self.allelic_expression *= (orig_allelic_expression_sum / self.allelic_expression.sum())


    def get_allelic_expression(self, at_group_level: bool = False) -> np.ndarray:
        """
        Get the current allele-specific expression estimates. Returns the current allele-specific
        expression estimates, optionally aggregated to the gene level if transcript groups are
        defined.
        
        Args:
            at_group_level: True to return the expression at the gene level, False at transcript
                level. If True, transcript-level estimates are aggregated using the group conversion
                matrix.
                
        Returns:
            Allele-specific expression matrix with shape (num_haplotypes, num_loci) or
            (num_haplotypes, num_groups) if at_group_level=True.
            
        Note:
            The returned matrix contains expression estimates in depth units (reads per effective
            length) if length normalization was applied, or raw read counts otherwise.
        """
        if at_group_level:
            return self.allelic_expression * self.grp_conv_mat
        else:
            return self.allelic_expression.copy()


    def update_probability_at_read_level(self, model: int = 3) -> None:
        """
        Update alignment probabilities at the read level (E-step of EM algorithm).
        
        This method implements the Expectation step of the EM algorithm, updating the probability
        that each read originated from each possible transcript-haplotype combination given the
        current expression estimates.
        
        The method supports four different normalization models for handling multi-mapping reads
        and transcript isoforms, each representing different assumptions about the hierarchical
        structure of gene expression:
        
        Model 1: Gene->Allele->Isoform
        - Normalize first across haplotypes, then across isoforms within each haplotype
        - Assumes hierarchical structure: gene → allele → isoform
        - Best for cases where allele-specific effects dominate isoform-specific effects
        
        Model 2: Gene->Isoform->Allele
        - Normalize first across isoforms, then across haplotypes within each isoform
        - Assumes hierarchical structure: gene → isoform → allele
        - Best for cases where isoform-specific effects dominate allele-specific effects
        
        Model 3: Gene->Isoform*Allele (default)
        - Normalize across both isoforms and haplotypes simultaneously
        - Assumes independent effects of isoform and allele
        - Recommended for most applications as it provides a good balance between
          biological realism and computational efficiency
        
        Model 4: Gene*Isoform*Allele
        - No hierarchical normalization, treats all combinations equally
        - Assumes complete independence between gene, isoform, and allele
        - Most computationally efficient but may oversimplify biological relationships
        
        Args:
            model: Normalization model to use (1, 2, 3, or 4). Defaults to 3.
                
        Note:
            This method modifies the probability matrix in-place. The choice of normalization model
            can significantly affect results, especially for genes with multiple isoforms or complex
            allele-specific patterns.
            
            Model 3 is recommended for most applications as it provides a good balance between
            biological realism and computational efficiency.
            
            The normalization process ensures that probabilities sum to 1 across all possible
            origins for each read, maintaining the probabilistic interpretation of the model.
        """
        # reset to alignment incidence matrix
        self.probability.reset()

        if model == 1:
            # Gene->Allele->Isoform normalization
            self.probability.multiply(self.allelic_expression, axis=APM.Axis.READ)
            self.probability.normalize_reads(axis=APM.Axis.HAPLOGROUP, grouping_mat=self.t2t_mat)
            haplogroup_sum_mat = self.allelic_expression * self.t2t_mat
            self.probability.multiply(haplogroup_sum_mat, axis=APM.Axis.READ)
            self.probability.normalize_reads(axis=APM.Axis.GROUP, grouping_mat=self.t2t_mat)
            self.probability.multiply(haplogroup_sum_mat.sum(axis=0), axis=APM.Axis.HAPLOTYPE)
            self.probability.normalize_reads(axis=APM.Axis.READ)
        elif model == 2:
            # Gene->Isoform->Allele normalization
            self.probability.multiply(self.allelic_expression, axis=APM.Axis.READ)
            self.probability.normalize_reads(axis=APM.Axis.LOCUS)
            self.probability.multiply(self.allelic_expression.sum(axis=0), axis=APM.Axis.HAPLOTYPE)
            self.probability.normalize_reads(axis=APM.Axis.GROUP, grouping_mat=self.t2t_mat)
            self.probability.multiply(
                (self.allelic_expression * self.t2t_mat).sum(axis=0),
                axis=APM.Axis.HAPLOTYPE,
            )
            self.probability.normalize_reads(axis=APM.Axis.READ)
        elif model == 3:
            # Gene->Isoform*Allele normalization (recommended)
            self.probability.multiply(self.allelic_expression, axis=APM.Axis.READ)
            self.probability.normalize_reads(axis=APM.Axis.GROUP, grouping_mat=self.t2t_mat)
            self.probability.multiply(
                (self.allelic_expression * self.t2t_mat).sum(axis=0),
                axis=APM.Axis.HAPLOTYPE,
            )
            self.probability.normalize_reads(axis=APM.Axis.READ)
        elif model == 4:
            # Gene*Isoform*Allele normalization (simplest)
            self.probability.multiply(self.allelic_expression, axis=APM.Axis.READ)
            self.probability.normalize_reads(axis=APM.Axis.READ)
        else:
            raise RuntimeError('The read normalization model should be 1, 2, 3, or 4.')


    def update_allelic_expression(self, model: int = 3) -> None:
        """
        Perform a single EM step: Update probability at read level and then re-estimate allelic
        expression.  This implements one complete iteration of the EM algorithm, combining the
        Expectation step (update_probability_at_read_level) and Maximization step (re-estimation
        of expression levels) in a single call.
        
        The method updates the allele-specific expression estimates based on the current alignment
        probabilities, which are then used in the next iteration to refine the probability
        estimates.
        
        Args:
            model: Normalization model to use (1, 2, 3, or 4):
                - 1: Gene->Allele->Isoform
                - 2: Gene->Isoform->Allele
                - 3: Gene->Isoform*Allele (recommended)
                - 4: Gene*Isoform*Allele
                
        Note:
            This method modifies both the probability matrix and allelic_expression in-place. It
            represents one iteration of the EM algorithm and should be called repeatedly until
            convergence is achieved.
        """
        self.update_probability_at_read_level(model)
        self.allelic_expression = self.probability.sum(axis=APM.Axis.READ)
        if self.target_lengths is not None:
            # normalize by transcript length to maintain depth-level expression
            self.allelic_expression = np.divide(self.allelic_expression, self.target_lengths)


    def run(
        self,
        model: int,
        tol: float = 0.001,
        max_iters: int = 999,
        verbose: bool = True
    ) -> None:
        """
        Run the complete EM algorithm until convergence.
        
        This method executes the full EM algorithm, iterating between the Expectation and
        Maximization steps until convergence or maximum iterations is reached.
        
        The algorithm monitors convergence by tracking changes in TPM (Transcripts Per Million)
        values between iterations. Convergence is achieved when the total absolute change in TPM
        values across all transcripts falls below the specified tolerance.
        
        The convergence criterion is based on TPM values rather than raw expression values to
        ensure scale-invariant convergence, making the algorithm robust to differences in
        sequencing depth between samples.
        
        Args:
            model: Normalization model to use (1, 2, 3, or 4):
                - 1: Gene->Allele->Isoform
                - 2: Gene->Isoform->Allele  
                - 3: Gene->Isoform*Allele (recommended)
                - 4: Gene*Isoform*Allele
            tol: Convergence tolerance. Algorithm stops when total TPM change < tol * 1,000,000.
                Default: 0.001 (1 TPM unit change per million transcripts)
            max_iters: Maximum number of EM iterations.
            verbose: If True, print progress information including iteration number, elapsed time,
                and convergence metric.
                
        Note:
            The convergence criterion is based on TPM values rather than raw expression values to
            ensure scale-invariant convergence.
            
            Typical convergence requires 10-50 iterations depending on data complexity and tolerance
            settings. Complex datasets with many multi-mapping reads may require more iterations.
            
            If the algorithm doesn't converge within max_iters, the final estimates may not be
            optimal. Consider increasing max_iters or adjusting the tolerance.
            
            The algorithm uses numpy's error handling to catch numerical issues and provides
            detailed progress information when verbose=True.
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
            # calculate TPM values for current iteration
            prev_isoform_expression = self.get_allelic_expression().sum(axis=0)
            prev_isoform_expression *= (1000000.0 / prev_isoform_expression.sum())
            
            # perform one EM step
            self.update_allelic_expression(model=model)
            
            # calculate TPM values after update
            curr_isoform_expression = self.get_allelic_expression().sum(axis=0)
            curr_isoform_expression *= (1000000.0 / curr_isoform_expression.sum())
            
            # calculate convergence metric
            err = np.abs(curr_isoform_expression - prev_isoform_expression)
            err_sum = err.sum()
            num_iters += 1

            if verbose:
                time1 = time.time()
                delmin, s = divmod(int(time1 - time0), 60)
                h, m = divmod(delmin, 60)
                print(' %5d      %4d:%02d:%02d     %9.1f / 1000000' % (num_iters, h, m, s, err_sum))


    def report_read_counts(
        self,
        filename: str,
        grp_wise: bool = False,
        reorder: str = 'as-is',
        notes: dict[str, str] | None = None
    ) -> None:
        """
        Export expected read counts to a tab-separated file.
        
        This method generates a comprehensive report of expected read counts for each
        transcript/haplotype combination. The expected read counts represent the posterior
        expectations of how many reads originated from each transcript given the final expression
        estimates.
        
        The output file contains expected read counts for each haplotype, plus a total column, and
        optionally includes notes for each transcript/gene.
        
        Args:
            filename: Path for the output file (TSV format).
            grp_wise: Whether to report at gene level (True) or transcript level (False).
                If True, transcript-level counts are aggregated to gene level.
            reorder: Sorting order for the output:
                - 'decreasing': Sort by total counts in descending order
                - 'increasing': Sort by total counts in ascending order
                - 'as-is': Maintain original order
            notes: Optional dictionary mapping transcript/gene names to notes.
                These notes are included as an additional column in the output.
                
        Output Format:
            The output file is tab-separated with the following columns:
            - locus: Transcript or gene identifier
            - haplotype columns: Expected read counts for each haplotype
            - total: Sum of read counts across all haplotypes
            - notes: Optional notes column (if notes parameter provided)
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

        count_data = np.vstack((expected_read_counts, total_read_counts))
        fh_out = open(filename, 'w')

        fh_out.write('locus\t' + '\t'.join(self.probability.hname) + '\ttotal')

        if notes is not None:
            fh_out.write('\tnotes')

        fh_out.write('\n')

        for locus_id in report_order:
            lname_cur = lname[locus_id]
            lout = [lname_cur]
            lout.extend(list(map(str, count_data[:, locus_id].ravel())))
            fh_out.write('\t'.join(lout))

            if notes is not None:
                fh_out.write('\t%s' % notes[lname_cur])

            fh_out.write('\n')
        fh_out.close()


    def report_depths(
        self,
        filename: str,
        tpm: bool = True,
        grp_wise: bool = False,
        reorder: str = 'as-is',
        notes: dict[str, str] | None = None
    ) -> None:
        """
        Export expression depth estimates to a tab-separated file.
        
        This method generates a comprehensive report of expression depth estimates for each
        transcript/haplotype combination. The depth estimates represent the final expression levels
        from the EM algorithm, optionally normalized to TPM (Transcripts Per Million) units.
        
        The output file contains depth estimates for each haplotype, plus a total column, and
        optionally includes notes for each transcript/gene.
        
        Args:
            filename: Path for the output file (TSV format).
            tpm: Whether to normalize to TPM units (True) or report raw depth values (False).
                TPM normalization accounts for transcript length and sequencing depth.
            grp_wise: Whether to report at gene level (True) or transcript level (False).
                If True, transcript-level estimates are aggregated to gene level.
            reorder: Sorting order for the output:
                - 'decreasing': Sort by total depth in descending order
                - 'increasing': Sort by total depth in ascending order
                - 'as-is': Maintain original order
            notes: Optional dictionary mapping transcript/gene names to notes.
                These notes are included as an additional column in the output.
                
        Output Format:
            The output file is tab-separated with the following columns:
            - locus: Transcript or gene identifier
            - haplotype columns: Expression depth for each haplotype
            - total: Sum of expression depth across all haplotypes
            - notes: Optional notes column (if notes parameter provided)
            
        Note:
            When tpm=True, the depth values are normalized to Transcripts Per Million units,
            which accounts for transcript length and sequencing depth, making them comparable
            across samples and transcripts of different lengths.
            
            When tpm=False, the depth values represent raw expression estimates in units
            of reads per effective length (if length normalization was applied) or raw
            read counts (otherwise).
        """
        if grp_wise:
            lname = self.probability.gname
            depths = self.allelic_expression * self.grp_conv_mat
        else:
            lname = self.probability.lname
            depths = self.allelic_expression

        if tpm:
            # normalize to TPM units
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

        count_data = np.vstack((depths, total_depths))
        fh_out = open(filename, 'w')

        fh_out.write('locus\t' + '\t'.join(self.probability.hname) + '\ttotal')

        if notes is not None:
            fh_out.write('\tnotes')

        fh_out.write('\n')

        for locus_id in report_order:
            lname_cur = lname[locus_id]
            fh_out.write(
                '\t'.join(
                    [lname_cur] + list(map(str, count_data[:, locus_id].ravel()))
                )
            )

            if notes is not None:
                fh_out.write(f'\t{notes[lname_cur]}')

            fh_out.write('\n')
        fh_out.close()


    def export_posterior_probability(
        self,
        filename: str,
        title: str = 'Posterior Probability'
    ) -> None:
        """
        Export the posterior probability matrix to an HDF5 file.
        
        This method saves the final posterior probability matrix to an HDF5 file, which contains
        the probability that each read originated from each possible transcript-haplotype
        combination given the final expression estimates.
        
        The posterior probability matrix is a key output of the EM algorithm and provides detailed
        information about the uncertainty in read assignments. This information can be used for
        downstream analyses that require uncertainty quantification or for debugging the EM algorithm.
        
        Args:
            filename: Path for the output HDF5 file.
            title: Title/description for the probability matrix in the HDF5 file.
                
        Note:
            The posterior probability matrix has the same structure as the input
            alignment matrix but contains probability values (summing to 1 across
            all possible origins for each read) rather than binary alignment indicators.
            
            This file can be quite large for datasets with many reads and transcripts,
            as it stores the full probability matrix in sparse format.
        """
        self.probability.save(h5_file=filename, title=title)
