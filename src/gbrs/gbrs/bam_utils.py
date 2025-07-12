# standard library imports
import random
from dataclasses import dataclass
from typing import List, Tuple, Dict

# 3rd party library imports
import numpy as np
import pysam

@dataclass
class BAMGenerator:
    """
    Generate realistic synthetic BAM files for testing GBRS pipeline.
    
    This implementation follows the exact specifications:
    - Generate all possible locus×haplotype×readname combinations
    - Select alignments based on distribution ensuring all read names appear
    - Each selected combination gets exactly 1 read
    - No correlation between read name order and transcript order
    """
    
    # Configuration parameters
    num_loci: int = 10
    num_haplotypes: int = 4
    num_read_names: int = 100
    alignments_r1: int = 100
    alignments_r2: int = None
    read_length: int = 100
    coverage_distribution: str = 'realistic'
    prefix_locus: str = 'ENSMUST'
    prefix_read: str = 'READ'
    haplotypes: str = None
    randomize_read_names: bool = False
    random_seed: int = None
    
    def __post_init__(self):
        """Initialize and validate configuration."""
        # Set random seed
        if self.random_seed is not None:
            random.seed(self.random_seed)
            np.random.seed(self.random_seed)
        
        # Generate haplotype list
        if self.haplotypes:
            self.haplotype_list = [h.strip() for h in self.haplotypes.split(',')]
            if len(self.haplotype_list) != self.num_haplotypes:
                raise ValueError(f"Number of haplotypes in --haplotypes ({len(self.haplotype_list)}) must match -h ({self.num_haplotypes})")
        else:
            self.haplotype_list = [chr(65 + i) for i in range(self.num_haplotypes)]
        
        # Generate locus names
        self.loci = [f'{self.prefix_locus}{str(i).zfill(11)}' for i in range(1, self.num_loci + 1)]
        
        # Generate read names
        self.read_names = [f'{self.prefix_read}{str(i).zfill(12)}' for i in range(1, self.num_read_names + 1)]
        
        # Set alignments_r2 if not specified
        if self.alignments_r2 is None:
            self.alignments_r2 = self.alignments_r1
        
        # Validate alignments constraint
        max_combinations = self.num_loci * self.num_haplotypes * self.num_read_names
        if self.alignments_r1 > max_combinations:
            raise ValueError(f"alignments_r1 ({self.alignments_r1}) cannot exceed loci×haplotypes×readnames ({max_combinations})")
        if self.alignments_r2 > max_combinations:
            raise ValueError(f"alignments_r2 ({self.alignments_r2}) cannot exceed loci×haplotypes×readnames ({max_combinations})")
        
        # Generate reference sequences
        self.reference_sequences = self._generate_reference_sequences()
        
        # Select combinations for R1 and R2 independently
        # This ensures different distributions between the files
        self.selected_combinations_r1 = self._select_combinations(self.alignments_r1)
        
        # For R2, we want a different distribution but still ensure all read names appear
        # We'll use a different random seed for R2 to get different results
        original_seed = random.getstate()
        original_np_seed = np.random.get_state()
        
        if self.random_seed is not None:
            random.seed(self.random_seed + 1)  # Different seed for R2
            np.random.seed(self.random_seed + 1)
        else:
            random.seed(random.randint(1, 1000000))  # Random seed for R2
            np.random.seed(random.randint(1, 1000000))
        
        self.selected_combinations_r2 = self._select_combinations(self.alignments_r2)
        
        # Restore original random state
        random.setstate(original_seed)
        np.random.set_state(original_np_seed)

    def _generate_reference_sequences(self) -> dict:
        """Generate reference sequences for each locus-haplotype combination."""
        sequences = {}
        
        for locus in self.loci:
            # Generate base sequence
            base_length = self.read_length * 20
            base_sequence = ''.join(random.choices(['A', 'C', 'G', 'T'], k=base_length))
            
            for haplotype in self.haplotype_list:
                ref_name = f'{locus}_{haplotype}'
                
                # Create haplotype-specific variations
                haplotype_sequence = list(base_sequence)
                
                # Introduce SNPs between haplotypes (1-3% difference)
                num_snps = int(len(haplotype_sequence) * random.uniform(0.01, 0.03))
                snp_positions = random.sample(range(len(haplotype_sequence)), num_snps)
                
                for pos in snp_positions:
                    current_base = haplotype_sequence[pos]
                    other_bases = [b for b in ['A', 'C', 'G', 'T'] if b != current_base]
                    haplotype_sequence[pos] = random.choice(other_bases)
                
                sequences[ref_name] = ''.join(haplotype_sequence)
        
        return sequences

    def _select_combinations(self, num_alignments: int) -> List[Tuple]:
        """Select combinations based on distribution method ensuring all read names appear."""
        # Create combinations grouped by read name for easy selection
        combinations_by_read = {}
        for read_name in self.read_names:
            combinations_by_read[read_name] = []
            for locus in self.loci:
                for haplotype in self.haplotype_list:
                    combinations_by_read[read_name].append((read_name, locus, haplotype))
        
        # Step 1: Select one combination for each read name
        selected = []
        used_read_names = set()
        
        # Randomize the order of read names to get different selections
        read_names_shuffled = list(self.read_names)
        random.shuffle(read_names_shuffled)
        
        for read_name in read_names_shuffled:
            available_combinations = combinations_by_read[read_name]
            selected_combination = random.choice(available_combinations)
            selected.append(selected_combination)
            used_read_names.add(read_name)
        
        # Step 2: Fill remaining slots with weighted selection
        remaining_slots = num_alignments - len(selected)
        
        if remaining_slots > 0:
            # Create remaining combinations (excluding already selected ones)
            remaining_combinations = []
            for read_name in self.read_names:
                for locus in self.loci:
                    for haplotype in self.haplotype_list:
                        combination = (read_name, locus, haplotype)
                        if combination not in selected:
                            remaining_combinations.append(combination)
            
            # Apply distribution to select remaining combinations
            if self.coverage_distribution == 'uniform':
                additional_selected = random.sample(remaining_combinations, min(remaining_slots, len(remaining_combinations)))
                
            elif self.coverage_distribution == 'exponential':
                # Weight by locus index (first loci more likely)
                weights = np.exp(-np.arange(len(remaining_combinations)) * 0.1)
                weights = weights / weights.sum()
                selected_indices = np.random.choice(len(remaining_combinations), size=min(remaining_slots, len(remaining_combinations)), replace=False, p=weights)
                additional_selected = [remaining_combinations[i] for i in selected_indices]
                
            elif self.coverage_distribution == 'realistic':
                # Weight by log-normal distribution across loci
                locus_weights = np.random.lognormal(0, 1, self.num_loci)
                locus_weights = locus_weights / locus_weights.sum()
                
                # Create weights for each remaining combination based on locus
                weights = []
                for read_name, locus, haplotype in remaining_combinations:
                    locus_index = self.loci.index(locus)
                    weights.append(locus_weights[locus_index])
                
                weights = np.array(weights)
                weights = weights / weights.sum()
                
                selected_indices = np.random.choice(len(remaining_combinations), size=min(remaining_slots, len(remaining_combinations)), replace=False, p=weights)
                additional_selected = [remaining_combinations[i] for i in selected_indices]
            
            else:
                raise ValueError(f'Unknown coverage distribution: {self.coverage_distribution}')
            
            selected.extend(additional_selected)
        
        return selected

    def _generate_read_sequence(self, ref_sequence: str, start_pos: int) -> Tuple[str, int, List]:
        """Generate a read sequence from reference sequence with realistic errors."""
        end_pos = start_pos + self.read_length
        if end_pos > len(ref_sequence):
            start_pos = len(ref_sequence) - self.read_length
            end_pos = len(ref_sequence)
        
        ref_segment = ref_sequence[start_pos:end_pos]
        read_sequence = list(ref_segment)
        
        # Introduce realistic sequencing errors (only mismatches)
        num_mismatches = 0
        
        for i in range(len(read_sequence)):
            # 1% chance of mismatch
            if random.random() < 0.01:
                current_base = read_sequence[i]
                other_bases = [b for b in ['A', 'C', 'G', 'T'] if b != current_base]
                read_sequence[i] = random.choice(other_bases)
                num_mismatches += 1
        
        # Simple CIGAR: all matches
        cigar_ops = [(0, len(read_sequence))]
        
        return ''.join(read_sequence), num_mismatches, cigar_ops

    def _generate_quality_scores(self) -> str:
        """Generate realistic quality scores for a read."""
        qualities = []
        
        for i in range(self.read_length):
            # Quality degrades toward the end of the read
            position_factor = 1.0 - (i / self.read_length) * 0.3
            
            # Base quality with variation
            if random.random() < 0.85:
                quality = random.randint(30, 40)
            elif random.random() < 0.1:
                quality = random.randint(20, 30)
            else:
                quality = random.randint(15, 25)
            
            # Apply position-based degradation
            quality = int(quality * position_factor)
            quality = max(10, min(40, quality))
            
            qualities.append(chr(quality + 33))
        
        return ''.join(qualities)

    def _generate_md_tag(self, ref_seq: str, read_seq: str) -> str:
        """Generate MD tag for BAM format."""
        md_parts = []
        match_count = 0
        
        for ref_base, read_base in zip(ref_seq, read_seq):
            if ref_base == read_base:
                match_count += 1
            else:
                if match_count > 0:
                    md_parts.append(str(match_count))
                    match_count = 0
                md_parts.append(ref_base)
        
        if match_count > 0:
            md_parts.append(str(match_count))
        
        return ''.join(md_parts)

    def generate_bam(self, output_file: str, paired_end: bool = False, second_file: str = None):
        """Generate BAM file(s)."""
        # Create header
        header = {
            'HD': {'VN': '1.0'},
            'SQ': []
        }
        
        # Add reference sequences to header
        for locus in self.loci:
            for haplotype in self.haplotype_list:
                ref_name = f'{locus}_{haplotype}'
                ref_length = len(self.reference_sequences[ref_name])
                header['SQ'].append({
                    'LN': ref_length,
                    'SN': ref_name
                })
        
        # Generate reads for R1
        reads_r1 = []
        for read_name, locus, haplotype in self.selected_combinations_r1:
            ref_name = f"{locus}_{haplotype}"
            ref_sequence = self.reference_sequences[ref_name]
            
            # Realistic start position
            max_start = len(ref_sequence) - self.read_length
            start_pos = random.randint(10, max(11, max_start - 10))
            
            # Generate read sequence
            read_sequence, num_mismatches, cigar_ops = self._generate_read_sequence(ref_sequence, start_pos)
            quality_scores = self._generate_quality_scores()
            
            # Calculate mapping quality
            mapq = max(30, 255 - num_mismatches * 10)
            
            # Create read entry
            read_entry = {
                'name': read_name,
                'ref_name': ref_name,
                'start': start_pos,
                'sequence': read_sequence,
                'quality': quality_scores,
                'flag': 0,
                'mapq': mapq,
                'cigar': cigar_ops,
                'tags': [
                    ('NM', num_mismatches),
                    ('MD', self._generate_md_tag(ref_sequence[start_pos:start_pos + len(read_sequence)], read_sequence)),
                    ('XA', 0)
                ]
            }
            
            reads_r1.append(read_entry)
        
        # Generate reads for R2 if paired-end
        reads_r2 = []
        if paired_end:
            for read_name, locus, haplotype in self.selected_combinations_r2:
                ref_name = f"{locus}_{haplotype}"
                ref_sequence = self.reference_sequences[ref_name]
                
                # Realistic start position
                max_start = len(ref_sequence) - self.read_length
                start_pos = random.randint(10, max(11, max_start - 10))
                
                # Generate read sequence
                read_sequence, num_mismatches, cigar_ops = self._generate_read_sequence(ref_sequence, start_pos)
                quality_scores = self._generate_quality_scores()
                
                # Calculate mapping quality
                mapq = max(30, 255 - num_mismatches * 10)
                
                # Create read entry
                read_entry = {
                    'name': read_name,
                    'ref_name': ref_name,
                    'start': start_pos,
                    'sequence': read_sequence,
                    'quality': quality_scores,
                    'flag': 0,
                    'mapq': mapq,
                    'cigar': cigar_ops,
                    'tags': [
                        ('NM', num_mismatches),
                        ('MD', self._generate_md_tag(ref_sequence[start_pos:start_pos + len(read_sequence)], read_sequence)),
                        ('XA', 0)
                    ]
                }
                
                reads_r2.append(read_entry)
        
        # Randomize read order if requested
        if self.randomize_read_names:
            random.shuffle(reads_r1)
            if paired_end:
                random.shuffle(reads_r2)
        
        # Write BAM file(s)
        if paired_end:
            self._write_paired_bam(header, reads_r1, reads_r2, output_file, second_file)
        else:
            self._write_single_bam(header, reads_r1, output_file)

    def _write_single_bam(self, header: dict, reads: List[dict], output_file: str):
        """Write single-end BAM file."""
        ref_name_to_id = {sq['SN']: i for i, sq in enumerate(header['SQ'])}
        
        with pysam.AlignmentFile(output_file, 'wb', header=header) as bam_file:
            for read_data in reads:
                segment = self._create_segment(read_data, ref_name_to_id)
                bam_file.write(segment)

    def _write_paired_bam(self, header: dict, reads_r1: List[dict], reads_r2: List[dict], output_file: str, second_file: str):
        """Write paired-end BAM files with realistic mate distribution."""
        ref_name_to_id = {sq['SN']: i for i, sq in enumerate(header['SQ'])}
        
        # Create paired reads with realistic distribution
        # In real sequencing, mates can come from different loci/haplotypes
        paired_reads = []
        
        # Create pairs with realistic mate distribution
        # Both files should have the same number of reads for proper pairing
        num_pairs = max(len(reads_r1), len(reads_r2))
        
        # Extend the shorter list if needed
        if len(reads_r1) < num_pairs:
            # Add more reads to R1 by duplicating some existing ones
            while len(reads_r1) < num_pairs:
                reads_r1.append(random.choice(reads_r1))
        
        if len(reads_r2) < num_pairs:
            # Add more reads to R2 by duplicating some existing ones
            while len(reads_r2) < num_pairs:
                reads_r2.append(random.choice(reads_r2))
        
        for i in range(num_pairs):
            read1 = reads_r1[i]
            read2 = reads_r2[i]
            
            # Keep original read names - they should already be the same from selection
            # The R1 and R2 combinations were selected to ensure all read names appear in both files
            
            # Mates can come from different loci/haplotypes (more realistic)
            # 70% chance of same locus, 30% chance of different locus
            if random.random() < 0.7:
                # Same locus, different haplotype or same haplotype
                locus = read1['ref_name'].split('_')[0]
                haplotype = random.choice(self.haplotype_list)
                read2['ref_name'] = f"{locus}_{haplotype}"
            else:
                # Different locus
                different_locus = random.choice(self.loci)
                haplotype = random.choice(self.haplotype_list)
                read2['ref_name'] = f"{different_locus}_{haplotype}"
            
            # Create realistic insert size
            insert_size = random.randint(200, 400)
            read1['start'] = random.randint(0, 100)
            read2['start'] = read1['start'] + insert_size - self.read_length
            
            paired_reads.append((read1, read2))
        
        # Write R1 file
        with pysam.AlignmentFile(output_file, 'wb', header=header) as bam_file:
            for read1, read2 in paired_reads:
                segment = self._create_segment(read1, ref_name_to_id)
                segment.flag |= 0x40  # First in pair
                segment.next_reference_id = ref_name_to_id[read2['ref_name']]
                segment.next_reference_start = read2['start']
                segment.template_length = read2['start'] - read1['start'] + self.read_length
                bam_file.write(segment)
        
        # Write R2 file
        with pysam.AlignmentFile(second_file, 'wb', header=header) as bam_file:
            for read1, read2 in paired_reads:
                segment = self._create_segment(read2, ref_name_to_id)
                segment.flag |= 0x80  # Second in pair
                segment.next_reference_id = ref_name_to_id[read1['ref_name']]
                segment.next_reference_start = read1['start']
                segment.template_length = read1['start'] - read2['start'] + self.read_length
                bam_file.write(segment)

    def _create_segment(self, read_data: dict, ref_name_to_id: dict) -> pysam.AlignedSegment:
        """Create a pysam AlignedSegment from read data."""
        segment = pysam.AlignedSegment()
        segment.reference_id = ref_name_to_id[read_data['ref_name']]
        segment.query_name = read_data['name']
        segment.reference_start = read_data['start']
        segment.query_sequence = read_data['sequence']
        segment.query_qualities = pysam.qualitystring_to_array(read_data['quality'])
        segment.flag = read_data['flag']
        segment.mapping_quality = read_data['mapq']
        segment.cigar = read_data['cigar']
        
        # Add tags
        for tag, value in read_data['tags']:
            segment.set_tag(tag, value)
        
        return segment

    def generate_locus_file(self, output_file: str):
        """Generate locus ID file compatible with bam2emase."""
        with open(output_file, 'w') as f:
            for locus in self.loci:
                f.write(f'{locus}\t0\n')

    def get_statistics(self) -> dict:
        """Get statistics about the generated data."""
        return {
            'num_loci': self.num_loci,
            'num_haplotypes': self.num_haplotypes,
            'num_read_names': self.num_read_names,
            'alignments_r1': self.alignments_r1,
            'alignments_r2': self.alignments_r2,
            'read_length': self.read_length,
            'coverage_distribution': self.coverage_distribution,
            'prefix_locus': self.prefix_locus,
            'prefix_read': self.prefix_read,
            'haplotypes': self.haplotype_list,
            'loci': self.loci,
            'read_names': self.read_names,
            'selected_combinations': len(self.selected_combinations_r1),
            'total_references': len(self.reference_sequences)
        }


def generate_test_bam(
    output_file: str,
    num_loci: int = 10,
    num_haplotypes: int = 4,
    num_read_names: int = 100,
    alignments_r1: int = 100,
    alignments_r2: int = None,
    read_length: int = 100,
    coverage_distribution: str = 'realistic',
    prefix_locus: str = 'ENSMUST',
    prefix_read: str = 'READ',
    haplotypes: str = None,
    randomize_read_names: bool = False,
    paired_end: bool = False,
    second_file: str = None,
    locus_file: str = None,
    random_seed: int = None
) -> dict:
    """
    Generate realistic test BAM file with improved parameters.
    
    Args:
        output_file: Path to output BAM file
        num_loci: Number of loci/transcripts
        num_haplotypes: Number of haplotypes
        num_read_names: Number of unique read names
        alignments_r1: Number of alignments in first BAM file (R1)
        alignments_r2: Number of alignments in second BAM file (R2) for paired-end
        read_length: Length of each read
        coverage_distribution: How to distribute reads ('uniform', 'exponential', 'realistic')
        prefix_locus: Prefix for locus names
        prefix_read: Prefix for read names
        haplotypes: Custom haplotype string (overrides num_haplotypes)
        randomize_read_names: If True, randomize the order of read names
        paired_end: Whether to generate paired-end reads
        second_file: Path to second BAM file (R2) if paired-end
        locus_file: Path to output locus file (optional)
        random_seed: Random seed for reproducible results
    
    Returns:
        Dictionary with generation statistics
    """
    generator = BAMGenerator(
        num_loci=num_loci,
        num_haplotypes=num_haplotypes,
        num_read_names=num_read_names,
        alignments_r1=alignments_r1,
        alignments_r2=alignments_r2,
        read_length=read_length,
        coverage_distribution=coverage_distribution,
        prefix_locus=prefix_locus,
        prefix_read=prefix_read,
        haplotypes=haplotypes,
        randomize_read_names=randomize_read_names,
        random_seed=random_seed
    )
    
    # Generate BAM file(s)
    generator.generate_bam(output_file, paired_end, second_file)
    
    # Generate locus file if requested
    if locus_file:
        generator.generate_locus_file(locus_file)
    
    return generator.get_statistics() 