#!/usr/bin/env python3
"""Thin wrapper forwarding to gbrs.tools.explore_outputs.cli_main."""

import argparse            # Command line argument parsing - provides user-friendly CLI interface
import os

from gbrs.tools.explore_outputs import analyze_genoprobs_file
from gbrs.tools.explore_outputs import analyze_genotypes_file
from gbrs.tools.explore_outputs import analyze_genotypes_file_tsv
from gbrs.tools.explore_outputs import analyze_alignment_counts_file
from gbrs.tools.explore_outputs import analyze_expected_read_counts_file
from gbrs.tools.explore_outputs import analyze_tpm_file


def main():
    # ===== ARGUMENT PARSER SETUP SECTION =====
    # Create an argument parser with comprehensive help information
    # This provides users with clear guidance on how to use the script
    parser = argparse.ArgumentParser(
        description="Explore GBRS file types and their structure",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # ===== INPUT FILE ARGUMENT =====
    # Required argument specifying the file to explore
    # This is the core input that determines what will be analyzed
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Path to the input file to explore'
    )

    # ===== FILE TYPE ARGUMENT =====
    # Required argument specifying the type of file being explored
    # This determines which exploration function will be called
    # The choices are limited to the file types the script can handle
    parser.add_argument(
        '--type',
        required=True,
        choices=[
            'tpm', 'expected_read_counts', 'alignment_counts',
            'genoprobs_npz', 'genotypes_npz', 'genotypes_tsv'
        ],
        help='Type of file to explore'
    )

    # ===== ARGUMENT PARSING SECTION =====
    # Parse the command line arguments provided by the user
    # This will exit with an error if required arguments are missing
    args = parser.parse_args()

    # ===== FILE EXISTENCE VALIDATION SECTION =====
    # Check if the specified input file actually exists
    # This prevents confusing errors when trying to analyze non-existent files
    if not os.path.exists(args.input):
        print(f"ERROR: File {args.input} does not exist!")
        return 1

    # ===== FUNCTION DISPATCH SECTION =====
    # Route the request to the appropriate exploration function
    # Each file type has its own specialized analysis function
    try:
        if args.type == 'tpm':
            # Explore TPM expression files (normalized expression levels)
            analyze_tpm_file(args.input)
        elif args.type == 'expected_read_counts':
            # Explore expected read count files (EMASE algorithm output)
            analyze_expected_read_counts_file(args.input)
        elif args.type == 'alignment_counts':
            # Explore alignment count files (multi-mapping statistics)
            analyze_alignment_counts_file(args.input)
        elif args.type == 'genoprobs_npz':
            # Explore genotype probability files (HMM output)
            analyze_genoprobs_file(args.input)
        elif args.type == 'genotypes_npz':
            # Explore genotype NPZ files (Viterbi algorithm output)
            analyze_genotypes_file(args.input)
        elif args.type == 'genotypes_tsv':
            # Explore genotype TSV files (human-readable format)
            analyze_genotypes_file_tsv(args.input)
        else:
            # This should never happen due to argparse choices, but included for safety
            print(f"Unknown file type: {args.type}")
            return 1

    except Exception as e:
        # ===== ERROR HANDLING SECTION =====
        # Catch and report any errors that occur during file exploration
        # This provides users with helpful error messages instead of crashes
        print(f"ERROR exploring file: {e}")
        return 1

    # ===== SUCCESS EXIT SECTION =====
    # Return success code if everything completed without errors
    return 0


# ===== SCRIPT ENTRY POINT =====
# This block ensures the script only runs when executed directly
# It allows the functions to be imported and used in other scripts if needed
if __name__ == "__main__":
    # Call the main function and exit with the appropriate code
    # This ensures proper error reporting and exit status
    exit(main())