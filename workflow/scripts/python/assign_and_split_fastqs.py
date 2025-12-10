#!/usr/bin/env python3
"""
Assign protein labels to clusters from FASTQ files and split by assignment.

This script implements the assign_label logic on unaligned FASTQ files:
1. Parse barcodes from DPM and BPM FASTQ read names
2. Group reads by cluster barcode
3. Count BPM antibody representations per cluster
4. Apply thresholds to assign protein labels
5. Split DPM reads into per-target FASTQ files with reformatted headers

Headers are formatted with cluster barcode (CB), read group (RG), and pseudo-UMI:
    @ReadID:CB:ActualBarcode:RG:TargetID_PseudoUMI

Example:
    @A00123:45:H3GYYDRX3:1:1101:1234:5678:CB:Y10.E4.O62.sample1:RG:IgG_ATCGTAGCATCGATCG
"""

import argparse
import collections
import gzip
import hashlib
import os
import re
import sys
from typing import Dict, List, Tuple

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import helpers


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Assign protein labels to clusters and split FASTQ files by assignment"
    )
    parser.add_argument(
        "--dpm_fastqs_r1",
        nargs="+",
        required=True,
        help="R1 DPM FASTQ files (barcoded_dpm.fastq.gz)"
    )
    parser.add_argument(
        "--dpm_fastqs_r2",
        nargs="+",
        required=True,
        help="R2 barcoded FASTQ files (mate pairs)"
    )
    parser.add_argument(
        "--bpm_fastqs",
        nargs="+",
        required=True,
        help="BPM FASTQ files (barcoded_bpm.fastq.gz)"
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory for split FASTQ files"
    )
    parser.add_argument(
        "--sample",
        required=True,
        help="Sample name"
    )
    parser.add_argument(
        "--min_oligos",
        type=int,
        required=True,
        help="Minimum number of BPM oligos to call a cluster"
    )
    parser.add_argument(
        "--proportion",
        type=float,
        required=True,
        help="Minimum representation proportion of oligos to call a cluster"
    )
    parser.add_argument(
        "--max_size",
        type=int,
        required=True,
        help="Maximum cluster size (DPM count) to keep"
    )
    parser.add_argument(
        "--num_tags",
        type=int,
        required=True,
        help="Number of tags in barcode (including DPM/BEAD tags)"
    )
    parser.add_argument(
        "--pseudo_umi_length",
        type=int,
        default=16,
        help="Length of pseudo-UMI in bp"
    )
    parser.add_argument(
        "--targets",
        nargs="+",
        required=True,
        help="Expected target names (ensures all output files are created)"
    )
    parser.add_argument(
        "--output_stats",
        required=True,
        help="Output statistics file"
    )
    return parser.parse_args()


def generate_pseudo_umi(cluster_barcode: str, length: int = 16) -> str:
    """
    Generate IUPAC-compliant pseudo-UMI from cluster barcode hash.

    Args:
        cluster_barcode: Full cluster barcode string (e.g., "Y10.E4.O62.sample1")
        length: Length of pseudo-UMI in bp

    Returns:
        Pseudo-UMI string using only ATGC alphabet
    """
    hash_obj = hashlib.sha256(cluster_barcode.encode())
    hash_bytes = hash_obj.digest()
    # Convert bytes to DNA alphabet (ATGC only)
    bases = ['A', 'T', 'C', 'G']
    umi = ''.join(bases[b % 4] for b in hash_bytes[:length])
    return umi


def parse_barcode_from_fastq_name(read_name: str, num_tags: int) -> Tuple[str, str]:
    """
    Parse barcode from FASTQ read name.

    Args:
        read_name: FASTQ read name (may include @ prefix)
        num_tags: Number of barcode tags expected

    Returns:
        (read_type, cluster_barcode) tuple
        Example: 'read1::[DPM1][Y10][E4][O62]' -> ('DPM1', 'Y10.E4.O62')

    Returns:
        (None, None) if no barcode pattern found
    """
    # Remove @ prefix if present
    if read_name.startswith('@'):
        read_name = read_name[1:]

    # Pattern from extract_barcode_to_tags.py line 148
    pattern = re.compile("::" + num_tags * r"\[([a-zA-Z0-9_-]+)\]")
    match = pattern.search(read_name)
    if match:
        read_type = match.group(1)
        barcode = '.'.join(match.groups()[1:])
        return read_type, barcode
    return None, None


def assign_clusters_from_fastq(
    dpm_fastqs_r1: List[str],
    bpm_fastqs: List[str],
    min_oligos: int,
    proportion: float,
    max_size: int,
    num_tags: int,
    sample: str
) -> Dict[str, str]:
    """
    Assign protein labels to clusters based on BPM read counts.

    This implements the logic from assign_label.py lines 280-339.

    Args:
        dpm_fastqs_r1: List of R1 DPM FASTQ file paths
        bpm_fastqs: List of BPM FASTQ file paths
        min_oligos: Minimum BPM count for assignment
        proportion: Minimum proportion of dominant antibody
        max_size: Maximum DPM count per cluster
        num_tags: Number of barcode tags
        sample: Sample name to append to cluster barcodes

    Returns:
        Dictionary mapping cluster_barcode -> assigned_label
    """
    # Read all BPM reads, group by CB
    bpm_counts = {}  # {cluster_barcode: {antibody_id: count}}
    bpm_read_count = 0
    bpm_parsed_count = 0

    print(f"Processing {len(bpm_fastqs)} BPM FASTQ files...", file=sys.stderr)
    for bpm_fastq in bpm_fastqs:
        with helpers.file_open(bpm_fastq, 'rt') as f:
            for name, seq, thrd, qual in helpers.fastq_parse(f):
                bpm_read_count += 1
                rt, cb = parse_barcode_from_fastq_name(name, num_tags)
                if rt and rt.startswith("BEAD_"):
                    antibody_id = rt  # e.g., "BEAD_AB1-A1"
                    cb_with_sample = f"{cb}.{sample}"
                    if cb_with_sample not in bpm_counts:
                        bpm_counts[cb_with_sample] = {}
                    bpm_counts[cb_with_sample][antibody_id] = bpm_counts[cb_with_sample].get(antibody_id, 0) + 1
                    bpm_parsed_count += 1
                if bpm_read_count % 100000 == 0:
                    print(f"  Processed {bpm_read_count} BPM reads...", file=sys.stderr)

    print(f"Total BPM reads: {bpm_read_count}, Parsed: {bpm_parsed_count}", file=sys.stderr)
    print(f"Unique clusters with BPM reads: {len(bpm_counts)}", file=sys.stderr)

    # Count DPM reads per cluster
    dpm_counts = {}  # {cluster_barcode: count}
    dpm_read_count = 0
    dpm_parsed_count = 0

    print(f"Processing {len(dpm_fastqs_r1)} R1 DPM FASTQ files...", file=sys.stderr)
    for dpm_fastq in dpm_fastqs_r1:
        with helpers.file_open(dpm_fastq, 'rt') as f:
            for name, seq, thrd, qual in helpers.fastq_parse(f):
                dpm_read_count += 1
                rt, cb = parse_barcode_from_fastq_name(name, num_tags)
                if rt and rt.startswith("DPM"):
                    cb_with_sample = f"{cb}.{sample}"
                    dpm_counts[cb_with_sample] = dpm_counts.get(cb_with_sample, 0) + 1
                    dpm_parsed_count += 1
                if dpm_read_count % 100000 == 0:
                    print(f"  Processed {dpm_read_count} DPM reads...", file=sys.stderr)

    print(f"Total DPM reads: {dpm_read_count}, Parsed: {dpm_parsed_count}", file=sys.stderr)
    print(f"Unique clusters with DPM reads: {len(dpm_counts)}", file=sys.stderr)

    # Apply assignment logic (from assign_label.py lines 280-339)
    assignments = {}  # {cluster_barcode: assigned_label}
    assignment_counts = collections.Counter()

    print("Assigning protein labels to clusters...", file=sys.stderr)
    for cb in dpm_counts:
        dpm_count = dpm_counts[cb]

        # Check max_size threshold
        if dpm_count > max_size:
            assignments[cb] = "filtered"
            assignment_counts['filtered'] += 1
            continue

        # Get BPM counts for this cluster
        if cb not in bpm_counts or len(bpm_counts[cb]) == 0:
            assignments[cb] = "none"
            assignment_counts['none'] += 1
            continue

        # Find most represented antibody
        antibody_counts = bpm_counts[cb]
        max_antibody = max(antibody_counts, key=antibody_counts.get)
        max_count = antibody_counts[max_antibody]
        total_bpm = sum(antibody_counts.values())

        # Check min_oligos threshold
        if max_count < min_oligos:
            assignments[cb] = "uncertain"
            assignment_counts['uncertain'] += 1
            continue

        # Check proportion threshold
        max_proportion = max_count / total_bpm if total_bpm > 0 else 0
        if max_proportion < proportion:
            assignments[cb] = "ambiguous"
            assignment_counts['ambiguous'] += 1
            continue

        # Assign to antibody (strip "BEAD_" prefix)
        assigned_label = max_antibody.replace("BEAD_", "")
        assignments[cb] = assigned_label
        assignment_counts[assigned_label] += 1

    print("Assignment counts:", file=sys.stderr)
    for label, count in sorted(assignment_counts.items()):
        print(f"  {label}: {count} clusters", file=sys.stderr)

    return assignments


def split_fastqs_by_assignment(
    dpm_fastqs_r1: List[str],
    dpm_fastqs_r2: List[str],
    assignments: Dict[str, str],
    output_dir: str,
    sample: str,
    num_tags: int,
    pseudo_umi_length: int,
    targets: List[str]
) -> Dict[str, int]:
    """
    Split DPM FASTQ files by protein assignment.

    R1 files contain barcodes and determine the assignment.
    R2 files are mate pairs that must be matched by read name to their R1 counterparts.

    Args:
        dpm_fastqs_r1: List of R1 DPM FASTQ file paths (have barcodes)
        dpm_fastqs_r2: List of R2 barcoded FASTQ file paths (mate pairs, no DPM barcodes)
        assignments: Dictionary mapping cluster_barcode -> assigned_label
        output_dir: Output directory for split FASTQ files
        sample: Sample name
        num_tags: Number of barcode tags
        pseudo_umi_length: Length of pseudo-UMI
        targets: List of expected targets (to create empty files if needed)

    Returns:
        Dictionary with read counts per target
    """
    os.makedirs(output_dir, exist_ok=True)

    print(f"R1 files: {len(dpm_fastqs_r1)}, R2 files: {len(dpm_fastqs_r2)}", file=sys.stderr)

    # Process R1 and R2 separately
    read_counts = {}

    # Process R1 - these have barcodes
    for read_num, fastq_list in [(1, dpm_fastqs_r1), (2, dpm_fastqs_r2)]:
        print(f"Processing R{read_num} files...", file=sys.stderr)

        # Open output file handles for each target
        output_handles = {}  # {target: file_handle}
        target_counts = collections.Counter()

        for target in targets:
            filename = f"{sample}.{target}_R{read_num}.fastq.gz"
            filepath = os.path.join(output_dir, filename)
            output_handles[target] = gzip.open(filepath, 'wt')

        reads_processed = 0
        reads_written = 0

        for fastq_file in fastq_list:
            with helpers.file_open(fastq_file, 'rt') as f:
                for name, seq, thrd, qual in helpers.fastq_parse(f):
                    reads_processed += 1

                    # Parse barcode from read name
                    # R1 DPM files: read type is "DPM1", "DPM2", etc.
                    # R2 barcoded files: read type could be anything with barcode
                    rt, cb = parse_barcode_from_fastq_name(name, num_tags)

                    # Skip if no barcode found
                    if not rt or not cb:
                        continue

                    # For R1: only process DPM reads (already filtered by split_bpm_dpm)
                    # For R2: process all reads with valid barcodes (these are paired mates)
                    if read_num == 1 and not rt.startswith("DPM"):
                        # R1 but not DPM - skip (shouldn't happen in barcoded_dpm.fastq.gz files)
                        continue

                    cb_with_sample = f"{cb}.{sample}"
                    assigned_label = assignments.get(cb_with_sample, "none")

                    # Generate pseudo-UMI
                    pseudo_umi = generate_pseudo_umi(cb_with_sample, pseudo_umi_length)

                    # Format header
                    # Original: @read_id::[TAG1][Y10][E4]...
                    # Remove barcode part, add CB:RG:UMI
                    clean_read_id = name.split("::")[0]
                    if clean_read_id.startswith('@'):
                        clean_read_id = clean_read_id[1:]
                    new_header = f"{clean_read_id}:CB:{cb_with_sample}:RG:{assigned_label}_{pseudo_umi}"

                    # Write to appropriate file
                    if assigned_label in output_handles:
                        fh = output_handles[assigned_label]
                        fh.write(f"@{new_header}\n{seq}\n+\n{qual}\n")
                        target_counts[assigned_label] += 1
                        reads_written += 1
                    else:
                        print(f"Warning: Assigned label '{assigned_label}' not in targets list", file=sys.stderr)

                    if reads_processed % 100000 == 0:
                        print(f"  Processed {reads_processed} R{read_num} reads, written {reads_written}...", file=sys.stderr)

        # Close all file handles
        for fh in output_handles.values():
            fh.close()

        print(f"R{read_num} complete: {reads_processed} processed, {reads_written} written", file=sys.stderr)
        print(f"R{read_num} counts by target:", file=sys.stderr)
        for target in sorted(targets):
            count = target_counts[target]
            print(f"  {target}: {count} reads", file=sys.stderr)
            read_counts[f"{target}_R{read_num}"] = count

    return read_counts


def main():
    args = parse_args()

    print("=" * 80, file=sys.stderr)
    print("assign_and_split_fastqs.py", file=sys.stderr)
    print("=" * 80, file=sys.stderr)
    print(f"Sample: {args.sample}", file=sys.stderr)
    print(f"DPM R1 FASTQ files: {len(args.dpm_fastqs_r1)}", file=sys.stderr)
    print(f"DPM R2 FASTQ files: {len(args.dpm_fastqs_r2)}", file=sys.stderr)
    print(f"BPM FASTQ files: {len(args.bpm_fastqs)}", file=sys.stderr)
    print(f"Output directory: {args.output_dir}", file=sys.stderr)
    print(f"Targets: {', '.join(args.targets)}", file=sys.stderr)
    print(f"Parameters: min_oligos={args.min_oligos}, proportion={args.proportion}, max_size={args.max_size}", file=sys.stderr)
    print(f"Pseudo-UMI length: {args.pseudo_umi_length}", file=sys.stderr)
    print("=" * 80, file=sys.stderr)

    # Step 1: Assign clusters based on BPM/DPM counts
    assignments = assign_clusters_from_fastq(
        args.dpm_fastqs_r1,
        args.bpm_fastqs,
        args.min_oligos,
        args.proportion,
        args.max_size,
        args.num_tags,
        args.sample
    )

    # Step 2: Split DPM FASTQ files by assignment
    read_counts = split_fastqs_by_assignment(
        args.dpm_fastqs_r1,
        args.dpm_fastqs_r2,
        assignments,
        args.output_dir,
        args.sample,
        args.num_tags,
        args.pseudo_umi_length,
        args.targets
    )

    # Step 3: Write statistics
    print(f"Writing statistics to {args.output_stats}", file=sys.stderr)
    with open(args.output_stats, 'w') as f:
        f.write("target\tread_count\n")
        for target_read, count in sorted(read_counts.items()):
            f.write(f"{target_read}\t{count}\n")

    print("=" * 80, file=sys.stderr)
    print("Complete!", file=sys.stderr)
    print("=" * 80, file=sys.stderr)


if __name__ == "__main__":
    main()