#!/usr/bin/env python3
"""
population_target_primers.py

This identifies shared genomic regions of desired genomes while ensuring they are absent
from desired excluded genomes. Then primers are designed on these regions.
"""
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .population_target_primers_utils import (
    shred_genome_parallel,
    map_shreds_minimap2,
    compute_candidate_regions)


def parse_genome_list(file_path):
    """
    Parse genome list file and determine reference, include and exclude genomes for primer design
    """
    reference = None
    include = []
    exclude = []
    with open(file_path, encoding='utf-8') as f:
        for line in f:
            if not line.strip():
                continue
            path, status = line.strip().split("\t")
            status = status.lower()
            if status == "reference":
                reference = path
            elif status == "include":
                include.append(path)
            elif status == "exclude":
                exclude.append(path)
    if reference is None:
        raise ValueError("No reference genome specified in genome_list.")
    return reference, include, exclude


def merge_contiguous_candidates(depth_table):
    """
    Merge candidate windows that are contiguous AND have identical include/exclude proportions.
    Updated to handle the new sorting criteria including region_size.
    """
    if not depth_table:
        return []

    merged = []
    current = depth_table[0].copy()

    for entry in depth_table[1:]:
        # Check if same contig, contiguous, and identical proportions
        if (entry["contig"] == current["contig"] and entry["start"] == current["end"] + 1 and entry["include_genome_proportion_present"] == current["include_genome_proportion_present"] and entry["exclude_genome_proportion_present"] == current["exclude_genome_proportion_present"]):
            # Extend current region
            current["end"] = entry["end"]
            current["region_size"] = current["end"] - current["start"] + 1  # Update region size
        else:
            merged.append(current)
            current = entry.copy()

    merged.append(current)

    # Re-sort after merging to ensure proper ordering with updated region sizes
    merged.sort(key=lambda x: (
        -x["include_genome_proportion_present"],
        x["exclude_genome_proportion_present"],
        -x["region_size"]
    ))

    return merged


def write_outputs(depth_table, detailed_table, out_prefix):
    """
    Write tsv outputs for identifying regions of the reference genome that are suitable for primer design
    """
    candidate_tsv = out_prefix + "_candidates.tsv"
    detailed_tsv = out_prefix + "_detailed.tsv"

    print(f"[INFO] Writing candidate table to {candidate_tsv}")
    with open(candidate_tsv, "w", encoding='utf-8') as c_out:
        headers = ["reference_genome", "contig", "contig_length", "start", "end", "region_size",
                   "include_genome_proportion_present", "exclude_genome_proportion_present"]
        c_out.write("\t".join(headers) + "\n")
        for row in depth_table:
            c_out.write("\t".join(str(row[h]) for h in headers) + "\n")

    print(f"[INFO] Writing detailed binary table to {detailed_tsv}")
    if detailed_table:
        genome_columns = [k for k in detailed_table[0].keys() if k not in ("contig", "start", "end")]
        with open(detailed_tsv, "w", encoding='utf-8') as d_out:
            d_out.write("\t".join(["contig", "start", "end"] + genome_columns) + "\n")
            for row in detailed_table:
                d_out.write("\t".join(str(row.get(c, 0)) for c in ["contig", "start", "end"] + genome_columns) + "\n")


def main(args=None):
    """
    Main runner for identifying population level target primers - Parallelized version
    """

    parser = argparse.ArgumentParser(description="Population-targeted primer candidate identification (shred-based, parallelized)")
    parser.add_argument("--genome_list", required=True, help="TSV: genome_path \\t reference/include/exclude")
    parser.add_argument("--shred_length", type=int, default=1000, help="Length of shredded fragments to consider")
    parser.add_argument("--step_size", type=int, default=None, help="Step size for shredding (default = shred_length)")
    parser.add_argument("--max_candidate_seqs", type=int, default=10, help="Maximum number of candidate sequences to return")
    parser.add_argument("--min_include_prop", type=float, default=0.95, help="Minimum proportion of include genomes present")
    parser.add_argument("--min_contig_size", type=int, default=2000, help="Minimum contig size to consider")
    parser.add_argument("--threads", type=int, default=2, help="Number of threads for parallel processing")
    parser.add_argument("--output", required=True, help="Prefix for output candidate files")

    if args is None or isinstance(args, list):
        args = parser.parse_args(args)

    tmp_dir = os.path.join(os.path.dirname(args.output), "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    reference, include, exclude = parse_genome_list(args.genome_list)

    print(f"[INFO] Starting parallelized shred-based primer candidates with {args.threads} threads")
    print(f"[INFO] Reference: {reference}")
    print(f"[INFO] Include genomes: {len(include)}")
    print(f"[INFO] Exclude genomes: {len(exclude)}")

    # Load reference genome and filter small contigs
    ref_seq_records = [rec for rec in SeqIO.parse(reference, "fasta") if len(rec.seq) >= args.min_contig_size]
    if not ref_seq_records:
        raise ValueError("No contigs meet the --min_contig_size requirement")

    print(f"[INFO] Reference genome has {len(ref_seq_records)} contigs >= {args.min_contig_size} bp")

    # Shred include + exclude genomes in parallel
    print("[INFO] Starting parallel genome shredding...")
    all_genomes = include + exclude
    shredded_genomes = shred_genome_parallel(
        all_genomes,
        args.shred_length,
        args.step_size or args.shred_length,
        args.threads
    )

    total_shreds = len(shredded_genomes)
    print(f"[INFO] Total shreds generated: {total_shreds}")

    # Map shredded fragments to the reference genome
    print("[INFO] Starting parallel mapping with minimap2...")
    paf_results = map_shreds_minimap2(shredded_genomes, [reference], args.threads, tmp_dir)

    # Compute candidate regions based on mapping
    print("[INFO] Computing candidate regions...")
    depth_table, detailed_table = compute_candidate_regions(
        ref_seq_records,
        paf_results,
        reference,
        include,
        exclude,
        shredded_genomes=shredded_genomes,
        shred_length=args.shred_length,
        min_include_prop=args.min_include_prop,
        min_contig_size=args.min_contig_size
    )

    print(f"[INFO] Found {len(depth_table)} candidate regions before merging")

    # Merge contiguous regions with identical proportions
    print("[INFO] Merging contiguous regions...")
    depth_table = merge_contiguous_candidates(depth_table)

    print(f"[INFO] Found {len(depth_table)} candidate regions after merging")

    # Limit number of candidate sequences
    depth_table = depth_table[:args.max_candidate_seqs]

    print(f"[INFO] Returning top {len(depth_table)} candidate regions")

    # Write outputs
    write_outputs(depth_table, detailed_table, args.output)

    print("[INFO] Primer regions calculated successfully!")

    # Print summary of top candidates
    if depth_table:
        print("\n[INFO] Top candidate regions:")
        for i, candidate in enumerate(depth_table[:5], 1):
            print(f"  {i}. {candidate['contig']}:{candidate['start']}-{candidate['end']} "
                  f"(size: {candidate['region_size']}, "
                  f"include: {candidate['include_genome_proportion_present']:.3f}, "
                  f"exclude: {candidate['exclude_genome_proportion_present']:.3f})")


if __name__ == "__main__":
    main()
