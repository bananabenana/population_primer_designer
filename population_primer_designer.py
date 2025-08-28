#!/usr/bin/env python3
"""
population_primer_designer.py

Controller script which runs submodule scripts to generate primers against a population of genomes
and performs in silico PCR as a confirmation
"""

import argparse

# Import submodules from the submodules package
from submodules import design, population_target_primers, in_silico_pcr
from submodules.autopilot import run_autopilot


def main():
    """
    Runs the main controller script and submodule processes
    """
    parser = argparse.ArgumentParser(description="Primer Designer")
    subparsers = parser.add_subparsers(dest="command")

    # -----------------------------
    # design subcommand
    # -----------------------------
    parser_design = subparsers.add_parser(
        "design", help="Design primers for a target region"
    )
    parser_design.add_argument(
        "--input_genome", help="Input genome file (FASTA or GenBank). Required if no --genome_list."
    )
    parser_design.add_argument(
        "--genome_list",
        help="Optional: tab-delimited file with genome path and label (reference/include/exclude)",
        default=None,
    )
    parser_design.add_argument(
        "--gff", help="Optional GFF3 annotation file (needed if using FASTA + --gene/--locus_tag)."
    )
    parser_design.add_argument("--output", required=True, help="Output file for designed primers")
    parser_design.add_argument("--primer3", default="primer3_core", help="Path to primer3_core binary")
    parser_design.add_argument("--bwa", default="bwa", help="Path to bwa binary")
    parser_design.add_argument("--relax", action="store_true", help="Relax constraints")
    parser_design.add_argument("--num_return", type=int, default=5, help="Number of primer pairs to return")
    parser_design.add_argument("--min_tm", type=float, help="Minimum melting temperature")
    parser_design.add_argument("--max_tm", type=float, help="Maximum melting temperature")
    parser_design.add_argument("--start_coord", type=int, help="Start coordinate (1-based)")
    parser_design.add_argument("--stop_coord", type=int, help="Stop coordinate (1-based)")
    parser_design.add_argument("--gene", help="Target by gene name (contig inferred from annotation)")
    parser_design.add_argument("--locus_tag", help="Target by locus_tag (contig inferred from annotation)")
    parser_design.add_argument("--max_amplicon_diff", type=int, help="Maximum amplicon size difference when providing input genome_list - prevents disparate amplicons")
    parser_design.add_argument("--threads", type=int, help="Number of threads to use")
    parser_design.add_argument(
        "--contig",
        help="Contig name (required if using --start_coord/--stop_coord, optional otherwise)"
    )

    # -----------------------------
    # in_silico_pcr subcommand
    # -----------------------------
    parser_pcr = subparsers.add_parser("in_silico_pcr", help="Run in silico PCR using primer FASTA")
    parser_pcr.add_argument("--input_directory", required=True)
    parser_pcr.add_argument("--primers", required=True)
    parser_pcr.add_argument("--output", required=True)
    parser_pcr.add_argument("--tmp_dir", default="tmp_in_silico_pcr")
    parser_pcr.add_argument("--threads", type=int, default=4, help="Number of parallel processes [default: 4]")
    parser_pcr.add_argument("--bwa", default="bwa")

    # -----------------------------
    # population_target_primers subcommand
    # -----------------------------
    parser_pop = subparsers.add_parser(
        "population_target_primers",
        help="Identify candidate population-targeted primer sites using k-mer analysis"
    )
    parser_pop.add_argument("--genome_list", required=True, help="TSV: genome_path \\t reference/include/exclude")
    parser_pop.add_argument("--shred_length", type=int, default=1000, help="Length of shredded fragments to consider")
    parser_pop.add_argument("--step_size", type=int, default=None, help="Step size for shredding (default = shred_length)")
    parser_pop.add_argument("--min_include_prop", type=float, default=0.95, help="Minimum proportion of include genomes present")
    parser_pop.add_argument("--min_contig_size", type=int, default=2000, help="Minimum contig size to consider")
    parser_pop.add_argument("--max_candidate_seqs", type=int, default=20, help="Maximum number of candidate sequences to return")
    parser_pop.add_argument("--threads", type=int, default=2, help="Number of threads for parallel processing")
    parser_pop.add_argument("--output", required=True, help="Prefix for output candidate files")

    # -----------------------------
    # autopilot subcommand
    # -----------------------------
    parser_auto = subparsers.add_parser(
        "autopilot", help="Automatically run population_target_primers then design"
    )
    parser_auto.add_argument("--genome_list", required=True, help="TSV: genome_path \\t reference/include/exclude")
    parser_auto.add_argument("--output_dir", required=True, help="Directory to write all outputs")
    parser_auto.add_argument("--amplicon_length", type=int, help="Target amplicon length")
    parser_auto.add_argument("--threads", type=int, default=4, help="Threads to use")
    parser_auto.add_argument("--max_amplicon_diff", type=int, help="Max allowed difference from target amplicon length")
    parser_auto.add_argument("--num_return", type=int, default=20, help="Number of primer pairs to return")
    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        exit(1)

    # -----------------------------
    # Enforce conditional requirements and genome_list handling
    # -----------------------------
    if args.command == "design":
        if not args.input_genome and not args.genome_list:
            parser.error("Must provide --input_genome if no --genome_list is supplied")
        if (args.start_coord or args.stop_coord) and not args.contig:
            parser.error("--contig is required when using --start_coord/--stop_coord")
        if (args.gene or args.locus_tag):
            fmt = design.detect_format(args.input_genome) if args.input_genome else "fasta"
            if fmt == "fasta" and not args.gff:
                parser.error("--gff is required when using --gene or --locus_tag on a FASTA genome")
        design.main(args)

    elif args.command == "in_silico_pcr":
        in_silico_pcr.main(args)

    elif args.command == "population_target_primers":
        population_target_primers.main(args)

    elif args.command == "autopilot":
        run_autopilot(args)

    else:
        print(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()
