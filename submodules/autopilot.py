#!/usr/bin/env python3
"""
autopilot.py

Automatically runs population_target_primers and primer design
with sliding-window logic for optimal primer selection.
"""

import os
import subprocess
# import argparse
import sys
import glob
# import shutil
from Bio import SeqIO


def get_genome_length(fasta_file):
    """Return total length of all sequences in a FASTA file."""
    mb_length = 0
    with open(fasta_file, encoding="utf-8") as f:
        for record in SeqIO.parse(f, "fasta"):
            mb_length += len(record.seq) / 1e6
    return mb_length


def cleanup_failed_design(design_out, output_dir):
    """
    Remove failed PCR primer files on the quest for correct primers
    """
    # Delete main TSV
    if os.path.exists(design_out):
        os.remove(design_out)

    # Delete associated primer FASTA files
    for fasta_file in glob.glob(f"{design_out}_pair*.fasta"):
        os.remove(fasta_file)

    # Delete associated in_silico_pcr outputs
    in_silico_dir = os.path.join(output_dir, "in_silico_pcr")
    if os.path.isdir(in_silico_dir):
        for tsv_file in glob.glob(os.path.join(in_silico_dir, "pair*_in_silico_pcr.tsv")):
            os.remove(tsv_file)

    # Optional: clean tmp files for this design
    tmp_dir = os.path.join(output_dir, "tmp")
    if os.path.isdir(tmp_dir):
        for tmp_file in glob.glob(os.path.join(tmp_dir, f"{os.path.basename(design_out)}*")):
            os.remove(tmp_file)


def run_autopilot(args):
    """
    Autopilot main runner
    """
    print("""
          ==================================================================================================================================
          Welcome to autopilot mode - They asked me how well I understood theoretical biology. I said I had a theoretical degree in biology.
          ==================================================================================================================================
          """)
    os.makedirs(args.output_dir, exist_ok=True)

    genome_list_path = os.path.abspath(args.genome_list)
    if not os.path.exists(genome_list_path):
        raise FileNotFoundError(f"Genome list file not found: {genome_list_path}")

    with open(genome_list_path, encoding="utf-8") as f:
        first_line = f.readline().strip().split("\t")
        ref_genome_path = first_line[0]  # first genome in list assumed reference
    ref_genome_length = get_genome_length(ref_genome_path)

    # Run population_target_primers
    pop_prefix = os.path.join(args.output_dir, "population")
    pop_cmd = [
        sys.executable, "population_primer_designer.py", "population_target_primers",
        "--genome_list", genome_list_path,
        "--output", pop_prefix,
        "--threads", str(args.threads),
    ]
    subprocess.run(pop_cmd, check=True)

    # Read population output and sort contigs by stretch length
    pop_out = f"{pop_prefix}_candidates.tsv"
    candidates = []
    with open(pop_out, encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
        col_idx = {name: i for i, name in enumerate(header)}
        for line in f:
            parts = line.strip().split("\t")
            start = int(parts[col_idx["start"]])
            end = int(parts[col_idx["end"]])
            stretch_len = end - start
            candidates.append({
                "contig": parts[col_idx["contig"]],
                "start": start,
                "end": end,
                "length": stretch_len,
                "reference_genome": parts[col_idx["reference_genome"]],
            })

    if not candidates:
        raise ValueError("No candidates found in population_target_primers output")

    # Sort descending by length
    candidates.sort(key=lambda x: x["length"], reverse=True)

    # Try each contig until a design succeeds
    success = False
    for candidate in candidates:
        contig = candidate["contig"]
        contig_start = candidate["start"]
        contig_end = candidate["end"]
        start_coord = contig_start

        print(f"Trying contig {contig} from {candidate['reference_genome']} (stretch length {candidate['length']})")

        while start_coord < contig_end:
            stop_coord = start_coord + args.amplicon_length

            # Safety check: cap stop_coord to contig_end
            if stop_coord > contig_end:
                stop_coord = contig_end
                # If remaining region is too short, break
                if stop_coord - start_coord < 50:  # arbitrary minimal amplicon length
                    print(f"Remaining region too short ({stop_coord - start_coord} bp). Moving to next contig.")
                    break

            design_out = os.path.join(
                args.output_dir, f"primer_design_{contig}_{start_coord}_{stop_coord}"
            )
            design_cmd = [
                sys.executable, "population_primer_designer.py", "design",
                "--genome_list", genome_list_path,
                "--contig", contig,
                "--start_coord", str(start_coord),
                "--stop_coord", str(stop_coord),
                "--output", design_out,
                "--threads", str(args.threads),
            ]

            if args.max_amplicon_diff is not None:
                design_cmd += ["--max_amplicon_diff", str(args.max_amplicon_diff)]
            if args.num_return is not None:
                design_cmd += ["--num_return", str(args.num_return)]

            print(f"Running design: {start_coord}-{stop_coord}")
            subprocess.run(design_cmd, check=True)

            # Check if design produced primers (TSV has more than header)
            with open(f'{design_out}.tsv', encoding="utf-8") as f:
                lines = f.readlines()
            if len(lines) > 1:
                print(f"Success! Primers found at {contig}:{start_coord}-{stop_coord}")
                print("""
          ==================================================================================================================================
          Yeah! Who won the lottery? I did!
          ==================================================================================================================================
          """)

                success = True
                break
            else:
                # Delete the failed TSV
                cleanup_failed_design(design_out, args.output_dir)
                start_coord += 100
                print(f"No primers found. Sliding to {start_coord}")
                print("""
          ==================================================================================================================================
          Patrolling these genomes almost makes you wish for a nuclear winter
          ==================================================================================================================================
          """)

        if success:
            break
        else:
            print(f"Could not find primers on contig {contig}. Trying next candidate.")

    if not success:
        print(f"""
        ==================================================================================================================================
        From where you're kneeling it must seem like a {ref_genome_length:.2f}-Mb run of bad luck.

        Truth is... the genome was rigged from the start.
        ==================================================================================================================================
        """)
        raise RuntimeError("Autopilot could not find primers on any contig.")
