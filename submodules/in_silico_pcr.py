#!/usr/bin/env python3
"""
in_silico_pcr.py

Performs in silico PCR.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
import os
from Bio import SeqIO
from .in_silico_pcr_utils import perform_in_silico_pcr_bwa, write_amplicon_sequences_to_fasta, write_amplicons_tsv


def parse_primer_fasta(primer_fasta):
    """
    Parses primer fasta file with _F and _R for in silico PCR
    """
    primers = []
    for rec in SeqIO.parse(primer_fasta, "fasta"):
        name = rec.id
        seq = str(rec.seq)
        if name.endswith("_F"):
            pair_id = "_".join(name.split("_")[:-1])
            primers.append({"pair": pair_id, "forward_name": name, "forward_seq": seq})
        elif name.endswith("_R"):
            pair_id = "_".join(name.split("_")[:-1])
            for p in primers:
                if p["pair"] == pair_id:
                    p.update({"reverse_name": name, "reverse_seq": seq})
                    break
            else:
                primers.append({"pair": pair_id, "reverse_name": name, "reverse_seq": seq})
    return [p for p in primers if "forward_seq" in p and "reverse_seq" in p]


def process_genome(fasta, primer_pairs, bwa_cmd, output_dir):
    """
    Loads genome and runs bwa by mapping primer pairs to perform PCR
    """
    fasta_name = os.path.basename(fasta).rsplit(".", 1)[0]

    # Load sequences
    primers_results = []
    for pair in primer_pairs:
        amplicons = perform_in_silico_pcr_bwa(
            fasta, pair["forward_seq"], pair["reverse_seq"], bwa_cmd=bwa_cmd
        )
        primers_results.append({
            "pair": pair["pair"],
            "forward_name": pair["forward_name"],
            "forward_seq": pair["forward_seq"],
            "reverse_name": pair["reverse_name"],
            "reverse_seq": pair["reverse_seq"],
            "amplicons": amplicons
        })

    # Write amplicons for this genome
    amplicon_dir = os.path.join(output_dir, "amplicons", fasta_name)
    os.makedirs(amplicon_dir, exist_ok=True)
    write_amplicon_sequences_to_fasta(primers_results, fasta, amplicon_dir)

    return fasta_name, primers_results


def main(args=None):
    """
    Main for in silico PCR. Controls script and callout to in_silico_utils.py
    """
    if args is None:
        parser = argparse.ArgumentParser(description="Run in silico PCR with BWA-backtrack using primer FASTA.")
        parser.add_argument("--input_directory", required=True)
        parser.add_argument("--primers", required=True)
        parser.add_argument("--output", required=True)
        parser.add_argument("--bwa", default="bwa")
        parser.add_argument("--threads", type=int, default=4, help="Number of parallel processes [default: 4]")
        args = parser.parse_args()

    fasta_files = [os.path.join(args.input_directory, f)
                   for f in os.listdir(args.input_directory)
                   if f.lower().endswith((".fasta", ".fa", ".fna"))]
    if not fasta_files:
        raise ValueError(f"No FASTA files found in {args.input_directory}")

    primer_pairs = parse_primer_fasta(args.primers)

    all_amplicons_dict = {}
    out_dir = os.path.dirname(args.output) or "."

    # Run in parallel across genomes
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(process_genome, fasta, primer_pairs, args.bwa, out_dir): fasta for fasta in fasta_files}
        for future in as_completed(futures):
            fasta_name, primers_results = future.result()
            all_amplicons_dict[fasta_name] = primers_results

    # Save summary TSV
    output_file = args.output
    if not output_file.endswith(".tsv"):
        output_file += ".tsv"

    write_amplicons_tsv(all_amplicons_dict, output_file)
    print(f"In silico PCR complete! Results saved to {output_file}")


if __name__ == "__main__":
    main()
