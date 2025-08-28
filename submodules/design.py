#!/usr/bin/env python3
"""
design.py

Designs primers against single genomes and populations and performs in silico PCR as a confirmation
"""

import concurrent.futures
import multiprocessing

import argparse
import os
import tempfile
import subprocess
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from .primer3_utils import (
    generate_primer3_input,
    parse_primer3_output_text,
    write_primer_results_tsv,
    write_primer_fasta_pairs
)
from .in_silico_pcr_utils import (
    perform_in_silico_pcr_bwa,
    write_amplicon_sequences_to_fasta,
    write_amplicons_tsv
)


def detect_format(filepath):
    """
    Detect format of input file, either fasta or genbank.
    """
    with open(filepath, encoding="utf-8") as f:
        first_line = f.readline().strip()
    if first_line.startswith(">"):
        return "fasta"
    elif first_line.startswith("LOCUS"):
        return "genbank"
    else:
        raise ValueError(f"Could not detect format of {filepath}")


def resolve_reference_genome(args, tmp_dir):
    """
    Determine reference genome and return (ref_genome_path, ref_fasta_path, genome_list_dict)

    - Handles single input genome (FASTA/GenBank) via --input_genome
    - Handles genome list via --genome_list, detecting 'reference' label
    - Converts GenBank to FASTA in tmp_dir
    """
    genome_list_dict = {}

    # Process genome list if provided
    if args.genome_list:
        with open(args.genome_list, encoding="utf-8") as f:
            for line in f:
                if not line.strip():
                    continue
                path, label = line.strip().split("\t")
                genome_list_dict[path] = label.lower()

        # Determine reference genome
        ref_candidates = [g for g, l in genome_list_dict.items() if l == "reference"]
        if ref_candidates:
            ref_genome = ref_candidates[0]
        else:
            ref_genome = list(genome_list_dict.keys())[0]  # default first genome
    elif args.input_genome:
        ref_genome = args.input_genome
    else:
        raise ValueError("Must provide --input_genome or a --genome_list with a reference")

    # Ensure FASTA for Primer3 mapping
    if ref_genome.endswith((".gbk", ".gbff")):
        ref_fasta = genbank_to_fasta(ref_genome, tmp_dir)
    else:
        ref_fasta = ref_genome

    return ref_genome, ref_fasta, genome_list_dict


def parse_gff(gff_file):
    """
    Parse a GFF3 file into a dict: contig -> list of features (gene/CDS),
    normalizing qualifiers to lowercase keys and stripping quotes.
    """

    gene_dict = {}
    with open(gff_file, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue
            seqid, _source, ftype, start, end, _score, strand, _phase, attributes = parts

            if ftype not in ("gene", "CDS"):
                continue

            quals = {}
            for attr in attributes.split(";"):
                if "=" in attr:
                    key, val = attr.split("=", 1)
                    key = key.lower().strip()          # normalize to lowercase
                    val = val.strip('"')               # strip quotes if present
                    quals[key] = [val]

            start = int(start) - 1  # convert 1-based to 0-based
            end = int(end)
            loc = FeatureLocation(start, end, strand=1 if strand == "+" else -1)
            feature = SeqFeature(location=loc, type=ftype, qualifiers=quals)

            gene_dict.setdefault(seqid, []).append(feature)
    return gene_dict


def prepare_primer3_input(seq_id, subsequence, target, num_return, output_path, relax=False, min_tm=None, max_tm=None, max_product_size=None):
    """
    Prepares input file for Primer3 using user defined inputs on command line
    """
    primer3_input_path = output_path
    primer3_in = generate_primer3_input(
        seq_id, subsequence, target,
        num_return,
        relax=relax,
        min_tm=min_tm,
        max_tm=max_tm,
        max_product_size=max_product_size,
        save_path=primer3_input_path,
    )
    return primer3_input_path, primer3_in


def genbank_to_fasta(gb_file, tmp_dir):
    """
    Saves genbank as a temp fasta file for mapping of primers with bwa
    """
    fasta_path = os.path.join(tmp_dir, os.path.basename(gb_file).rsplit(".", 1)[0] + ".fasta")
    with open(fasta_path, "w", encoding="utf-8") as out_f:
        for rec in SeqIO.parse(gb_file, "genbank"):
            out_f.write(f">{rec.id}\n{str(rec.seq)}\n")
    return fasta_path


def get_target_region(records, args, gff_features=None):
    """
    Determine target region from args.
    Returns (seq_record, start, stop).
    Supports:
      - explicit coordinates + contig
      - gene name / locus_tag in GenBank features
      - gene name / locus_tag in GFF features
    """
    # Case 1: explicit coords
    if args.start_coord and args.stop_coord:
        if not args.contig:
            raise ValueError("--contig is required when using --start_coord/--stop_coord")
        if args.contig not in records:
            raise ValueError(f"Contig {args.contig} not found in genome")
        seq_record = records[args.contig]
        return seq_record, int(args.start_coord), int(args.stop_coord)

    # Case 2: gene / locus_tag
    if args.gene or args.locus_tag:
        matches = []

        # Check GenBank features first
        for rec in records.values():
            if getattr(rec, "features", None):
                for feature in rec.features:
                    if feature.type != "gene":
                        continue
                    q = {k.lower(): v for k, v in feature.qualifiers.items()}
                    if args.gene and "gene" in q and args.gene in q["gene"]:
                        matches.append((rec, feature))
                    if args.locus_tag and "locus_tag" in q and args.locus_tag in q["locus_tag"]:
                        return rec, int(feature.location.start) + 1, int(feature.location.end)

        # If nothing in GenBank, check GFF features
        if gff_features:
            for contig, features in gff_features.items():
                for feature in features:
                    quals = {k.lower(): v for k, v in getattr(feature, "qualifiers", {}).items()}
                    if args.gene and "gene" in quals and args.gene in quals["gene"]:
                        matches.append((records[contig], feature))
                    if args.locus_tag and "locus_tag" in quals and args.locus_tag in quals["locus_tag"]:
                        return records[contig], int(feature.location.start) + 1, int(feature.location.end)

        if args.gene:
            if len(matches) == 0:
                raise ValueError(f"No matches found for gene '{args.gene}'")
            elif len(matches) > 1:
                raise ValueError(
                    f"Multiple matches ({len(matches)}) found for gene '{args.gene}'. "
                    f"Use --locus_tag for unambiguous targeting."
                )
            rec, feature = matches[0]
            return rec, int(feature.location.start) + 1, int(feature.location.end)

        if args.locus_tag:
            raise ValueError(f"Locus_tag {args.locus_tag} not found in genome / annotation")

    raise ValueError("Must supply either (start/stop+contig) OR --gene OR --locus_tag.")


def map_primer_to_genome(primer_pair, genome_path, tmp_dir, bwa_cmd):
    """
    Map primers to genome with bwa for in silico PCR
    """
    fasta_for_bwa = genbank_to_fasta(genome_path, tmp_dir) if genome_path.endswith((".gbk", ".gbff")) else genome_path
    amplicons = perform_in_silico_pcr_bwa(fasta_for_bwa, primer_pair["forward"], primer_pair["reverse"], bwa_cmd=bwa_cmd)
    for a in amplicons:
        a["contig"] = a.get("contig", "")
    return amplicons


def run_in_silico_pcr_on_design(results, genome_list_dict, tmp_dir, bwa_cmd, output_dir):
    """
    Run in-silico PCR on designed primers across included/excluded genomes
    and write one TSV per primer pair.
    """
    os.makedirs(output_dir, exist_ok=True)

    for r in results:
        pair_num = r.get('pair', '')
        pair_name = f"pair{pair_num}"
        all_amplicons_dict = {}

        for genome_path in genome_list_dict:
            # genome_label = genome_list_dict[genome_path]
            amplicons = map_primer_to_genome(r, genome_path, tmp_dir, bwa_cmd)
            forward_name = f"{pair_name}_F"
            reverse_name = f"{pair_name}_R"

            all_amplicons_dict[os.path.basename(genome_path).rsplit('.', 1)[0]] = [{
                "forward_name": forward_name,
                "forward_seq": r.get("forward_seq") or r.get("forward") or "",
                "reverse_name": reverse_name,
                "reverse_seq": r.get("reverse_seq") or r.get("reverse") or "",
                "amplicons": amplicons
            }]

        out_file = os.path.join(output_dir, f"{pair_name}_in_silico_pcr.tsv")
        write_amplicons_tsv(all_amplicons_dict, out_file)
        print(f"Wrote in-silico PCR results for {pair_name} -> {out_file}")


def filter_and_validate_primers(results, genome_list_dict, ref_genome, tmp_dir, bwa_cmd, threads, max_diff=200):
    """
    Filter primers for those that have a reasonably similar amplicon length, and if they are present in the "include"
    genomes while being absent in the "exclude" genomes
    """
    validated_results = []
    included = [g for g, l in genome_list_dict.items() if l == "include"]
    excluded = [g for g, l in genome_list_dict.items() if l == "exclude"]

    for r in results:
        include_ok = True
        exclude_fail = False
        genome_amplicons = {}

        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            futures = {executor.submit(map_primer_to_genome, r, g, tmp_dir, bwa_cmd): g for g in genome_list_dict}
            print(f"Validating primer pair {r.get('pair', '')}...")
            for fut in concurrent.futures.as_completed(futures):
                genome = futures[fut]
                amplicons = fut.result()
                print(f"  Mapped to {genome}: {len(amplicons)} amplicons found")
                genome_amplicons[genome] = amplicons
                label = genome_list_dict[genome]
                if label == "include" and not amplicons:
                    include_ok = False
                if label == "exclude" and amplicons:
                    exclude_fail = True

        if include_ok and not exclude_fail:
            # Check amplicon lengths across included genomes
            included_amplicons = [genome_amplicons[g] for g in included if genome_amplicons[g]]
            included_lengths = [a[0]["amplicon_length"] for a in included_amplicons if a]
            if included_lengths and max(included_lengths) - min(included_lengths) > max_diff:
                print(f"  Skipping primer {r.get('pair', '')} due to amplicon length variation: {max(included_lengths) - min(included_lengths)} bp")
                continue

            prop_incl = sum(bool(genome_amplicons[g]) for g in included) / max(len(included), 1)
            prop_excl = sum(bool(genome_amplicons[g]) for g in excluded) / max(len(excluded), 1)

            r["amplicons"] = genome_amplicons[ref_genome]
            r["prop_included"] = prop_incl
            r["prop_excluded"] = prop_excl
            validated_results.append(r)

    return validated_results


def extract_subsequence(seq_record, start, stop, flank=500):
    """Return subsequence with flanking regions and target coordinates relative to subsequence."""
    contig_len = len(seq_record.seq)
    subseq_start = max(1, start - flank)
    subseq_end = min(contig_len, stop + flank)
    subsequence = str(seq_record.seq[subseq_start - 1:subseq_end])
    target_start_rel = start - subseq_start
    target_length = stop - start + 1
    target = f"{target_start_rel},{target_length}"
    return subsequence, target


def main(args=None):
    """
    Main for designing PCR primers
    """

    if args is None:
        parser = argparse.ArgumentParser(description="Design primers with Primer3")
        parser.add_argument("--input_genome")
        parser.add_argument("--genome_list")
        parser.add_argument("--gff")
        parser.add_argument("--output", required=True)
        parser.add_argument("--primer3", default="primer3_core")
        parser.add_argument("--bwa", default="bwa")
        parser.add_argument("--relax", action="store_true")
        parser.add_argument("--num_return", type=int, default=20)
        parser.add_argument("--min_tm", type=float)
        parser.add_argument("--max_tm", type=float)
        parser.add_argument("--start_coord", type=int)
        parser.add_argument("--stop_coord", type=int)
        parser.add_argument("--gene")
        parser.add_argument("--locus_tag")
        parser.add_argument("--contig")
        parser.add_argument("--max_amplicon_diff", type=int, default=500,
                            help="Maximum allowed difference in amplicon lengths across included genomes (default: 500 bp)")
        parser.add_argument("--threads", type=int, default=multiprocessing.cpu_count())
        args = parser.parse_args()

    out_dir = os.path.dirname(args.output) or "."
    tmp_dir = os.path.join(out_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    # Determine reference genome & FASTA
    ref_genome, ref_fasta, genome_list_dict = resolve_reference_genome(args, tmp_dir)

    # Load reference genome records
    fmt = detect_format(ref_genome)
    records = {rec.id.split()[0]: rec for rec in SeqIO.parse(ref_genome, fmt)}

    # Parse GFF if provided
    gff_features = parse_gff(args.gff) if args.gff else None

    seq_record, start, stop = get_target_region(records, args, gff_features=gff_features)
    print(f"Targeting {seq_record.id}:{start}-{stop} in {ref_genome}")

    # Extract subsequence ±200bp
    subsequence, target = extract_subsequence(seq_record, start, stop, flank=500)

    _, primer3_in = prepare_primer3_input(
        seq_record.id, subsequence, target,
        num_return=args.num_return * 2,
        output_path=os.path.join(tmp_dir, os.path.basename(args.output) + "_primer3_input.txt"),
        relax=args.relax,
        min_tm=getattr(args, "min_tm", None),
        max_tm=getattr(args, "max_tm", None),
        max_product_size=(stop - start + 1) + 2 * 200
    )

    # Run Primer3
    with tempfile.NamedTemporaryFile("w", delete=False) as temp_in:
        temp_in.write(primer3_in)
        temp_in_path = temp_in.name
    try:
        with open(temp_in_path, "r", encoding="utf-8") as temp_in_file:
            # pylint: disable=subprocess-run-check
            result = subprocess.run([args.primer3],
                                    stdin=temp_in_file,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    universal_newlines=True)
    finally:
        os.remove(temp_in_path)
    if result.returncode != 0:
        print("Primer3 failed:\n", result.stderr)
        return

    results = parse_primer3_output_text(result.stdout)
    print(f"Primer3 returned {len(results)} primer pairs.")

    # Map and validate primers across genomes
    if genome_list_dict:
        results = filter_and_validate_primers(
            results, genome_list_dict, ref_fasta, tmp_dir, args.bwa, args.threads,
            max_diff=args.max_amplicon_diff
        )
        results = results[:args.num_return]
    else:
        for r in results:
            r["amplicons"] = map_primer_to_genome(r, ref_genome, tmp_dir, args.bwa)
            r["prop_included"] = None
            r["prop_excluded"] = None
            # store list of genomes checked to ensure even missing hits are written
            r["genomes_checked"] = [ref_genome]

    # Write in silico PCR TSVs per primer pair
    in_silico_outdir = os.path.join(out_dir, "in_silico_pcr")
    run_in_silico_pcr_on_design(results, genome_list_dict, tmp_dir, args.bwa, in_silico_outdir)

    # Write outputs
    write_primer_results_tsv(results, f"{args.output}.tsv", include_props=True)
    write_amplicon_sequences_to_fasta(results, ref_fasta, tmp_dir)
    write_primer_fasta_pairs(results, args.output)

    print(f"Done!\nTSV: {args.output}.tsv\nPrimer FASTAs: {args.output}_pairX.fasta")
    print(f"Temporary files: {tmp_dir}")


if __name__ == "__main__":
    main()
