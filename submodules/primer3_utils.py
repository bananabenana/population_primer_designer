#!/usr/bin/env python3
"""
primer3_utils.py

Utilities script for primer3. Assists with design.py.
"""

import re
from Bio import SeqIO


def generate_primer3_input(seq_id, sequence, target=None, num_return=5, relax=False,
                           min_tm=None, max_tm=None, max_product_size=None, save_path=None):
    """Create Primer3 input text for a given sequence and target region."""
    primer3_input = [
        f"SEQUENCE_ID={seq_id}",
        f"SEQUENCE_TEMPLATE={sequence}",
    ]
    if target:
        primer3_input.append(f"SEQUENCE_TARGET={target}")
    primer3_input.extend([
        "PRIMER_TASK=pick_pcr_primers",
        "PRIMER_PICK_LEFT_PRIMER=1",
        "PRIMER_PICK_RIGHT_PRIMER=1",
        "PRIMER_EXPLAIN_FLAG=1",
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1",
        "PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=1",
        f"PRIMER_NUM_RETURN={num_return}",
    ])
    if relax:
        primer3_input.extend([
            "PRIMER_MIN_TM=55",
            "PRIMER_MAX_TM=65",
            "PRIMER_MIN_SIZE=18",
            "PRIMER_MAX_SIZE=25",
            "PRIMER_MAX_SELF_ANY=8",
        ])
    if min_tm is not None:
        primer3_input.append(f"PRIMER_MIN_TM={min_tm}")
    if max_tm is not None:
        primer3_input.append(f"PRIMER_MAX_TM={max_tm}")
    if max_product_size is not None and target:
        region_length = int(target.split(",")[1])
        primer3_input.append(f"PRIMER_PRODUCT_SIZE_RANGE={region_length}-{max_product_size}")
    primer3_input.append("=")

    content = "\n".join(primer3_input)
    if save_path:
        with open(save_path, "w", encoding="utf-8") as f:
            f.write(content)
        print(f"Saved Primer3 input to {save_path}")
    return content


def parse_primer3_output_text(output_text):
    """Parse Primer3 output text into structured primer pair dicts."""
    results = []
    i = 0
    while True:
        left = re.search(f"PRIMER_LEFT_{i}_SEQUENCE=(.+)", output_text)
        right = re.search(f"PRIMER_RIGHT_{i}_SEQUENCE=(.+)", output_text)
        if not left or not right:
            break

        def grab(pattern, cast=float):
            m = re.search(pattern, output_text)
            return cast(m.group(1)) if m else None
        results.append({
            "pair": i,
            "forward": left.group(1),
            "reverse": right.group(1),
            "forward_tm": grab(f"PRIMER_LEFT_{i}_TM=(.+)"),
            "reverse_tm": grab(f"PRIMER_RIGHT_{i}_TM=(.+)"),
            "product_size": grab(f"PRIMER_PAIR_{i}_PRODUCT_SIZE=(.+)", int),
            "forward_self_any": grab(f"PRIMER_LEFT_{i}_SELF_ANY_TH=(.+)"),
            "reverse_self_any": grab(f"PRIMER_RIGHT_{i}_SELF_ANY_TH=(.+)"),
            "cross_dimer": grab(f"PRIMER_PAIR_{i}_COMPL_ANY_TH=(.+)")
        })
        i += 1
    return results


def calculate_gc_percent(seq):
    seq = seq.upper()
    gc_count = seq.count("G") + seq.count("C")
    return round(100 * gc_count / len(seq), 2) if len(seq) > 0 else 0


def write_primer_results_tsv(results, tsv_file, include_props=False):
    header = [
        "Primer_pair", "Forward_primer", "Forward_Tm", "Forward_GC%",
        "Reverse_primer", "Reverse_Tm", "Reverse_GC%",
        "Product_Size", "Forward_self_dimer_temperature",
        "Reverse_self_dimer_temperature", "Forward_reverse_dimer_temperature",
        "Contig_amplified_from", "Amplicon_start_in_silico",
        "Amplicon_end_in_silico", "Amplicon_length_in_silico",
        "Amplicon_count_in_silico"
    ]
    if include_props:
        header.extend(["Proportion_included_genomes_with_amplicon", "Proportion_excluded_genomes_with_amplicon"])

    with open(tsv_file, "w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in results:
            amplicons = r.get("amplicons", [])
            if amplicons:
                contigs = [a.get("contig", "") for a in amplicons]
                starts = [str(a.get("amplicon_start", "")) for a in amplicons]
                ends = [str(a.get("amplicon_end", "")) for a in amplicons]
                lengths = [str(a.get("amplicon_length", "")) for a in amplicons]
                contig_str = ";".join(contigs)
                start_str = ";".join(starts)
                end_str = ";".join(ends)
                length_str = ";".join(lengths)
                count = len(amplicons)
            else:
                contig_str = start_str = end_str = length_str = ""
                count = 0

            row = [
                str(r.get("pair", "")),
                r.get("forward", ""), str(r.get("forward_tm", "")), str(calculate_gc_percent(r.get("forward", ""))),
                r.get("reverse", ""), str(r.get("reverse_tm", "")), str(calculate_gc_percent(r.get("reverse", ""))),
                str(r.get("product_size", "")),
                str(r.get("forward_self_any", "")),
                str(r.get("reverse_self_any", "")),
                str(r.get("cross_dimer", "")),
                contig_str, start_str, end_str, length_str, str(count)
            ]

            if include_props:
                prop_incl = r.get("prop_included")
                prop_excl = r.get("prop_excluded")
                row.extend([
                    str(round(prop_incl, 3)) if prop_incl is not None else "",
                    str(round(prop_excl, 3)) if prop_excl is not None else ""
                ])
            f.write("\t".join(row) + "\n")


def write_primer_fasta_pairs(results, output_prefix):
    """Write primer sequences to separate FASTA files for each pair."""
    for r in results:
        fasta_file = f"{output_prefix}_pair{r['pair']}.fasta"
        with open(fasta_file, "w", encoding="utf-8") as f:
            f.write(f">pair{r['pair']}_F\n{r['forward']}\n")
            f.write(f">pair{r['pair']}_R\n{r['reverse']}\n")
