#!/usr/bin/env python3
"""
in_silico_pcr_utils.py

Utilities script for in silico PCR. Supports design.py and in_silico_pcr.py
with sliding-window logic for optimal primer selection.
"""

import os
import subprocess
import tempfile
from Bio import SeqIO


def perform_in_silico_pcr_bwa(fasta_reference, forward_primer, reverse_primer, bwa_cmd="bwa"):
    """
    Run in silico PCR using BWA-backtrack (aln + samse).
    Returns list of amplicons with:
        contig, amplicon_start, amplicon_end, amplicon_length
    """
    amplicons = []

    # Index reference if needed
    index_files = [fasta_reference + ext for ext in [".bwt", ".pac", ".ann", ".amb", ".sa"]]
    if not all(os.path.exists(f) for f in index_files):
        print(f"Indexing reference {fasta_reference} with BWA...")
        subprocess.run([bwa_cmd, "index", fasta_reference], check=True)

    # Temporary directory for primers
    with tempfile.TemporaryDirectory() as tmp_dir:
        fwd_path = os.path.join(tmp_dir, "primer_F.fasta")
        rev_path = os.path.join(tmp_dir, "primer_R.fasta")
        with open(fwd_path, "w", encoding="utf-8") as f:
            f.write(f">F\n{forward_primer}\n")
        with open(rev_path, "w", encoding="utf-8") as f:
            f.write(f">R\n{reverse_primer}\n")

        # Align primers
        fwd_sam = _bwa_aln_sam(fasta_reference, fwd_path, tmp_dir, bwa_cmd, "F")
        rev_sam = _bwa_aln_sam(fasta_reference, rev_path, tmp_dir, bwa_cmd, "R")

        fwd_hits = _parse_sam_hits(fwd_sam)
        rev_hits = _parse_sam_hits(rev_sam)

        # Match hits: one + and one - on same contig
        for f_hit in fwd_hits:
            for r_hit in rev_hits:
                if f_hit["contig"] != r_hit["contig"]:
                    continue
                if {f_hit["strand"], r_hit["strand"]} == {"+", "-"}:
                    start = min(f_hit["start"], r_hit["start"])
                    end = max(f_hit["start"] + len(f_hit["seq"]) - 1,
                              r_hit["start"] + len(r_hit["seq"]) - 1)
                    amplicons.append({
                        "contig": f_hit["contig"],
                        "amplicon_start": start,
                        "amplicon_end": end,
                        "amplicon_length": end - start + 1
                    })
    return amplicons


def _bwa_aln_sam(fasta_reference, primer_fasta, tmp_dir, bwa_cmd, prefix):
    """Run BWA aln + samse for a single primer and return SAM path."""
    sai_file = os.path.join(tmp_dir, f"{prefix}.sai")
    sam_file = os.path.join(tmp_dir, f"{prefix}.sam")
    subprocess.run([bwa_cmd, "aln", fasta_reference, primer_fasta, "-f", sai_file],
                   check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run([bwa_cmd, "samse", fasta_reference, sai_file, primer_fasta, "-f", sam_file],
                   check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return sam_file


def _parse_sam_hits(sam_file):
    """Parse SAM file to extract primer hits with contig, start, strand, and seq."""
    hits = []
    with open(sam_file, encoding="utf-8") as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.split("\t")
            if len(fields) < 11 or fields[2] == "*":
                continue
            contig = fields[2].split()[0]
            pos = int(fields[3])
            strand = "-" if int(fields[1]) & 16 else "+"
            hits.append({"contig": contig, "start": pos, "strand": strand, "seq": fields[9]})
    return hits


def safe_str(x):
    """Convert None to empty string, otherwise str(x)."""
    return "" if x is None else str(x)


def write_amplicons_tsv(all_amplicons_dict, output_file):
    """
    Write in silico PCR results to TSV.
    """
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    header = [
        "Genome",
        "Contig",
        "Forward_primer_name",
        "Forward_primer_seq",
        "Reverse_primer_name",
        "Reverse_primer_seq",
        "Amplicon_start_in_silico",
        "Amplicon_end_in_silico",
        "Amplicon_length_in_silico",
        "Amplicon_count_in_silico"
    ]

    with open(output_file, "w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for genome_name, primer_dicts in all_amplicons_dict.items():
            for p in primer_dicts:
                amplicons = p.get("amplicons", [])
                count = len(amplicons)
                if count == 0:
                    row = [
                        str(genome_name),
                        "",
                        str(p.get("forward_name") or ""),
                        str(p.get("forward_seq") or ""),
                        str(p.get("reverse_name") or ""),
                        str(p.get("reverse_seq") or ""),
                        "", "", "", "0"
                    ]
                    f.write("\t".join(row) + "\n")
                else:
                    for a in amplicons:
                        row = [
                            str(genome_name),
                            str(a.get("contig") or ""),
                            str(p.get("forward_name") or ""),
                            str(p.get("forward_seq") or ""),
                            str(p.get("reverse_name") or ""),
                            str(p.get("reverse_seq") or ""),
                            str(a.get("amplicon_start") or ""),
                            str(a.get("amplicon_end") or ""),
                            str(a.get("amplicon_length") or ""),
                            str(count)
                        ]
                        f.write("\t".join(row) + "\n")


def write_amplicon_sequences_to_fasta(results, reference_fasta, output_dir):
    """
    Write predicted amplicon sequences to FASTA files with robust contig handling.
    """

    os.makedirs(output_dir, exist_ok=True)
    ref_dict = {rec.id.split()[0]: rec for rec in SeqIO.parse(reference_fasta, "fasta")}
    print(f"[write_amplicon_sequences_to_fasta] Loaded {len(ref_dict)} contigs from reference.")

    for r in results:
        amplicons = r.get("amplicons", [])
        if not amplicons:
            print(f"[write_amplicon_sequences_to_fasta] Primer pair {r['pair']} has no amplicons. Skipping.")
            continue

        fasta_path = os.path.join(output_dir, f"pair{r['pair']}_amplicons.fasta")
        with open(fasta_path, "w", encoding="utf-8") as f:
            for i, a in enumerate(amplicons):
                norm_contig = a["contig"].strip()
                # remove trailing .fasta if present
                if norm_contig.endswith(".fasta"):
                    norm_contig = norm_contig.rsplit(".", 1)[0]

                if norm_contig not in ref_dict:
                    print(f"[WARNING] Contig '{norm_contig}' not found in reference. Skipping amplicon {i}.")
                    continue

                seq = ref_dict[norm_contig].seq[a["amplicon_start"] - 1: a["amplicon_end"]]
                f.write(f">pair{r['pair']}_amplicon{i}_{norm_contig}:{a['amplicon_start']}-{a['amplicon_end']}\n")
                f.write(str(seq) + "\n")
        print(f"[write_amplicon_sequences_to_fasta] Wrote {len(amplicons)} amplicons to {fasta_path}")
