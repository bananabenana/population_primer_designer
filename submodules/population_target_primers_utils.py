#!/usr/bin/env python3
"""
population_target_primers_utils.py

Tools and functions for population_target_primers.py
"""
import subprocess
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def shred_genome(fasta_path, shred_length, step_size):
    """Shred genome into fixed-length overlapping fragments."""
    shreds = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq = str(rec.seq)
        for start in range(0, len(seq) - shred_length + 1, step_size):
            end = start + shred_length
            frag = seq[start:end]
            shreds.append((rec.id, start, end, frag, fasta_path))
    print(f"[INFO] Shredded {fasta_path} into {len(shreds)} fragments")
    return shreds


def shred_genome_parallel(genome_paths, shred_length, step_size, threads):
    """Shred multiple genomes in parallel."""
    shredded_genomes = []

    with ThreadPoolExecutor(max_workers=threads) as executor:
        # Submit shredding tasks
        future_to_genome = {
            executor.submit(shred_genome, genome_path, shred_length, step_size): genome_path
            for genome_path in genome_paths
        }

        # Collect results as they complete
        for future in as_completed(future_to_genome):
            genome_path = future_to_genome[future]
            try:
                shreds = future.result()
                shredded_genomes.extend(shreds)
            except Exception as exc:
                print(f'[ERROR] Genome {genome_path} generated an exception: {exc}')
                raise

    return shredded_genomes


def run_minimap2_single(ref_genome, shreds_fasta, tmp_dir, threads_per_job):
    """Run minimap2 for a single reference genome."""
    paf_file = os.path.join(tmp_dir, f"{os.path.basename(ref_genome)}.paf")

    cmd = ["minimap2", "-x", "sr", "-t", str(threads_per_job), ref_genome, shreds_fasta, "-o", paf_file]

    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

        # Parse results
        hits = {}
        with open(paf_file, encoding='utf-8') as f:
            for line in f:
                fields = line.strip().split("\t")
                if not fields:
                    continue
                qname = fields[0]
                hits[qname] = hits.get(qname, 0) + 1

        return ref_genome, hits
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] minimap2 failed for {ref_genome}: {e}")
        print(f"[ERROR] Command: {' '.join(cmd)}")
        print(f"[ERROR] stderr: {e.stderr}")
        raise


def map_shreds_minimap2(shreds, reference_genomes, threads, tmp_dir):
    """Map shredded fragments to reference genome(s) using minimap2 (PAF output) - parallelized."""
    fasta_path = os.path.join(tmp_dir, "shreds.fasta")

    # Include genome name in the FASTA ID for traceability
    SeqIO.write(
        [SeqRecord(Seq(seq), id=f"{contig}_{start}_{end}_{os.path.basename(genome_path)}")
         for contig, start, end, seq, genome_path in shreds],
        fasta_path, "fasta"
    )

    paf_results = {}

    # For single reference genome, just run it with all threads
    if len(reference_genomes) == 1:
        ref_genome, hits = run_minimap2_single(reference_genomes[0], fasta_path, tmp_dir, threads)
        paf_results[ref_genome] = hits
    else:
        # For multiple reference genomes, parallelize across them
        threads_per_job = max(1, threads // len(reference_genomes))

        with ThreadPoolExecutor(max_workers=len(reference_genomes)) as executor:
            future_to_ref = {
                executor.submit(run_minimap2_single, ref, fasta_path, tmp_dir, threads_per_job): ref
                for ref in reference_genomes
            }

            for future in as_completed(future_to_ref):
                ref_genome = future_to_ref[future]
                try:
                    ref_genome, hits = future.result()
                    paf_results[ref_genome] = hits
                except Exception as exc:
                    print(f'[ERROR] Reference {ref_genome} generated an exception: {exc}')
                    raise

    return paf_results


def parse_paf_file(paf_file):
    """
    Parse a PAF file and return a list of hits as dictionaries.

    Each dictionary contains:
        - query_name
        - query_length
        - query_start
        - query_end
        - strand
        - target_name
        - target_length
        - target_start
        - target_end
        - num_matching_bases
        - alignment_block_length
        - mapping_quality
    """
    hits = []
    with open(paf_file, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue  # ignore malformed lines

            hit = {
                "query_name": parts[0],
                "query_length": int(parts[1]),
                "query_start": int(parts[2]),
                "query_end": int(parts[3]),
                "strand": parts[4],
                "target_name": parts[5],
                "target_length": int(parts[6]),
                "target_start": int(parts[7]),
                "target_end": int(parts[8]),
                "num_matching_bases": int(parts[9]),
                "alignment_block_length": int(parts[10]),
                "mapping_quality": int(parts[11])
            }
            hits.append(hit)
    return hits


def compute_candidate_regions(ref_seq_records, paf_results, reference_genome,
                              include_genomes, exclude_genomes, shredded_genomes,
                              shred_length=1000, min_include_prop=0.95, min_contig_size=1000):
    """
    Compute candidate regions from reference genome and mapping results.
    Each shred-length window is treated independently.
    """

    depth_table = []
    detailed_table = []

    # Build mapping from shred_id -> genome
    shred_genome_map = {f"{contig}_{start}_{end}_{os.path.basename(genome_path)}": genome_path
                        for contig, start, end, seq, genome_path in shredded_genomes}

    ref_hits = paf_results[reference_genome]

    for rec in ref_seq_records:
        contig = rec.id
        contig_len = len(rec.seq)
        if contig_len < min_contig_size:
            continue

        for start in range(0, contig_len, shred_length):
            end = min(start + shred_length, contig_len)

            # Find shreds mapping to this window
            shreds_in_window = [sid for sid in ref_hits
                                if sid.startswith(f"{contig}_") and int(sid.split("_")[1]) < end and int(sid.split("_")[2]) > start]

            # Map to originating genome
            genome_in_window = [shred_genome_map[sid] for sid in shreds_in_window if sid in shred_genome_map]

            include_flags = [1 if g in genome_in_window else 0 for g in include_genomes]
            exclude_flags = [1 if g in genome_in_window else 0 for g in exclude_genomes]

            incl_prop = sum(include_flags) / len(include_flags) if include_flags else 0
            excl_prop = sum(exclude_flags) / len(exclude_flags) if exclude_flags else 0

            # Detailed table (always)
            detailed_table.append({
                "contig": contig,
                "start": start + 1,
                "end": end,
                **{g: include_flags[i] for i, g in enumerate(include_genomes)},
                **{g: exclude_flags[i] for i, g in enumerate(exclude_genomes)}
            })

            # Depth table: only windows meeting min_include_prop
            if incl_prop >= min_include_prop:
                depth_table.append({
                    "reference_genome": reference_genome,
                    "contig": contig,
                    "contig_length": contig_len,
                    "start": start + 1,
                    "end": end,
                    "region_size": end - start,  # Add region size for sorting
                    "include_genome_proportion_present": incl_prop,
                    "exclude_genome_proportion_present": excl_prop
                })

    # Sort depth_table: highest include, lowest exclude, largest size
    depth_table.sort(key=lambda x: (
        -x["include_genome_proportion_present"],
        x["exclude_genome_proportion_present"],
        -x["region_size"]
    ))

    return depth_table, detailed_table


def write_fasta(record, out_path):
    """
    Writes a fasta file output
    """
    with open(out_path, "w", encoding="utf-8") as f:
        SeqIO.write(record, f, "fasta")
