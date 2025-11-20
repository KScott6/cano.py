#!/usr/bin/env python3

import os
import math
import argparse
import time
import subprocess
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO
import pandas as pd
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Canopy phylogenomic pipeline")
    parser.add_argument("--omes", required=True, help="Path to file with OME prefixes (one OME per line)")
    parser.add_argument("--busco", required=True, help="Path to BUSCO output parent directory")
    parser.add_argument("--threads", "-t", type=int, default=4, help="Number of threads to use")
    parser.add_argument(
        "--trimal_strategy",
        default="automated1",
        choices=["automated1", "gappyout", "strict", "strictplus"],
        help="TrimAl strategy"
    )
    parser.add_argument(
        "--bootstraps",
        "-B",
        type=int,
        default=1000,
        help="Number of IQ-TREE ULTRAFAST bootstraps for single-gene trees and final tree"
    )
    parser.add_argument(
        "--sco_threshold",
        type=float,
        default=1.0,
        help="Minimum proportion of genomes with a BUSCO to retain it as single-copy (default: 1.0)"
    )
    parser.add_argument("--output", "-o", required=True, help="Output folder")
    parser.add_argument(
        "--continue",
        dest="continue_run",
        action="store_true",
        help="Continue from an existing run by skipping steps where outputs already exist"
    )
    return parser.parse_args()

def find_busco_sc_path(ome, busco_base):
    p1 = Path(busco_base) / ome / ome
    p2 = Path(busco_base) / ome
    for base in [p1, p2]:
        if base.exists():
            for sub in base.glob("**/single_copy_busco_sequences"):
                return sub
    return None

def run_muscle(args):
    input_fasta, output_aln = args
    cmd = ["muscle", "-align", input_fasta, "-output", output_aln, "-threads", "1"]
    try:
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    except subprocess.CalledProcessError:
        return Path(input_fasta).stem, False
    return Path(input_fasta).stem, True

def run_trimal(args):
    input_aln, output_trimmed, strategy = args
    cmd = ["trimal", "-in", input_aln, "-out", output_trimmed, f"-{strategy}"]
    try:
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    except subprocess.CalledProcessError:
        return Path(input_aln).stem, False
    return Path(input_aln).stem, True

def run_iqtree(aln_file, output_dir, bootstraps=None):
    """
    Run IQ-TREE v3 on a single-gene alignment.

    aln_file   : Path to trimmed alignment
    output_dir : Directory for IQ-TREE outputs
    bootstraps : If not None and >0, run ultrafast bootstraps (-B)
    """
    busco_id = aln_file.stem.replace(".trimmed", "")
    prefix_path = output_dir / busco_id

    cmd = [
        "iqtree3",
        "-s", str(aln_file),
        "-m", "MFP+MERGE",
        "--prefix", str(prefix_path),
        "-nt", "AUTO",
    ]

    # Only add bootstraps if requested 
    if bootstraps is not None and int(bootstraps) > 0:
        cmd.extend(["-B", str(bootstraps)])

    try:
        subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )
    except subprocess.CalledProcessError:
        return None

    return prefix_path.with_suffix(".iqtree")

def extract_model_info(iqtree_file):
    if not iqtree_file.exists():
        return None
    model, length = None, None
    busco_id = iqtree_file.stem
    with iqtree_file.open() as f:
        for line in f:
            if "Input data:" in line and "amino-acid sites" in line:
                match = re.search(r"(\d+) amino-acid sites", line)
                if match:
                    length = int(match.group(1))
            elif "Best-fit model according to BIC:" in line:
                model = line.strip().split(":")[-1].strip()
            if model and length:
                break
    return (busco_id, length, model) if model and length else None

def main():
    args = parse_args()
    start_time = time.time()

    output_dir = Path(args.output)
    sc_seq_dir = output_dir / "single_copy_sequences"
    aln_dir = output_dir / "alignments"
    trimmed_dir = output_dir / "trimmed_alignments"
    single_tree_dir = output_dir / "single_gene_trees"
    concat_dir = output_dir / "concatenated_alignment"
    final_tree_dir = output_dir / "final_tree"
    log_dir = output_dir / "logs"
    for d in [sc_seq_dir, aln_dir, trimmed_dir, single_tree_dir, concat_dir, final_tree_dir, log_dir]:
        d.mkdir(parents=True, exist_ok=True)

    with open(args.omes) as f:
        prefixes = [line.strip() for line in f if line.strip()]

    all_buscos, genome_buscos, busco_presence = set(), defaultdict(set), defaultdict(set)
    missing_dirs, genomes_with_no_buscos = [], []

    # Scan BUSCO outputs
    for ome in prefixes:
        sc_dir = find_busco_sc_path(ome, args.busco)
        if sc_dir is None:
            missing_dirs.append(ome)
            continue
        faa_files = list(sc_dir.glob("*.faa"))
        if not faa_files:
            genomes_with_no_buscos.append(ome)
            continue
        ids = [f.stem for f in faa_files]
        genome_buscos[ome] = set(ids)
        for b in ids:
            busco_presence[b].add(ome)
        all_buscos.update(ids)

    # Log missing BUSCO dirs and genomes with no BUSCOs
    with open(log_dir / "missing_busco_dirs.txt", "w") as f:
        for ome in missing_dirs:
            f.write(f"{ome}\n")
    with open(log_dir / "genomes_with_no_buscos.txt", "w") as f:
        for ome in genomes_with_no_buscos:
            f.write(f"{ome}\n")

    # Full presence/absence matrix
    full_matrix = []
    for ome in prefixes:
        row = {"genome": ome}
        for b in sorted(all_buscos):
            row[b] = 1 if b in genome_buscos.get(ome, set()) else 0
        full_matrix.append(row)
    pd.DataFrame(full_matrix).to_csv(
        output_dir / "busco_genes_presence_absence_matrix_FULL.tsv",
        sep="\t",
        index=False
    )

    # Retained single-copy BUSCOs (using math.ceil for threshold, otherwise weird rounding errors)
    min_required = math.ceil(len(prefixes) * args.sco_threshold)
    retained_buscos = sorted([b for b in all_buscos if len(busco_presence[b]) >= min_required])

    retained_matrix = []
    for ome in prefixes:
        row = {"genome": ome}
        for b in retained_buscos:
            row[b] = 1 if b in genome_buscos.get(ome, set()) else 0
        retained_matrix.append(row)
    pd.DataFrame(retained_matrix).to_csv(
        output_dir / "busco_genes_presence_absence_matrix.tsv",
        sep="\t",
        index=False
    )

    if not retained_buscos:
        with open(output_dir / "final_report.txt", "w") as f:
            f.write(
                "No BUSCOs met the single-copy threshold. Cannot proceed. "
                "Consider lowering the threshold with --sco_threshold.\n"
            )
        return

    # Extract single-copy sequences per BUSCO and genome
    busco_records = defaultdict(list)
    for ome in prefixes:
        sc_dir = find_busco_sc_path(ome, args.busco)
        if ome not in genome_buscos or sc_dir is None:
            continue
        for busco_id in retained_buscos:
            fasta_path = sc_dir / f"{busco_id}.faa"
            if fasta_path.exists():
                try:
                    record = SeqIO.read(fasta_path, "fasta")
                    record.id = ome
                    record.description = ""
                    busco_records[busco_id].append(record)
                except Exception:
                    continue

    # Write combined per-BUSCO FASTAs (safe to overwrite / re-write)
    for busco_id, records in busco_records.items():
        SeqIO.write(records, sc_seq_dir / f"{busco_id}.faa", "fasta")

    # MUSCLE alignments (checkpoint-aware)
    muscle_jobs = []
    for busco_id in busco_records:
        input_fasta = sc_seq_dir / f"{busco_id}.faa"
        output_aln = aln_dir / f"{busco_id}.aln"
        if args.continue_run and output_aln.exists():
            continue  # skip existing alignments
        muscle_jobs.append((str(input_fasta), str(output_aln)))

    if muscle_jobs:
        with Pool(processes=args.threads) as pool:
            pool.map(run_muscle, muscle_jobs)

    # TrimAl (checkpoint-aware)
    trimal_jobs = []
    for busco_id in busco_records:
        input_aln = aln_dir / f"{busco_id}.aln"
        output_trimmed = trimmed_dir / f"{busco_id}.trimmed.aln"
        if not input_aln.exists():
            continue  # nothing to trim
        if args.continue_run and output_trimmed.exists():
            continue  # skip existing trimmed alignments
        trimal_jobs.append((str(input_aln), str(output_trimmed), args.trimal_strategy))

    if trimal_jobs:
        with Pool(processes=args.threads) as pool:
            pool.map(run_trimal, trimal_jobs)

    # Single-gene IQ-TREE (checkpoint-aware, with bootstraps)
    aln_files = sorted(trimmed_dir.glob("*.trimmed.aln"))
    iqtree_jobs = []
    for aln_file in aln_files:
        busco_id = aln_file.stem.replace(".trimmed", "")
        expected_iqtree = single_tree_dir / f"{busco_id}.iqtree"
        if args.continue_run and expected_iqtree.exists():
            continue  # skip if IQ-TREE output already exists
        iqtree_jobs.append(aln_file)

    if iqtree_jobs:
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            futures = {
                executor.submit(run_iqtree, aln_file, single_tree_dir, args.bootstraps): aln_file
                for aln_file in iqtree_jobs
            }
            for future in as_completed(futures):
                _ = future.result()

    # Collect all .iqtree files (from previous and current runs)
    iqtree_files = sorted(single_tree_dir.glob("*.iqtree"))

    # Summarize per-gene models
    model_summary = []
    for iqtree_file in iqtree_files:
        info = extract_model_info(iqtree_file)
        if info:
            model_summary.append(info)
    model_df = pd.DataFrame(model_summary, columns=["busco_id", "length", "model"])
    model_df.to_csv(single_tree_dir / "partition_model_summary.tsv", sep="\t", index=False)

    # Build concatenated alignment and partition info
    genome_seqs = defaultdict(list)
    partition_lines = []
    current_start = 1
    gene_order = [f.stem.split(".")[0] for f in aln_files]

    for aln_file in aln_files:
        busco_id = aln_file.stem.split(".")[0]
        records = {rec.id: str(rec.seq) for rec in SeqIO.parse(aln_file, "fasta")}
        if not records:
            continue
        aln_length = len(next(iter(records.values())))
        for genome in prefixes:
            genome_seqs[genome].append(records.get(genome, "-" * aln_length))
        partition_lines.append(f"PROT, {busco_id} = {current_start}-{current_start + aln_length - 1}")
        current_start += aln_length

    # Write concatenated alignment
    with open(concat_dir / "concatenated_alignment.fasta", "w") as f:
        for genome, seqs in genome_seqs.items():
            f.write(f">{genome}\n{''.join(seqs)}\n")

    # Write partition coordinates
    with open(concat_dir / "partition_coordinates.txt", "w") as f:
        for line in partition_lines:
            f.write(f"{line}\n")

    model_dict = dict(zip(model_df["busco_id"].astype(str), model_df["model"]))
    with open(concat_dir / "partition.nexus", "w") as f:
        f.write("#nexus\nbegin sets;\n")
        for line in partition_lines:
            busco_id = line.split(",")[1].split("=")[0].strip()
            coords = line.split("=")[1].strip()
            f.write(f"    charset {busco_id} = {coords};\n")
        f.write("    charpartition mine = ")
        f.write(", ".join(f"{model_dict[busco_id]}:{busco_id}" for busco_id in gene_order if busco_id in model_dict))
        f.write(";\nend;\n")

    # Final concatenated IQ-TREE run (checkpoint-aware)
    final_treefile = final_tree_dir / "final_tree.treefile"
    if not (args.continue_run and final_treefile.exists()):
        final_cmd = [
            "iqtree3",
            "-s", str(concat_dir / "concatenated_alignment.fasta"),
            "-p", str(concat_dir / "partition.nexus"),
            "-T", str(args.threads),
            "--ufboot", str(args.bootstraps),
            "--prefix", str(final_tree_dir / "final_tree")
        ]
        subprocess.run(final_cmd, check=True)

    # Final report
    with open(output_dir / "final_report.txt", "w") as f:
        f.write("Final report generated by canopy phylogenomic pipeline\n")
        f.write("-----------------------------------------------\n")
        f.write(f"Prefix file used: {args.omes}\n")
        f.write(f"Number of genomes selected: {len(prefixes)}\n")
        f.write(f"Number of BUSCOs retained: {len(retained_buscos)}\n")
        f.write(f"SCO threshold used: {args.sco_threshold} ({min_required} genomes required per BUSCO)\n")
        f.write(f"Bootstraps used: {args.bootstraps}\n")
        f.write(f"Continue mode: {args.continue_run}\n")
        f.write(f"Runtime (seconds): {round(time.time() - start_time, 2)}\n")

if __name__ == "__main__":
    main()
