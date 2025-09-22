#!/usr/bin/env python3
"""
Search a FASTA for occurrences of each k-mer window in a query and
evaluate local duplex stability (Tm) around the match, allowing a single
bulge on either side with a fixed penalty.

Notes:
- The scanning logic is unchanged from the original script; this is a
  readability refactor that extracts helpers and adds type hints.
- Output for '+' strand hits: Tm, contig, start, end, query_align, db_align
- Tm threshold can be configured with --min-tm (default: 30.0).
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Iterable, List, Tuple
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

# Load NN Tm model (original behavior preserved)
import types

_NN_PATH = 'Santalucia_NN_Tm.py'
_SRC = open(_NN_PATH, 'r').read().replace('xrange', 'range')
_MOD = types.ModuleType('Santalucia_NN_Tm_py3')
exec(compile(_SRC, _NN_PATH, 'exec'), _MOD.__dict__)
NN_Tm = _MOD.NN_Tm

# Thermodynamic parameters (kept as in original)
primer_conc = 100
Na = 0
K = 50
Tris = 10
Mg = 1.5
dNTPs = 0.2


# ------------------------------
# Basic utilities
# ------------------------------


def read_fasta(path: str) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    header: str | None = None
    seq_chunks: List[str] = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_chunks).upper()))
                header = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            records.append((header, "".join(seq_chunks).upper()))
    return records


def sliding_kmers(seq: str, k: int) -> List[str]:
    seq = seq.upper()
    if k <= 0 or k > len(seq):
        return []
    return [seq[i : i + k] for i in range(len(seq) - k + 1)]


def revcomp(seq: str) -> str:
    """Return reverse complement (IUPAC-aware) of ``seq``."""
    comp = str.maketrans(
        "ACGTRYMKBDHVNacgtrymkbdhvn",
        "TGCAYRKMVHDBNtgcayrkmvhdbn",
    )
    return seq.translate(comp)[::-1]


def compl(seq: str) -> str:
    """Return complement (IUPAC-aware) of ``seq`` (no reverse)."""
    return seq.translate(
        str.maketrans(
            "ACGTRYMKBDHVNacgtrymkbdhvn",
            "TGCAYRKMVHDBNtgcayrkmvhdbn",
        )
    )


def find_all(haystack: str, needle: str) -> Iterable[int]:
    """Yield all 0-based start positions (including overlaps)."""
    start = 0
    while True:
        idx = haystack.find(needle, start)
        if idx == -1:
            break
        yield idx
        start = idx + 1


# ------------------------------
# Flank optimization
# ------------------------------


@dataclass
class FlankResult:
    """Optimization result for one flank side.

    All fields mirror the original variables: ``*_seq1/_seq2`` are the
    strings used for Tm calculation; ``*_ori`` preserve alignment markup
    where a gap is represented by a single '-' in the opposite string.
    ``seq2_len`` is used to expand the DB coordinates.
    """

    seq1: str = ""
    seq2: str = ""
    seq1_ori: str = ""
    seq2_ori: str = ""
    seq2_len: int = 0
    penalty: int = 0
    buldge_pos: int = -1
    buldge_seq: int = -1
    desc: str = ""


def calc_tm(seq: str, comp: str) -> float:
    """Wrapper for NN_Tm with configured salts and primer concentration."""
    return float(
        NN_Tm(
            seq=seq,
            compl_seq=comp,
            primer_conc=primer_conc,
            Na=Na,
            K=K,
            Tris=Tris,
            Mg=Mg,
            dNTPs=dNTPs,
            ion_corr=True,
        )
    )


def optimize_left_flank(query_left: str, left_flank: str, kmer_u: str) -> FlankResult:
    """Find best left-flank alignment per original heuristic.

    Baseline uses ``left_flank[1:]``; then try one-base bulge in DB or query
    with a fixed 10°C penalty.
    """
    if not (query_left and left_flank[1:]):
        return FlankResult()

    best_tm = calc_tm(query_left + kmer_u, compl(left_flank[1:] + kmer_u))
    best = FlankResult(
        seq1=query_left,
        seq2=left_flank[1:],
        seq1_ori=query_left,
        seq2_ori=left_flank[1:],
        seq2_len=len(left_flank[1:]),
        penalty=0,
        desc=f"max_pos=0 tm_C={best_tm} seq1={query_left} seq2={left_flank[1:]}",
    )

    # DB bulge (delete one base from DB flank)
    seq1 = query_left
    seq2 = left_flank
    for i in range(len(seq2)):
        seq2_del = seq2[:i] + seq2[i + 1 :]
        if seq1 and seq2_del:
            tm_del = calc_tm(seq1 + kmer_u, compl(seq2_del + kmer_u)) - 10
            if tm_del > best_tm:
                best_tm = tm_del
                best = FlankResult(
                    seq1=seq1,
                    seq2=seq2_del,
                    seq1_ori=seq1[:i] + "-" + seq1[i:],
                    seq2_ori=seq2,
                    seq2_len=len(seq2_del) + 1,
                    penalty=10,
                    buldge_pos=i,
                    buldge_seq=1,
                    desc=f"del_pos={i+1:02d} base={seq2[i]} tm_C={tm_del} seq1={seq1[:i]+'-'+seq1[i:]} seq2={seq2}",
                )

    # Query bulge (delete one base from query flank) using left_flank[2:]
    seq1 = query_left
    seq2_del = left_flank[2:]
    for i in range(len(seq1)):
        seq1_del = seq1[:i] + seq1[i + 1 :]
        if seq1_del and seq2_del:
            tm_del = calc_tm(seq1_del + kmer_u, compl(seq2_del + kmer_u)) - 10
            if tm_del > best_tm:
                best_tm = tm_del
                best = FlankResult(
                    seq1=seq1_del,
                    seq2=seq2_del,
                    seq1_ori=seq1,
                    seq2_ori=seq2_del[:i] + "-" + seq2_del[i:],
                    seq2_len=len(seq2_del),
                    penalty=10,
                    buldge_pos=i,
                    buldge_seq=2,
                    desc=f"del_pos={i+1:02d} base={seq1[i]} tm_C={tm_del} seq1_del={seq1} seq2={seq2_del[:i]+'-'+seq2_del[i:]}",
                )

    return best


def optimize_right_flank(query_right: str, right_flank: str, kmer_u: str) -> FlankResult:
    """Find best right-flank alignment per original heuristic.

    Baseline uses ``right_flank[:-1]``; then try one-base bulge in DB or query
    with a fixed 10°C penalty.
    """
    if not (query_right and right_flank[:-1]):
        return FlankResult()

    best_tm = calc_tm(kmer_u + query_right, compl(kmer_u + right_flank[:-1]))
    best = FlankResult(
        seq1=query_right,
        seq2=right_flank[:-1],
        seq1_ori=query_right,
        seq2_ori=right_flank[:-1],
        seq2_len=len(right_flank[:-1]),
        penalty=0,
        desc=f"max_pos=0 tm_C={best_tm} seq1={query_right} seq2={right_flank[:-1]}",
    )

    # DB bulge (delete one base from DB flank)
    seq1 = query_right
    seq2 = right_flank
    for i in range(len(seq2)):
        seq2_del = seq2[:i] + seq2[i + 1 :]
        if seq1 and seq2_del:
            tm_del = calc_tm(kmer_u + seq1, compl(kmer_u + seq2_del)) - 10
            if tm_del > best_tm:
                best_tm = tm_del
                best = FlankResult(
                    seq1=seq1,
                    seq2=seq2_del,
                    seq1_ori=seq1[:i] + "-" + seq1[i:],
                    seq2_ori=seq2,
                    seq2_len=len(seq2_del) + 1,
                    penalty=10,
                    buldge_pos=i,
                    buldge_seq=1,
                    desc=f"del_pos={i+1:02d} base={seq2[i]} tm_C={tm_del} seq1={seq1[:i]+'-'+seq1[i:]} seq2={seq2}",
                )

    # Query bulge (delete one base from query flank) using right_flank[:-2]
    seq1 = query_right
    seq2_del = right_flank[:-2]
    for i in range(len(seq1)):
        seq1_del = seq1[:i] + seq1[i + 1 :]
        if seq1_del and seq2_del:
            tm_del = calc_tm(kmer_u + seq1_del, compl(kmer_u + seq2_del)) - 10
            if tm_del > best_tm:
                best_tm = tm_del
                best = FlankResult(
                    seq1=seq1_del,
                    seq2=seq2_del,
                    seq1_ori=seq1,
                    seq2_ori=seq2_del[:i] + "-" + seq2_del[i:],
                    seq2_len=len(seq2_del),
                    penalty=10,
                    buldge_pos=i,
                    buldge_seq=2,
                    desc=f"del_pos={i+1:02d} base={seq1[i]} tm_C={tm_del} seq1_del={seq1} seq2={seq2_del[:i]+'-'+seq2_del[i:]}",
                )

    return best


# ------------------------------
# Main search
# ------------------------------


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fasta", required=True, help="Input FASTA file path")
    p.add_argument(
        "--query",
        required=True,
        help="Query sequence from which to take k-mer windows",
    )
    p.add_argument("--k", type=int, default=4, help="Window size (k-mer length)")
    p.add_argument(
        "--min-tm",
        type=float,
        default=30.0,
        help="Minimum Tm to report a hit (default: 30.0)",
    )
    p.add_argument(
        "--min-identity",
        type=float,
        default=65.0,
        help="Minimum identity percent to report (default: 65)",
    )
    p.add_argument(
        "--cpus",
        type=int,
        default=0,
        help="CPU processes (0=all, 1=disable parallel)",
    )
    return p.parse_args()


def process_match(
    name: str,
    seq: str,
    kmer_u: str,
    start0: int,
    win_idx: int,
    k: int,
    query: str,
) -> Tuple[float, int, int, str, str]:
    """Compute best-flank Tm and alignment strings for one match.

    Returns a tuple of (Tm, db_start1, db_end1, q_align, db_align).
    """
    start1 = start0 + 1
    end1 = start0 + len(kmer_u)

    qlen = len(query)
    left_len = win_idx + 1
    right_len = (qlen - (win_idx + k)) + 1

    left_start = max(0, start0 - left_len)
    left_flank = seq[left_start:start0]
    right_end = min(len(seq), start0 + len(kmer_u) + right_len)
    right_flank = seq[start0 + len(kmer_u) : right_end]

    query_left = query[:win_idx]
    query_right = query[win_idx + len(kmer_u) :]

    left_res = optimize_left_flank(query_left, left_flank, kmer_u)
    right_res = optimize_right_flank(query_right, right_flank, kmer_u)

    seq1 = left_res.seq1 + kmer_u + right_res.seq1
    seq2 = left_res.seq2 + kmer_u + right_res.seq2
    tm = calc_tm(seq1, compl(seq2)) - left_res.penalty - right_res.penalty

    dbstart = start1 - left_res.seq2_len
    dbend = end1 + right_res.seq2_len
    qalign = left_res.seq1_ori + kmer_u + right_res.seq1_ori
    dbalign = left_res.seq2_ori + kmer_u + right_res.seq2_ori

    return tm, dbstart, dbend, qalign, dbalign


def identity_percent(qalign: str, dbalign: str) -> float:
    matches = 0
    comp = 0
    for a, b in zip(qalign, dbalign):
        comp += 1
        if a == b:
            matches += 1
    if comp == 0:
        return 0.0
    return (matches / comp) * 100.0


def search_and_print(records: List[Tuple[str, str]], query: str, k: int, revcomp_flag: bool, min_tm: float) -> None:
    kmers = sliding_kmers(query, k)
    if not kmers:
        raise SystemExit("No k-mers generated; check --query and --k")

    upper_records = [(name, seq.upper()) for name, seq in records]
    output_data = set()

    for win_idx, kmer in enumerate(kmers):
        kmer_u = kmer.upper()
        # '+' strand with Tm evaluation
        for name, seq in upper_records:
            for pos0 in find_all(seq, kmer_u):
                tm, dbstart, dbend, qalign, dbalign = process_match(
                    name, seq, kmer_u, pos0, win_idx, k, query
                )
                ident = identity_percent(qalign, dbalign)
                if tm > min_tm and ident > MIN_IDENTITY:
                    key = (tm, name, dbstart, dbend)
                    if key not in output_data:
                        output_data.add(key)
                        if revcomp_flag:
                            print(f"{tm}\t{name}\t{dbstart}\t{dbend}\t-\t{ident:.1f}\t{qalign}\t{dbalign}")
                        else:
                            print(f"{tm}\t{name}\t{dbstart}\t{dbend}\t+\t{ident:.1f}\t{qalign}\t{dbalign}")

def search_record(name: str, seq: str, query: str, k: int, revcomp_flag: bool, min_tm: float, min_identity: float) -> List[Tuple[float, str, int, int, str, float, str, str]]:
    kmers = sliding_kmers(query, k)
    if not kmers:
        return []
    seq = seq.upper()
    rows: List[Tuple[float, str, int, int, str, float, str, str]] = []
    output_data = set()
    for win_idx, kmer in enumerate(kmers):
        kmer_u = kmer.upper()
        for pos0 in find_all(seq, kmer_u):
            tm, dbstart, dbend, qalign, dbalign = process_match(
                name, seq, kmer_u, pos0, win_idx, k, query
            )
            ident = identity_percent(qalign, dbalign)
            if tm > min_tm and ident > min_identity:
                key = (tm, name, dbstart, dbend)
                if key not in output_data:
                    output_data.add(key)
                    strand = '-' if revcomp_flag else '+'
                    rows.append((tm, name, dbstart, dbend, strand, ident, qalign, dbalign))
    return rows


def process_entry(record: Tuple[str, str], query: str, k: int, min_tm: float, min_identity: float) -> List[str]:
    name, seq = record
    plus_rows = search_record(name, seq, query, k, False, min_tm, min_identity)
    minus_rows = search_record(name, seq, revcomp(query), k, True, min_tm, min_identity)
    rows = plus_rows + minus_rows
    # Sort by Tm descending within this FASTA entry
    rows.sort(key=lambda r: r[0], reverse=True)
    # Format lines with Tm to 2 decimals; identity unchanged (1 decimal)
    lines: List[str] = []
    for tm, nm, dbs, dbe, strand, ident, qalign, dbalign in rows:
        lines.append(f"{tm:.2f}\t{nm}\t{dbs}\t{dbe}\t{strand}\t{ident:.1f}\t{qalign}\t{dbalign}")
    return lines


def main() -> None:
    args = parse_args()
    records = read_fasta(args.fasta)
    if not records:
        raise SystemExit(f"No sequences found in FASTA: {args.fasta}")
    min_identity = float(args.min_identity)
    min_tm = float(args.min_tm)
    cpus = args.cpus if args.cpus != 0 else (os.cpu_count() or 1)

    if cpus <= 1:
        for rec in records:
            lines = process_entry(rec, args.query, args.k, min_tm, min_identity)
            for line in lines:
                print(line, flush=True)
    else:
        with ProcessPoolExecutor(max_workers=cpus) as ex:
            futures = {
                ex.submit(process_entry, rec, args.query, args.k, min_tm, min_identity): rec[0]
                for rec in records
            }
            for fut in as_completed(futures):
                for line in fut.result():
                    print(line, flush=True)

if __name__ == "__main__":
    main()
