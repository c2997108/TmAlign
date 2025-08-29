#!/usr/bin/env python3
"""
TmBLAST: Batch Tm-based local alignment search with degenerate queries.

This tool extends TmAlign to:
- Accept multiple query sequences (from CLI or file/FASTA).
- Allow IUPAC degenerate bases in queries and expand them to all A/C/G/T
  sequences for evaluation.

Output columns per hit:
  query_id  expanded_query  Tm  contig  start  end  strand  identity  q_align  db_align

Notes:
- Expansion can grow exponentially; use --max-expansions to guard.
- Results are deduplicated per expanded query and DB coordinates within a query.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed

# Reuse core logic from TmAlign
import TmAlign as tmalign


IUPAC_MAP: Dict[str, str] = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT',
    'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT',
    'H': 'ACT', 'V': 'ACG', 'N': 'ACGT',
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fasta', required=True, help='Target FASTA database')
    # Queries can be provided multiple times or via a file
    p.add_argument('--query', action='append', default=[], help='Query sequence (IUPAC allowed). Can repeat.')
    p.add_argument('--queries-file', help='File of queries (FASTA or lines).')
    p.add_argument('--k', type=int, default=4, help='k-mer window size (default: 4)')
    p.add_argument('--min-tm', type=float, default=30.0, help='Minimum Tm to report (default: 30.0)')
    p.add_argument('--min-identity', type=float, default=50.0, help='Minimum identity percent (default: 50)')
    p.add_argument('--cpus', type=int, default=0, help='CPU processes (0=all, 1=disable parallel)')
    p.add_argument('--max-expansions', type=int, default=100000, help='Guardrail for max expansions per query (default: 100000)')
    p.add_argument('--quiet', action='store_true', help='Suppress warnings to stderr')
    return p.parse_args()


def is_fasta_file(path: str) -> bool:
    _, ext = os.path.splitext(path)
    if ext.lower() in {'.fa', '.fasta', '.fna', '.ffn'}:
        return True
    # Fallback: sniff first non-empty char
    try:
        with open(path, 'r') as fh:
            for line in fh:
                if line.strip():
                    return line.startswith('>')
    except OSError:
        pass
    return False


def load_queries(args: argparse.Namespace) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    # From repeated --query
    for i, q in enumerate(args.query, start=1):
        q = q.strip().upper()
        if not q:
            continue
        records.append((f'q{i}', q))

    # From file
    if args.queries_file:
        if is_fasta_file(args.queries_file):
            records.extend(tmalign.read_fasta(args.queries_file))
        else:
            with open(args.queries_file, 'r') as fh:
                idx = 1
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    # Support optional "id<tab>seq" or just "seq"
                    if '\t' in line:
                        qid, seq = line.split('\t', 1)
                    elif ' ' in line:
                        parts = line.split(None, 1)
                        if len(parts) == 2 and all(c.isprintable() for c in parts[1]):
                            qid, seq = parts
                        else:
                            qid, seq = f'q{idx}', line
                    else:
                        qid, seq = f'q{idx}', line
                    records.append((qid, seq.upper()))
                    idx += 1

    if not records:
        raise SystemExit('No queries provided. Use --query and/or --queries-file.')
    return records


def expand_degenerate(seq: str, max_expansions: int, quiet: bool = False) -> List[str]:
    # Pre-map to lists of options per position; validate characters
    options: List[str] = []
    for ch in seq.upper():
        if ch not in IUPAC_MAP:
            raise SystemExit(f'Unsupported base in query: {ch!r}')
        options.append(IUPAC_MAP[ch])

    # Early exit: no degeneracy
    if all(len(opt) == 1 for opt in options):
        return [seq.upper()]

    # Iteratively build expansions; truncate if exceeding max_expansions
    expansions: List[str] = ['']
    for opt in options:
        new_exp: List[str] = []
        for prefix in expansions:
            for base in opt:
                new_exp.append(prefix + base)
                if len(new_exp) == max_expansions:
                    if not quiet:
                        print(f'[warn] expansions truncated at {max_expansions}', file=sys.stderr)
                    return new_exp
        expansions = new_exp
    return expansions


def search_one_query(records: List[Tuple[str, str]], query: str, k: int, min_tm: float, min_identity: float, cpus: int) -> List[Tuple[float, str, int, int, str, float, str, str]]:
    # Parallelize across FASTA records using tmalign.search_record
    rows: List[Tuple[float, str, int, int, str, float, str, str]] = []
    if cpus <= 1:
        for name, seq in records:
            rows.extend(tmalign.search_record(name, seq, query, k, False, min_tm, min_identity))
            rows.extend(tmalign.search_record(name, seq, tmalign.revcomp(query), k, True, min_tm, min_identity))
    else:
        with ProcessPoolExecutor(max_workers=cpus) as ex:
            futs = []
            for name, seq in records:
                futs.append(ex.submit(tmalign.search_record, name, seq, query, k, False, min_tm, min_identity))
                futs.append(ex.submit(tmalign.search_record, name, seq, tmalign.revcomp(query), k, True, min_tm, min_identity))
            for fut in as_completed(futs):
                rows.extend(fut.result())
    # Sort by Tm descending then contig/start
    rows.sort(key=lambda r: (r[0], r[1], r[2], r[3]), reverse=True)
    return rows


def main() -> None:
    args = parse_args()
    records = tmalign.read_fasta(args.fasta)
    if not records:
        raise SystemExit(f'No sequences found in FASTA: {args.fasta}')

    cpus = args.cpus if args.cpus != 0 else (os.cpu_count() or 1)
    min_tm = float(args.min_tm)
    min_identity = float(args.min_identity)

    queries = load_queries(args)

    for qid, qseq in queries:
        exps = expand_degenerate(qseq, args.max_expansions, args.quiet)
        # Deduplicate expansions if degenerate map overlaps
        seen = set()
        unique_exps = []
        for e in exps:
            if e not in seen:
                seen.add(e)
                unique_exps.append(e)

        for idx, exp in enumerate(unique_exps, start=1):
            rows = search_one_query(records, exp, args.k, min_tm, min_identity, cpus)
            # dedupe per expanded query to avoid duplicate DB coordinates
            emitted = set()
            for tm, name, dbs, dbe, strand, ident, qalign, dbalign in rows:
                key = (name, dbs, dbe, strand, qalign, dbalign)
                if key in emitted:
                    continue
                emitted.add(key)
                print(f"{qid}\t{exp}\t{tm:.2f}\t{name}\t{dbs}\t{dbe}\t{strand}\t{ident:.1f}\t{qalign}\t{dbalign}")


if __name__ == '__main__':
    main()
