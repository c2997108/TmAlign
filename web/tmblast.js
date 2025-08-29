import { calcTmSanta, complement } from './nn_tm.js';

export const IUPAC_MAP = {
  A: 'A', C: 'C', G: 'G', T: 'T',
  R: 'AG', Y: 'CT', S: 'GC', W: 'AT',
  K: 'GT', M: 'AC', B: 'CGT', D: 'AGT',
  H: 'ACT', V: 'ACG', N: 'ACGT',
};

export function parseFasta(text) {
  const lines = text.split(/\r?\n/);
  const out = [];
  let header = null;
  let seq = [];
  for (const raw of lines) {
    const line = raw.trim();
    if (!line) continue;
    if (line.startsWith('>')) {
      if (header) out.push([header, seq.join('').toUpperCase()]);
      header = line.slice(1).split(/[\s\t]/)[0];
      seq = [];
    } else {
      seq.push(line);
    }
  }
  if (header) out.push([header, seq.join('').toUpperCase()]);
  return out;
}

export function parseFastaOrRaw(text, defaultName='entry1') {
  const recs = parseFasta(text);
  if (recs.length) return recs;
  const seq = text.replace(/\s+/g, '').toUpperCase();
  if (/^[ACGTN-]+$/.test(seq) && seq.length > 0) {
    return [[defaultName, seq]];
  }
  return [];
}

export function parseQueries(text) {
  // Accept FASTA; else lines of seq or "id\tseq"
  if (/^>/.test(text.trim())) return parseFasta(text);
  const lines = text.split(/\r?\n/);
  const out = [];
  let i = 1;
  for (const raw of lines) {
    const line = raw.trim();
    if (!line) continue;
    if (line.includes('\t')) {
      const [id, seq] = line.split('\t');
      out.push([id || `q${i++}`, (seq || '').toUpperCase()]);
    } else {
      out.push([`q${i++}`, line.toUpperCase()]);
    }
  }
  return out;
}

export function revcomp(seq) {
  const map = { A:'T', C:'G', G:'C', T:'A',
                R:'Y', Y:'R', S:'S', W:'W', K:'M', M:'K',
                B:'V', D:'H', H:'D', V:'B', N:'N' };
  return seq.toUpperCase().split('').map(c => map[c] || 'N').reverse().join('');
}

export function compl(seq) {
  return revcomp(seq).split('').reverse().join('');
}

export function slidingKmers(seq, k) {
  const s = seq.toUpperCase();
  if (k <= 0 || k > s.length) return [];
  const arr = [];
  for (let i = 0; i <= s.length - k; i++) arr.push(s.slice(i, i + k));
  return arr;
}

export function findAll(hay, needle) {
  const out = [];
  let idx = 0;
  for (;;) {
    const found = hay.indexOf(needle, idx);
    if (found === -1) break;
    out.push(found);
    idx = found + 1;
  }
  return out;
}

export function expandDegenerate(seq, maxExpansions) {
  const opts = seq.toUpperCase().split('').map(c => {
    if (!IUPAC_MAP[c]) throw new Error(`Unsupported base: ${c}`);
    return IUPAC_MAP[c];
  });
  if (opts.every(s => s.length === 1)) return [seq.toUpperCase()];
  let exps = [''];
  for (const s of opts) {
    const next = [];
    for (const p of exps) {
      for (const b of s) {
        next.push(p + b);
        if (next.length === maxExpansions) return next;
      }
    }
    exps = next;
  }
  return exps;
}

export function calcTm(seq, dbseq, opts) {
  // Use exact SantaLucia port to align with Python
  return calcTmSanta(seq, complement(dbseq), { ...opts, ionCorr: true });
}

function optimizeLeftFlank(queryLeft, leftFlank, kmer, opts) {
  if (!queryLeft || leftFlank.length < 1) return { seq1:'', seq2:'', seq1_ori:'', seq2_ori:'', seq2_len:0, penalty:0 };
  const baseTm = calcTm(queryLeft + kmer, leftFlank.slice(1) + kmer, opts);
  let best = { tm: baseTm, seq1: queryLeft, seq2: leftFlank.slice(1), seq1_ori: queryLeft, seq2_ori: leftFlank.slice(1), seq2_len: leftFlank.length - 1, penalty: 0 };

  // DB bulge: delete one base from DB flank
  for (let i = 0; i < leftFlank.length; i++) {
    const seq2_del = leftFlank.slice(0, i) + leftFlank.slice(i+1);
    if (queryLeft && seq2_del) {
      const tm = calcTm(queryLeft + kmer, seq2_del + kmer, opts) - 10;
      if (tm > best.tm) {
        best = { tm, seq1: queryLeft, seq2: seq2_del, seq1_ori: queryLeft.slice(0,i)+'-'+queryLeft.slice(i), seq2_ori: leftFlank, seq2_len: seq2_del.length + 1, penalty:10 };
      }
    }
  }

  // Query bulge: delete one base from query flank; use leftFlank[2:]
  const seq2_del0 = leftFlank.slice(2);
  for (let i = 0; i < queryLeft.length; i++) {
    const seq1_del = queryLeft.slice(0,i) + queryLeft.slice(i+1);
    if (seq1_del && seq2_del0) {
      const tm = calcTm(seq1_del + kmer, seq2_del0 + kmer, opts) - 10;
      if (tm > best.tm) {
        best = { tm, seq1: seq1_del, seq2: seq2_del0, seq1_ori: queryLeft, seq2_ori: seq2_del0.slice(0,i)+'-'+seq2_del0.slice(i), seq2_len: seq2_del0.length, penalty:10 };
      }
    }
  }
  return best;
}

function optimizeRightFlank(queryRight, rightFlank, kmer, opts) {
  if (!queryRight || rightFlank.length < 1) return { seq1:'', seq2:'', seq1_ori:'', seq2_ori:'', seq2_len:0, penalty:0 };
  const baseTm = calcTm(kmer + queryRight, kmer + rightFlank.slice(0, -1), opts);
  let best = { tm: baseTm, seq1: queryRight, seq2: rightFlank.slice(0,-1), seq1_ori: queryRight, seq2_ori: rightFlank.slice(0,-1), seq2_len: rightFlank.length - 1, penalty: 0 };

  // DB bulge
  for (let i = 0; i < rightFlank.length; i++) {
    const seq2_del = rightFlank.slice(0,i) + rightFlank.slice(i+1);
    if (queryRight && seq2_del) {
      const tm = calcTm(kmer + queryRight, kmer + seq2_del, opts) - 10;
      if (tm > best.tm) {
        best = { tm, seq1: queryRight, seq2: seq2_del, seq1_ori: queryRight.slice(0,i)+'-'+queryRight.slice(i), seq2_ori: rightFlank, seq2_len: seq2_del.length + 1, penalty:10 };
      }
    }
  }

  // Query bulge: delete one base from query flank; use rightFlank[:-2]
  const seq2_del0 = rightFlank.slice(0, -2);
  for (let i = 0; i < queryRight.length; i++) {
    const seq1_del = queryRight.slice(0,i) + queryRight.slice(i+1);
    if (seq1_del && seq2_del0) {
      const tm = calcTm(kmer + seq1_del, kmer + seq2_del0, opts) - 10;
      if (tm > best.tm) {
        best = { tm, seq1: seq1_del, seq2: seq2_del0, seq1_ori: queryRight, seq2_ori: seq2_del0.slice(0,i)+'-'+seq2_del0.slice(i), seq2_len: seq2_del0.length, penalty:10 };
      }
    }
  }
  return best;
}

function processMatch(name, seq, kmer, start0, winIdx, k, query, opts) {
  const start1 = start0 + 1;
  const end1 = start0 + kmer.length;
  const qlen = query.length;
  const left_len = winIdx + 1;
  const right_len = (qlen - (winIdx + k)) + 1;

  const left_start = Math.max(0, start0 - left_len);
  const left_flank = seq.slice(left_start, start0);
  const right_end = Math.min(seq.length, start0 + kmer.length + right_len);
  const right_flank = seq.slice(start0 + kmer.length, right_end);

  const query_left = query.slice(0, winIdx);
  const query_right = query.slice(winIdx + kmer.length);

  const left_res = optimizeLeftFlank(query_left, left_flank, kmer, opts);
  const right_res = optimizeRightFlank(query_right, right_flank, kmer, opts);

  const seq1 = left_res.seq1 + kmer + right_res.seq1;
  const seq2 = left_res.seq2 + kmer + right_res.seq2;
  const tm = calcTm(seq1, seq2, opts) - left_res.penalty - right_res.penalty;

  const dbstart = start1 - left_res.seq2_len;
  const dbend = end1 + right_res.seq2_len;
  const qalign = (left_res.seq1_ori || '') + kmer + (right_res.seq1_ori || '');
  const dbalign = (left_res.seq2_ori || '') + kmer + (right_res.seq2_ori || '');

  return { tm, dbstart, dbend, qalign, dbalign };
}

function identityPercent(qalign, dbalign) {
  let match = 0, comp = 0;
  const n = Math.min(qalign.length, dbalign.length);
  for (let i = 0; i < n; i++) {
    comp++;
    if (qalign[i] === dbalign[i]) match++;
  }
  return comp ? (match / comp) * 100.0 : 0.0;
}

export function searchRecord(name, seq, query, k, revFlag, minTm, minIdentity, opts) {
  const kmers = slidingKmers(query, k);
  if (!kmers.length) return [];
  const seqU = seq.toUpperCase();
  const rows = [];
  const seen = new Set();
  kmers.forEach((kmer, winIdx) => {
    const kmerU = kmer.toUpperCase();
    const positions = findAll(seqU, kmerU);
    for (const pos0 of positions) {
      const { tm, dbstart, dbend, qalign, dbalign } = processMatch(name, seqU, kmerU, pos0, winIdx, k, query, opts);
      const ident = identityPercent(qalign, dbalign);
      if (tm > minTm && ident > minIdentity) {
        const key = `${tm.toFixed(2)}\t${name}\t${dbstart}\t${dbend}`;
        if (!seen.has(key)) {
          seen.add(key);
          rows.push([tm, name, dbstart, dbend, revFlag ? '-' : '+', ident, qalign, dbalign]);
        }
      }
    }
  });
  return rows;
}

export function tmblast(records, queries, k, minTm, minIdentity, maxExpansions, opts) {
  // records: [ [name, seq], ... ]
  // queries: [ [id, seq], ... ] (IUPAC allowed)
  const lines = [];
  for (const [qid, qseqRaw] of queries) {
    const exps = expandDegenerate(qseqRaw, maxExpansions);
    const uniq = Array.from(new Set(exps));
    for (const exp of uniq) {
      const rowsAll = [];
      for (const [name, seq] of records) {
        rowsAll.push(...searchRecord(name, seq, exp, k, false, minTm, minIdentity, opts));
        rowsAll.push(...searchRecord(name, seq, revcomp(exp), k, true, minTm, minIdentity, opts));
      }
      rowsAll.sort((a, b) => b[0] - a[0]);
      const emitted = new Set();
      for (const [tm, name, dbs, dbe, strand, ident, qalign, dbalign] of rowsAll) {
        const key = `${name}\t${dbs}\t${dbe}\t${strand}\t${qalign}\t${dbalign}`;
        if (emitted.has(key)) continue;
        emitted.add(key);
        lines.push(`${qid}\t${exp}\t${tm.toFixed(2)}\t${name}\t${dbs}\t${dbe}\t${strand}\t${ident.toFixed(1)}\t${qalign}\t${dbalign}`);
      }
    }
  }
  return lines;
}

export function tmblastRows(records, queries, k, minTm, minIdentity, maxExpansions, opts) {
  const rowsOut = [];
  for (const [qid, qseqRaw] of queries) {
    const exps = expandDegenerate(qseqRaw, maxExpansions);
    const uniq = Array.from(new Set(exps));
    for (const exp of uniq) {
      const rowsAll = [];
      for (const [name, seq] of records) {
        rowsAll.push(...searchRecord(name, seq, exp, k, false, minTm, minIdentity, opts));
        rowsAll.push(...searchRecord(name, seq, revcomp(exp), k, true, minTm, minIdentity, opts));
      }
      rowsAll.sort((a, b) => b[0] - a[0]);
      const emitted = new Set();
      for (const [tm, name, dbs, dbe, strand, ident, qalign, dbalign] of rowsAll) {
        const key = `${name}\t${dbs}\t${dbe}\t${strand}\t${qalign}\t${dbalign}`;
        if (emitted.has(key)) continue;
        emitted.add(key);
        rowsOut.push({ qid, exp, tm, name, start: dbs, end: dbe, strand, ident, qalign, dbalign });
      }
    }
  }
  return rowsOut;
}
