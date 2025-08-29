#!/usr/bin/env node
import fs from 'node:fs/promises';
import path from 'node:path';
import url from 'node:url';

// Import the browser-core module (ESM)
import { parseFasta, parseQueries, tmblast } from '../web/tmblast.js';

function parseArgs(argv) {
  const args = {};
  for (let i = 0; i < argv.length; i++) {
    const a = argv[i];
    if (!a.startsWith('--')) continue;
    const key = a.slice(2);
    const next = argv[i + 1];
    if (!next || next.startsWith('--')) {
      args[key] = true;
    } else {
      args[key] = next;
      i++;
    }
  }
  return args;
}

async function main() {
  const argv = process.argv.slice(2);
  const args = parseArgs(argv);

  if (!args.fasta) {
    console.error('Usage: node scripts/tmblast_cli.mjs --fasta <file.fa> (--queries-file <file>|--query <seq> ...) [--k 4] [--min-tm 30] [--min-identity 50] [--max-expansions 100000]');
    process.exit(1);
  }

  const fastaTxt = await fs.readFile(args.fasta, 'utf8');
  const records = parseFasta(fastaTxt);

  let queries = [];
  if (args['queries-file']) {
    const qTxt = await fs.readFile(args['queries-file'], 'utf8');
    queries = parseQueries(qTxt);
  } else if (args.query) {
    const arr = Array.isArray(args.query) ? args.query : [args.query];
    queries = arr.map((q, i) => [ `q${i+1}`, String(q).toUpperCase() ]);
  } else {
    console.error('Provide --queries-file or one/more --query');
    process.exit(1);
  }

  const k = args.k ? Number(args.k) : 4;
  const minTm = args['min-tm'] ? Number(args['min-tm']) : 30.0;
  const minIdentity = args['min-identity'] ? Number(args['min-identity']) : 50.0;
  const maxExp = args['max-expansions'] ? Number(args['max-expansions']) : 100000;

  const opts = {
    primerConc: args['primer-conc'] ? Number(args['primer-conc']) : 100,
    Na: args.Na ? Number(args.Na) : 0,
    K: args.K ? Number(args.K) : 50,
    Tris: args.Tris ? Number(args.Tris) : 10,
    Mg: args.Mg ? Number(args.Mg) : 1.5,
    dNTPs: args.dNTPs ? Number(args.dNTPs) : 0.2,
  };

  const lines = tmblast(records, queries, k, minTm, minIdentity, maxExp, opts);
  process.stdout.write(lines.join('\n') + (lines.length ? '\n' : ''));
}

main().catch((err) => {
  console.error('[tmblast-cli] error:', err);
  process.exit(2);
});

