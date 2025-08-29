// TmBLAST browser bundle (no modules) – defines window.TmBLAST
(function(){
  const T = {};

  // ===== nn_tm (exact SantaLucia port) =====
  const DNA_NN_table = { 'init':[0.2,-5.7], 'init_A/T':[2.2,6.9], 'init_G/C':[0,0], 'init_oneG/C':[0,0], 'init_allA/T':[0,0], 'init_5T/A':[0,0], 'sym':[0,-1.4], 'AA/TT':[-7.6,-21.3], 'AT/TA':[-7.2,-20.4], 'TA/AT':[-7.2,-20.4], 'CA/GT':[-8.5,-22.7], 'GT/CA':[-8.4,-22.4], 'CT/GA':[-7.8,-21.0], 'GA/CT':[-8.2,-22.2], 'CG/GC':[-10.6,-27.2], 'GC/CG':[-9.8,-24.4], 'GG/CC':[-8.0,-19.0] };
  const DNA_IMM_table = { 'AG/TT':[1.0,0.9],'AT/TG':[-2.5,-8.3],'CG/GT':[-4.1,-11.7],'CT/GG':[-2.8,-8.0],'GG/CT':[3.3,10.4],'GG/TT':[5.8,16.3],'GT/CG':[-4.4,-12.3],'GT/TG':[4.1,9.5],'TG/AT':[-0.1,-1.7],'TG/GT':[-1.4,-6.2],'TT/AG':[-1.3,-5.3],'AA/TG':[-0.6,-2.3],'AG/TA':[-0.7,-2.3],'CA/GG':[-0.7,-2.3],'CG/GA':[-4.0,-13.2],'GA/CG':[-0.6,-1.0],'GG/CA':[0.5,3.2],'TA/AG':[0.7,0.7],'TG/AA':[3.0,7.4],'AC/TT':[0.7,0.2],'AT/TC':[-1.2,-6.2],'CC/GT':[-0.8,-4.5],'CT/GC':[-1.5,-6.1],'GC/CT':[2.3,5.4],'GT/CC':[5.2,13.5],'TC/AT':[1.2,0.7],'TT/AC':[1.0,0.7],'AA/TC':[2.3,4.6],'AC/TA':[5.3,14.6],'CA/GC':[1.9,3.7],'CC/GA':[0.6,-0.6],'GA/CC':[5.2,14.2],'GC/CA':[-0.7,-3.8],'TA/AC':[3.4,8.0],'TC/AA':[7.6,20.2],'AA/TA':[1.2,1.7],'CA/GA':[-0.9,-4.2],'GA/CA':[-2.9,-9.8],'TA/AA':[4.7,12.9],'AC/TC':[0.0,-4.4],'CC/GC':[-1.5,-7.2],'GC/CC':[3.6,8.9],'TC/AC':[6.1,16.4],'AG/TG':[-3.1,-9.5],'CG/GG':[-4.9,-15.3],'GG/CG':[-6.0,-15.8],'TG/AG':[1.6,3.6],'AT/TT':[-2.7,-10.8],'CT/GT':[-5.0,-15.8],'GT/CT':[-2.2,-8.4],'TT/AT':[0.2,-1.5] };
  const DNA_TMM_table = { 'AA/TA':[-3.1,-7.8],'TA/AA':[-2.5,-6.3],'CA/GA':[-4.3,-10.7],'GA/CA':[-8.0,-22.5],'AC/TC':[-0.1,0.5],'TC/AC':[-0.7,-1.3],'CC/GC':[-2.1,-5.1],'GC/CC':[-3.9,-10.6],'AG/TG':[-1.1,-2.1],'TG/AG':[-1.1,-2.7],'CG/GG':[-3.8,-9.5],'GG/CG':[-0.7,-19.2],'AT/TT':[-2.4,-6.5],'TT/AT':[-3.2,-8.9],'CT/GT':[-6.1,-16.9],'GT/CT':[-7.4,-21.2],'AA/TC':[-1.6,-4.0],'AC/TA':[-1.8,-3.8],'CA/GC':[-2.6,-5.9],'CC/GA':[-2.7,-6.0],'GA/CC':[-5.0,-13.8],'GC/CA':[-3.2,-7.1],'TA/AC':[-2.3,-5.9],'TC/AA':[-2.7,-7.0],'AC/TT':[-0.9,-1.7],'AT/TC':[-2.3,-6.3],'CC/GT':[-3.2,-8.0],'CT/GC':[-3.9,-10.6],'GC/CT':[-4.9,-13.5],'GT/CC':[-3.0,-7.8],'TC/AT':[-2.5,-6.3],'TT/AC':[-0.7,-1.2],'AA/TG':[-1.9,-4.4],'AG/TA':[-2.5,-5.9],'CA/GG':[-3.9,-9.6],'CG/GA':[-6.0,-15.5],'GA/CG':[-4.3,-11.1],'GG/CA':[-4.6,-11.4],'TA/AG':[-2.0,-4.7],'TG/AA':[-2.4,-5.8],'AG/TT':[-3.2,-8.7],'AT/TG':[-3.5,-9.4],'CG/GT':[-3.8,-9.0],'CT/GG':[-6.6,-18.7],'GG/CT':[-5.7,-15.9],'GT/CG':[-5.9,-16.1],'TG/AT':[-3.9,-10.5],'TT/AG':[-3.6,-9.8] };
  const R = 1.9872;
  function complBase(b){ const m={A:'T',C:'G',G:'C',T:'A','-':'-','N':'N'}; return m[b]||'N'; }
  function complement(seq){ return seq.toUpperCase().split('').map(complBase).join(''); }
  function ion_correction(Na=0,K=0,Tris=0,Mg=0,dNTPs=0,seq_len=0){ const Monovalent_mmol = Number(Na||0)+Number(K||0)+(Number(Tris||0)/2); const mg=Number(Mg||0), dn=Number(dNTPs||0); const Na_eq_mmol = dn>=mg ? Monovalent_mmol : (Monovalent_mmol + (120*Math.sqrt(Math.max(0, mg-dn)))); const Na_eq_mol = Na_eq_mmol/1000.0; return 0.368 * (seq_len - 1) * Math.log(Na_eq_mol); }
  function reverseKey(NN){ return NN.split('').reverse().join(''); }
  function NN_Tm_JS({seq, compl_seq, primer_conc=0, Na=0, K=0, Tris=0, Mg=0, dNTPs=0, ion_corr=false}){
    if(!seq) throw new Error('Please provide an input sequence!');
    const clean = s=>String(s).replace(/\s+/g,'').toUpperCase().replace(/[^ACGTN\-]/g,'');
    seq = clean(seq);
    if(!compl_seq) compl_seq = complement(seq);
    compl_seq = clean(compl_seq);
    let temp_seq = String(seq), temp_compl_seq=String(compl_seq);
    let dH=0.0, dS=0.0; const dHi=0,dSi=1;
    const left_tmm = temp_compl_seq.slice(0,2).split('').reverse().join('') + '/' + temp_seq.slice(0,2).split('').reverse().join('');
    if (DNA_TMM_table[left_tmm]){ dH+=DNA_TMM_table[left_tmm][dHi]; dS+=DNA_TMM_table[left_tmm][dSi]; temp_compl_seq=temp_compl_seq.slice(1); temp_seq=temp_seq.slice(1); }
    const right_tmm = temp_seq.slice(-2) + '/' + temp_compl_seq.slice(-2);
    if (DNA_TMM_table[right_tmm]){ dH+=DNA_TMM_table[right_tmm][dHi]; dS+=DNA_TMM_table[right_tmm][dSi]; temp_compl_seq=temp_compl_seq.slice(0,-1); temp_seq=temp_seq.slice(0,-1); }
    dH += DNA_NN_table['init'][dHi]; dS += DNA_NN_table['init'][dSi];
    const terminals = seq[0] + '/' + compl_seq[0] + seq[seq.length-1] + '/' + compl_seq[compl_seq.length-1];
    const count_AT = ((terminals.match(/A\/T/g)||[]).length) + ((terminals.match(/T\/A/g)||[]).length);
    dH += DNA_NN_table['init_A/T'][dHi] * count_AT; dS += DNA_NN_table['init_A/T'][dSi] * count_AT;
    for(let i=0;i<temp_seq.length-1;i++){
      const NN = temp_seq.slice(i,i+2) + '/' + temp_compl_seq.slice(i,i+2);
      if(DNA_IMM_table[NN]){ dH+=DNA_IMM_table[NN][dHi]; dS+=DNA_IMM_table[NN][dSi]; }
      else if(DNA_IMM_table[reverseKey(NN)]){ const k=reverseKey(NN); dH+=DNA_IMM_table[k][dHi]; dS+=DNA_IMM_table[k][dSi]; }
      else if(DNA_NN_table[NN]){ dH+=DNA_NN_table[NN][dHi]; dS+=DNA_NN_table[NN][dSi]; }
      else if(DNA_NN_table[reverseKey(NN)]){ const k=reverseKey(NN); dH+=DNA_NN_table[k][dHi]; dS+=DNA_NN_table[k][dSi]; }
    }
    if(ion_corr){ dS += ion_correction(Na,K,Tris,Mg,dNTPs, seq.length); }
    const x=4; const Ct = (Number(primer_conc)||0)/1e9; const Tm = (1000*dH) / (dS + (R*Math.log(Ct/x))) - 273.15; return Tm;
  }
  T.complement = complement;

  // ===== tmblast core =====
  const IUPAC_MAP = {A:'A',C:'C',G:'G',T:'T',R:'AG',Y:'CT',S:'GC',W:'AT',K:'GT',M:'AC',B:'CGT',D:'AGT',H:'ACT',V:'ACG',N:'ACGT'};
  T.IUPAC_MAP = IUPAC_MAP;

  function parseFasta(text){
    const lines=text.split(/\r?\n/); const out=[]; let header=null, seq=[];
    for (const raw of lines){ const line=raw.trim(); if(!line) continue; if(line.startsWith('>')){ if(header) out.push([header, seq.join('').toUpperCase()]); header=line.slice(1).split(/[\s\t]/)[0]; seq=[]; } else seq.push(line); }
    if(header) out.push([header, seq.join('').toUpperCase()]);
    return out;
  }
  function parseFastaOrRaw(text, defaultName='entry1'){
    const recs = parseFasta(text);
    if (recs.length) return recs;
    const seq = text.replace(/\s+/g,'').toUpperCase();
    if (/^[ACGTN-]+$/.test(seq) && seq.length>0) return [[defaultName, seq]];
    return [];
  }
  function parseQueries(text){ if(/^>/.test(text.trim())) return parseFasta(text); const lines=text.split(/\r?\n/); const out=[]; let i=1; for(const raw of lines){ const line=raw.trim(); if(!line) continue; if(line.includes('\t')){ const [id,seq]=(line.split('\t')); out.push([id||`q${i++}`, (seq||'').toUpperCase()]); } else { out.push([`q${i++}`, line.toUpperCase()]); } } return out; }
  function revcomp(seq){ const m={A:'T',C:'G',G:'C',T:'A',R:'Y',Y:'R',S:'S',W:'W',K:'M',M:'K',B:'V',D:'H',H:'D',V:'B',N:'N'}; return seq.toUpperCase().split('').map(c=>m[c]||'N').reverse().join(''); }
  function compl(seq){ return revcomp(seq).split('').reverse().join(''); }
  function slidingKmers(seq,k){ const s=seq.toUpperCase(); if(k<=0||k>s.length) return []; const arr=[]; for(let i=0;i<=s.length-k;i++) arr.push(s.slice(i,i+k)); return arr; }
  function findAll(hay, needle){ const out=[]; let idx=0; for(;;){ const f=hay.indexOf(needle, idx); if(f===-1) break; out.push(f); idx=f+1; } return out; }
  function expandDegenerate(seq, maxExpansions){ const opts=seq.toUpperCase().split('').map(c=>{ if(!IUPAC_MAP[c]) throw new Error(`Unsupported base: ${c}`); return IUPAC_MAP[c]; }); if(opts.every(s=>s.length===1)) return [seq.toUpperCase()]; let exps=['']; for(const s of opts){ const next=[]; for(const p of exps){ for(const b of s){ next.push(p+b); if(next.length===maxExpansions) return next; } } exps=next; } return exps; }
  function calcTm(seq, dbseq, opts){ return NN_Tm_JS({ seq, compl_seq: complement(dbseq), primer_conc: (opts&&opts.primerConc)||100, Na:(opts&&opts.Na)||0, K:(opts&&opts.K)||50, Tris:(opts&&opts.Tris)||10, Mg:(opts&&opts.Mg)||1.5, dNTPs:(opts&&opts.dNTPs)||0.2, ion_corr:true }); }
  function optimizeLeftFlank(queryLeft, leftFlank, kmer, opts){ if(!queryLeft||leftFlank.length<1) return {seq1:'',seq2:'',seq1_ori:'',seq2_ori:'',seq2_len:0,penalty:0}; const baseTm=calcTm(queryLeft+kmer, leftFlank.slice(1)+kmer, opts); let best={ tm:baseTm, seq1:queryLeft, seq2:leftFlank.slice(1), seq1_ori:queryLeft, seq2_ori:leftFlank.slice(1), seq2_len:leftFlank.length-1, penalty:0}; for(let i=0;i<leftFlank.length;i++){ const seq2_del=leftFlank.slice(0,i)+leftFlank.slice(i+1); if(queryLeft&&seq2_del){ const tm=calcTm(queryLeft+kmer, seq2_del+kmer, opts)-10; if(tm>best.tm){ best={ tm, seq1:queryLeft, seq2:seq2_del, seq1_ori:queryLeft.slice(0,i)+'-'+queryLeft.slice(i), seq2_ori:leftFlank, seq2_len:seq2_del.length+1, penalty:10}; } } } const seq2_del0=leftFlank.slice(2); for(let i=0;i<queryLeft.length;i++){ const seq1_del=queryLeft.slice(0,i)+queryLeft.slice(i+1); if(seq1_del&&seq2_del0){ const tm=calcTm(seq1_del+kmer, seq2_del0+kmer, opts)-10; if(tm>best.tm){ best={ tm, seq1:seq1_del, seq2:seq2_del0, seq1_ori:queryLeft, seq2_ori:seq2_del0.slice(0,i)+'-'+seq2_del0.slice(i), seq2_len:seq2_del0.length, penalty:10}; } } } return best; }
  function optimizeRightFlank(queryRight, rightFlank, kmer, opts){ if(!queryRight||rightFlank.length<1) return {seq1:'',seq2:'',seq1_ori:'',seq2_ori:'',seq2_len:0,penalty:0}; const baseTm=calcTm(kmer+queryRight, kmer+rightFlank.slice(0,-1), opts); let best={ tm:baseTm, seq1:queryRight, seq2:rightFlank.slice(0,-1), seq1_ori:queryRight, seq2_ori:rightFlank.slice(0,-1), seq2_len:rightFlank.length-1, penalty:0}; for(let i=0;i<rightFlank.length;i++){ const seq2_del=rightFlank.slice(0,i)+rightFlank.slice(i+1); if(queryRight&&seq2_del){ const tm=calcTm(kmer+queryRight, kmer+seq2_del, opts)-10; if(tm>best.tm){ best={ tm, seq1:queryRight, seq2:seq2_del, seq1_ori:queryRight.slice(0,i)+'-'+queryRight.slice(i), seq2_ori:rightFlank, seq2_len:seq2_del.length+1, penalty:10}; } } } const seq2_del0=rightFlank.slice(0,-2); for(let i=0;i<queryRight.length;i++){ const seq1_del=queryRight.slice(0,i)+queryRight.slice(i+1); if(seq1_del&&seq2_del0){ const tm=calcTm(kmer+seq1_del, kmer+seq2_del0, opts)-10; if(tm>best.tm){ best={ tm, seq1:seq1_del, seq2:seq2_del0, seq1_ori:queryRight, seq2_ori:seq2_del0.slice(0,i)+'-'+seq2_del0.slice(i), seq2_len:seq2_del0.length, penalty:10}; } } } return best; }
  function processMatch(name, seq, kmer, start0, winIdx, k, query, opts){ const start1=start0+1; const end1=start0+kmer.length; const qlen=query.length; const left_len=winIdx+1; const right_len=(qlen-(winIdx+k))+1; const left_start=Math.max(0, start0-left_len); const left_flank=seq.slice(left_start,start0); const right_end=Math.min(seq.length, start0+kmer.length+right_len); const right_flank=seq.slice(start0+kmer.length, right_end); const query_left=query.slice(0,winIdx); const query_right=query.slice(winIdx+kmer.length); const left_res=optimizeLeftFlank(query_left,left_flank,kmer,opts); const right_res=optimizeRightFlank(query_right,right_flank,kmer,opts); const seq1=(left_res.seq1||'')+kmer+(right_res.seq1||''); const seq2=(left_res.seq2||'')+kmer+(right_res.seq2||''); const tm=calcTm(seq1, seq2, opts) - (left_res.penalty||0) - (right_res.penalty||0); const dbstart=start1 - (left_res.seq2_len||0); const dbend=end1 + (right_res.seq2_len||0); const qalign=(left_res.seq1_ori||'')+kmer+(right_res.seq1_ori||''); const dbalign=(left_res.seq2_ori||'')+kmer+(right_res.seq2_ori||''); return { tm, dbstart, dbend, qalign, dbalign }; }
  function identityPercent(qalign, dbalign){ let m=0,c=0; const n=Math.min(qalign.length, dbalign.length); for(let i=0;i<n;i++){ c++; if(qalign[i]===dbalign[i]) m++; } return c? (m/c)*100.0 : 0.0; }
  function searchRecord(name, seq, query, k, revFlag, minTm, minIdentity, opts){ const kmers=slidingKmers(query,k); if(!kmers.length) return []; const seqU=seq.toUpperCase(); const rows=[]; const seen=new Set(); kmers.forEach((kmer,winIdx)=>{ const kmerU=kmer.toUpperCase(); const positions=findAll(seqU,kmerU); for(const pos0 of positions){ const {tm,dbstart,dbend, qalign, dbalign}=processMatch(name,seqU,kmerU,pos0,winIdx,k,query,opts); const ident=identityPercent(qalign, dbalign); if(tm>minTm && ident>minIdentity){ const key=`${tm.toFixed(2)}\t${name}\t${dbstart}\t${dbend}`; if(!seen.has(key)){ seen.add(key); rows.push([tm,name,dbstart,dbend, revFlag?'-':'+', ident, qalign, dbalign]); } } } }); return rows; }
  function tmblast(records, queries, k, minTm, minIdentity, maxExpansions, opts){ const lines=[]; for(const [qid,qseqRaw] of queries){ const exps=expandDegenerate(qseqRaw, maxExpansions); const uniq=[...new Set(exps)]; for(const exp of uniq){ const rowsAll=[]; for(const [name,seq] of records){ rowsAll.push(...searchRecord(name, seq, exp, k, false, minTm, minIdentity, opts)); rowsAll.push(...searchRecord(name, seq, revcomp(exp), k, true, minTm, minIdentity, opts)); }
        rowsAll.sort((a,b)=>b[0]-a[0]); const emitted=new Set(); for(const [tm,name,dbs,dbe,strand,ident,qalign,dbalign] of rowsAll){ const key=`${name}\t${dbs}\t${dbe}\t${strand}\t${qalign}\t${dbalign}`; if(emitted.has(key)) continue; emitted.add(key); lines.push(`${qid}\t${exp}\t${tm.toFixed(2)}\t${name}\t${dbs}\t${dbe}\t${strand}\t${ident.toFixed(1)}\t${qalign}\t${dbalign}`); }
      }
    }
    return lines;
  }

  function tmblastRows(records, queries, k, minTm, minIdentity, maxExpansions, opts){ const rowsOut=[]; for(const [qid,qseqRaw] of queries){ const exps=expandDegenerate(qseqRaw, maxExpansions); const uniq=[...new Set(exps)]; for(const exp of uniq){ const rowsAll=[]; for(const [name,seq] of records){ rowsAll.push(...searchRecord(name, seq, exp, k, false, minTm, minIdentity, opts)); rowsAll.push(...searchRecord(name, seq, revcomp(exp), k, true, minTm, minIdentity, opts)); } rowsAll.sort((a,b)=>b[0]-a[0]); const emitted=new Set(); for(const [tm,name,dbs,dbe,strand,ident,qalign,dbalign] of rowsAll){ const key=`${name}\t${dbs}\t${dbe}\t${strand}\t${qalign}\t${dbalign}`; if(emitted.has(key)) continue; emitted.add(key); rowsOut.push({ qid, exp, tm, name, start: dbs, end: dbe, strand, ident, qalign, dbalign }); } } } return rowsOut; }

  T.parseFasta = parseFasta;
  T.parseQueries = parseQueries;
  T.parseFastaOrRaw = parseFastaOrRaw;
  T.tmblast = tmblast;
  T.tmblastRows = tmblastRows;
  window.TmBLAST = T;
})();
