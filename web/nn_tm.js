// Port of Santalucia_NN_Tm.py (Thermo-Align) to JS for parity with Python

// Thermodynamic tables
export const DNA_NN_table = {
  'init': [0.2, -5.7], 'init_A/T': [2.2, 6.9], 'init_G/C': [0, 0],
  'init_oneG/C': [0, 0], 'init_allA/T': [0, 0], 'init_5T/A': [0, 0],
  'sym': [0, -1.4],
  'AA/TT': [-7.6, -21.3], 'AT/TA': [-7.2, -20.4], 'TA/AT': [-7.2, -20.4],
  'CA/GT': [-8.5, -22.7], 'GT/CA': [-8.4, -22.4], 'CT/GA': [-7.8, -21.0],
  'GA/CT': [-8.2, -22.2], 'CG/GC': [-10.6, -27.2], 'GC/CG': [-9.8, -24.4],
  'GG/CC': [-8.0, -19.0]
};

export const DNA_IMM_table = {
  'AG/TT': [1.0, 0.9], 'AT/TG': [-2.5, -8.3], 'CG/GT': [-4.1, -11.7],
  'CT/GG': [-2.8, -8.0], 'GG/CT': [3.3, 10.4], 'GG/TT': [5.8, 16.3],
  'GT/CG': [-4.4, -12.3], 'GT/TG': [4.1, 9.5], 'TG/AT': [-0.1, -1.7],
  'TG/GT': [-1.4, -6.2], 'TT/AG': [-1.3, -5.3],

  'AA/TG': [-0.6, -2.3], 'AG/TA': [-0.7, -2.3], 'CA/GG': [-0.7, -2.3],
  'CG/GA': [-4.0, -13.2], 'GA/CG': [-0.6, -1.0], 'GG/CA': [0.5, 3.2],
  'TA/AG': [0.7, 0.7], 'TG/AA': [3.0, 7.4],

  'AC/TT': [0.7, 0.2], 'AT/TC': [-1.2, -6.2], 'CC/GT': [-0.8, -4.5],
  'CT/GC': [-1.5, -6.1], 'GC/CT': [2.3, 5.4], 'GT/CC': [5.2, 13.5],
  'TC/AT': [1.2, 0.7], 'TT/AC': [1.0, 0.7],

  'AA/TC': [2.3, 4.6], 'AC/TA': [5.3, 14.6], 'CA/GC': [1.9, 3.7],
  'CC/GA': [0.6, -0.6], 'GA/CC': [5.2, 14.2], 'GC/CA': [-0.7, -3.8],
  'TA/AC': [3.4, 8.0], 'TC/AA': [7.6, 20.2],

  'AA/TA': [1.2, 1.7], 'CA/GA': [-0.9, -4.2], 'GA/CA': [-2.9, -9.8],
  'TA/AA': [4.7, 12.9], 'AC/TC': [0.0, -4.4], 'CC/GC': [-1.5, -7.2],
  'GC/CC': [3.6, 8.9], 'TC/AC': [6.1, 16.4], 'AG/TG': [-3.1, -9.5],
  'CG/GG': [-4.9, -15.3], 'GG/CG': [-6.0, -15.8], 'TG/AG': [1.6, 3.6],
  'AT/TT': [-2.7, -10.8], 'CT/GT': [-5.0, -15.8], 'GT/CT': [-2.2, -8.4],
  'TT/AT': [0.2, -1.5]
};

export const DNA_TMM_table = {
  'AA/TA': [-3.1, -7.8], 'TA/AA': [-2.5, -6.3],
  'CA/GA': [-4.3, -10.7], 'GA/CA': [-8.0, -22.5],
  'AC/TC': [-0.1, 0.5], 'TC/AC': [-0.7, -1.3],
  'CC/GC': [-2.1, -5.1], 'GC/CC': [-3.9, -10.6],
  'AG/TG': [-1.1, -2.1], 'TG/AG': [-1.1, -2.7],
  'CG/GG': [-3.8, -9.5], 'GG/CG': [-0.7, -19.2],
  'AT/TT': [-2.4, -6.5], 'TT/AT': [-3.2, -8.9],
  'CT/GT': [-6.1, -16.9], 'GT/CT': [-7.4, -21.2],
  'AA/TC': [-1.6, -4.0], 'AC/TA': [-1.8, -3.8], 'CA/GC': [-2.6, -5.9],
  'CC/GA': [-2.7, -6.0], 'GA/CC': [-5.0, -13.8], 'GC/CA': [-3.2, -7.1],
  'TA/AC': [-2.3, -5.9], 'TC/AA': [-2.7, -7.0],
  'AC/TT': [-0.9, -1.7], 'AT/TC': [-2.3, -6.3], 'CC/GT': [-3.2, -8.0],
  'CT/GC': [-3.9, -10.6], 'GC/CT': [-4.9, -13.5], 'GT/CC': [-3.0, -7.8],
  'TC/AT': [-2.5, -6.3], 'TT/AC': [-0.7, -1.2],
  'AA/TG': [-1.9, -4.4], 'AG/TA': [-2.5, -5.9], 'CA/GG': [-3.9, -9.6],
  'CG/GA': [-6.0, -15.5], 'GA/CG': [-4.3, -11.1], 'GG/CA': [-4.6, -11.4],
  'TA/AG': [-2.0, -4.7], 'TG/AA': [-2.4, -5.8],
  'AG/TT': [-3.2, -8.7], 'AT/TG': [-3.5, -9.4], 'CG/GT': [-3.8, -9.0],
  'CT/GG': [-6.6, -18.7], 'GG/CT': [-5.7, -15.9], 'GT/CG': [-5.9, -16.1],
  'TG/AT': [-3.9, -10.5], 'TT/AG': [-3.6, -9.8]
};

const R = 1.9872;

function complBase(b) {
  const map = { A: 'T', C: 'G', G: 'C', T: 'A', '-': '-', N: 'N' };
  return map[b] || 'N';
}

export function complement(seq) {
  return seq.toUpperCase().split('').map(complBase).join('');
}

export function seq_qc(seq) {
  const allowed = new Set(['A','C','G','T','N','-']);
  return String(seq).replace(/\s+/g, '').toUpperCase().split('').filter(c => allowed.has(c)).join('');
}

export function ion_correction(Na=0, K=0, Tris=0, Mg=0, dNTPs=0, seq_len=0) {
  let NaNum = Number(Na)||0, KNum = Number(K)||0, TrisNum = Number(Tris)||0, MgNum = Number(Mg)||0, dNTPsNum = Number(dNTPs)||0;
  const Monovalent_mmol = NaNum + KNum + (TrisNum/2);
  let Na_eq_mmol = dNTPsNum >= MgNum ? Monovalent_mmol : (Monovalent_mmol + (120 * Math.sqrt(Math.max(0, MgNum - dNTPsNum))));
  const Na_eq_mol = Na_eq_mmol/1000.0;
  const corr = 0.368 * (seq_len - 1) * Math.log(Na_eq_mol);
  return corr;
}

export function NN_Tm_JS({ seq, compl_seq, primer_conc=0, Na=0, K=0, Tris=0, Mg=0, dNTPs=0, ion_corr=false }) {
  if (!seq) throw new Error('Please provide an input sequence!');
  seq = seq_qc(seq);
  if (!compl_seq) compl_seq = complement(seq);
  compl_seq = seq_qc(compl_seq);

  let temp_seq = String(seq);
  let temp_compl_seq = String(compl_seq);
  let dH = 0.0; // kcal/mol
  let dS = 0.0; // cal/mol*K
  const dH_index = 0, dS_index = 1;

  // Terminal mismatches
  const left_tmm = temp_compl_seq.slice(0,2).split('').reverse().join('') + '/' + temp_seq.slice(0,2).split('').reverse().join('');
  if (DNA_TMM_table[left_tmm]) {
    dH += DNA_TMM_table[left_tmm][dH_index];
    dS += DNA_TMM_table[left_tmm][dS_index];
    temp_compl_seq = temp_compl_seq.slice(1);
    temp_seq = temp_seq.slice(1);
  }
  const right_tmm = temp_seq.slice(-2) + '/' + temp_compl_seq.slice(-2);
  if (DNA_TMM_table[right_tmm]) {
    dH += DNA_TMM_table[right_tmm][dH_index];
    dS += DNA_TMM_table[right_tmm][dS_index];
    temp_compl_seq = temp_compl_seq.slice(0, -1);
    temp_seq = temp_seq.slice(0, -1);
  }

  // General initiation
  dH += DNA_NN_table['init'][dH_index];
  dS += DNA_NN_table['init'][dS_index];

  // A/T terminal bonuses (consider terminal mismatches in counting)
  const terminals = seq[0] + '/' + compl_seq[0] + seq[seq.length-1] + '/' + compl_seq[compl_seq.length-1];
  const count_AT = (terminals.match(/A\/T/g)||[]).length + (terminals.match(/T\/A/g)||[]).length;
  dH += DNA_NN_table['init_A/T'][dH_index] * count_AT;
  dS += DNA_NN_table['init_A/T'][dS_index] * count_AT;

  // NN loop with internal mismatches
  for (let base_index = 0; base_index < temp_seq.length - 1; base_index++) {
    const NN = temp_seq.slice(base_index, base_index+2) + '/' + temp_compl_seq.slice(base_index, base_index+2);
    if (DNA_IMM_table[NN]) {
      dH += DNA_IMM_table[NN][dH_index];
      dS += DNA_IMM_table[NN][dS_index];
    } else if (DNA_IMM_table[reverseKey(NN)]) {
      const key = reverseKey(NN);
      dH += DNA_IMM_table[key][dH_index];
      dS += DNA_IMM_table[key][dS_index];
    } else if (DNA_NN_table[NN]) {
      dH += DNA_NN_table[NN][dH_index];
      dS += DNA_NN_table[NN][dS_index];
    } else if (DNA_NN_table[reverseKey(NN)]) {
      const key = reverseKey(NN);
      dH += DNA_NN_table[key][dH_index];
      dS += DNA_NN_table[key][dS_index];
    }
  }

  if (ion_corr) {
    const seq_len = seq.length;
    const cf = ion_correction(Na, K, Tris, Mg, dNTPs, seq_len);
    dS += cf;
  }

  const x = 4; // non-self-complementary
  const primerConcM = (Number(primer_conc)||0) / 1e9; // nM -> M
  const Tm = (1000.0 * dH) / (dS + (R * Math.log(primerConcM / x))) - 273.15;
  return Tm;
}

function reverseKey(NN) {
  // reverse the entire key string as in Python NN[::-1]
  return NN.split('').reverse().join('');
}

export function calcTmSanta(seq, complSeq, opts = {}) {
  const {
    primerConc = 100, Na = 0, K = 50, Tris = 10, Mg = 1.5, dNTPs = 0.2,
    ionCorr = true,
  } = opts || {};
  return NN_Tm_JS({ seq, compl_seq: complSeq, primer_conc: primerConc, Na, K, Tris, Mg, dNTPs, ion_corr: ionCorr });
}
