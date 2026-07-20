#!/usr/bin/env python3
"""
Fill in missing junction/junction_aa calls in an AIRR TSV via J-gene anchor scanning.

IgBLAST's junction-detection heuristic silently leaves junction_aa blank for some
bovine heavy chain sequences (notably ultra-long CDR3H3, which exceed its internal
window) even though V/D/J genes were assigned. This backfills those rows by scanning
the row's own 'sequence' field for the conserved J-region Trp codon (via IGHJ germline
anchor matching) and walking back to the nearest upstream Cys codon, the same method
used in find_ultralong_cdrh3.py. Rows IgBLAST already annotated are left untouched.

Only applied to rows with an IGHV v_call (heavy chain) — light chain CDR3s are short
and are already reliably called by IgBLAST.
"""

import argparse
import csv
import re

_ACGT = frozenset('ACGT')

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def revcomp(seq):
    comp = str.maketrans('ACGTacgtNnRrYyKkMmSsWwBbDdHhVv',
                          'TGCAtgcaNnYyRrMmKkSsWwVvHhDdBb')
    return seq.translate(comp)[::-1]


def translate(nt_seq):
    aa = []
    for i in range(0, len(nt_seq) - 2, 3):
        codon = nt_seq[i:i + 3].upper()
        if 'N' in codon or len(codon) < 3:
            aa.append('X')
        else:
            aa.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(aa)


def parse_fasta(path):
    seqs = {}
    current_id, current_seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if current_id:
                    seqs[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id:
        seqs[current_id] = ''.join(current_seq)
    return seqs


def load_j_germlines(fasta_path):
    raw = parse_fasta(fasta_path)
    j_genes = []
    for header, seq in raw.items():
        name = header.split('|')[1] if '|' in header else header
        clean = re.sub(r'[.\-\s]', '', seq).upper()
        trp_pos = clean.find('TGG')
        if trp_pos == -1:
            continue
        j_genes.append((name, clean, trp_pos))
    return j_genes


def find_j_in_sequence(sequence, j_genes, anchor_len=20, min_identity=0.75):
    seq_u = sequence.upper()
    best = (None, None, 0.0)
    for j_name, j_seq, trp_pos in j_genes:
        anchor = j_seq[trp_pos: trp_pos + anchor_len]
        if len(anchor) < anchor_len:
            continue
        for i in range(len(seq_u) - len(anchor) + 1):
            window = seq_u[i: i + len(anchor)]
            matches = sum(1 for a, b in zip(anchor, window) if b in _ACGT and a == b)
            unambiguous = sum(1 for b in window if b in _ACGT)
            if unambiguous == 0:
                continue
            identity = matches / unambiguous
            if identity > best[2]:
                best = (i, j_name, identity)
    if best[2] >= min_identity:
        return best
    return (None, None, 0.0)


def find_cys_before_trp(sequence, trp_pos, max_cdr3_nt=600):
    seq_u = sequence.upper()
    for pos in range(trp_pos - 3, max(0, trp_pos - max_cdr3_nt) - 1, -3):
        codon = seq_u[pos: pos + 3]
        if codon in ('TGT', 'TGC'):
            return pos
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, help='AIRR TSV from MakeDb.py')
    parser.add_argument('--germline', required=True, help='IGHJ germline FASTA')
    parser.add_argument('--output', required=True, help='Output AIRR TSV')
    parser.add_argument('--min-identity', type=float, default=0.75,
                         help='Minimum J anchor match identity (default: 0.75)')
    args = parser.parse_args()

    j_genes = load_j_germlines(args.germline)

    with open(args.input) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    if 'junction_source' not in fieldnames:
        fieldnames.append('junction_source')

    filled = 0
    for row in rows:
        v_call = row.get('v_call', '') or ''
        jaa = row.get('junction_aa', '') or ''

        if jaa:
            row['junction_source'] = 'igblast'
            continue

        row['junction_source'] = ''
        if 'IGHV' not in v_call or not j_genes:
            continue

        sequence = row.get('sequence', '') or ''
        if not sequence:
            continue

        for test_seq in (sequence, revcomp(sequence)):
            trp_pos, j_name, identity = find_j_in_sequence(
                test_seq, j_genes, min_identity=args.min_identity
            )
            if trp_pos is None:
                continue
            cys_pos = find_cys_before_trp(test_seq, trp_pos)
            if cys_pos is None:
                continue

            junction_nt = test_seq[cys_pos: trp_pos + 3]
            row['junction'] = junction_nt
            row['junction_aa'] = translate(junction_nt)
            row['junction_length'] = str(len(junction_nt))
            if not row.get('j_call'):
                row['j_call'] = j_name
            row['junction_source'] = 'anchor'
            filled += 1
            break

    print(f'Filled {filled} missing junction call(s) via J-anchor fallback '
          f'(of {len(rows)} total rows)')

    with open(args.output, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t',
                                 extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()
