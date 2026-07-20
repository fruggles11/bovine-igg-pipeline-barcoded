#!/usr/bin/env python3
"""
Parse IgBLAST output and extract V/D/J gene assignments and CDR3 sequences.
Designed for bovine immunoglobulin analysis.
"""

import argparse
import os
import re
from pathlib import Path


def parse_igblast_output(igblast_file):
    """Parse IgBLAST tabular output and extract annotations."""
    annotations = []

    if not os.path.exists(igblast_file) or os.path.getsize(igblast_file) == 0:
        return annotations

    current_query = None
    current_hits = []

    with open(igblast_file, 'r') as f:
        for line in f:
            line = line.strip()

            # Skip comments and empty lines
            if line.startswith('#') or not line:
                if 'Query:' in line:
                    # Save previous query results
                    if current_query and current_hits:
                        annotations.append({
                            'sequence_id': current_query,
                            'hits': current_hits
                        })
                    current_query = line.split('Query:')[-1].strip()
                    current_hits = []
                continue

            # Parse hit table rows. IgBLAST's "-outfmt 7 std qseq sseq" hit
            # table prepends a chain-type column (V/D/J) ahead of the
            # documented "Fields:" list, and inserts an extra "gaps" column
            # between "gap opens" and "q. start" -- so the real layout is:
            # hit_type, query id, subject id, % identity, alignment length,
            # mismatches, gap opens, gaps, q.start, q.end, s.start, s.end,
            # evalue, bit score, query seq, subject seq.
            fields = line.split('\t')
            if len(fields) >= 14:
                hit_type = fields[0]
                hit = {
                    'query_id': fields[1],
                    'subject_id': fields[2],
                    'identity': float(fields[3]) if fields[3] else 0,
                    'alignment_length': int(fields[4]) if fields[4] else 0,
                    'mismatches': int(fields[5]) if fields[5] else 0,
                    'gap_opens': int(fields[6]) if fields[6] else 0,
                    'q_start': int(fields[8]) if fields[8] else 0,
                    'q_end': int(fields[9]) if fields[9] else 0,
                    's_start': int(fields[10]) if fields[10] else 0,
                    's_end': int(fields[11]) if fields[11] else 0,
                    'evalue': float(fields[12]) if fields[12] else 0,
                    'bit_score': float(fields[13]) if fields[13] else 0,
                }

                # IgBLAST's hit table labels each row's chain type directly
                # (see "# Hit table (the first field indicates the chain
                # type of the hit)"), which is authoritative -- no need to
                # guess from the subject id.
                hit['gene_type'] = hit_type if hit_type in ('V', 'D', 'J') else 'unknown'

                if not current_query:
                    current_query = hit['query_id']
                current_hits.append(hit)

    # Don't forget the last query
    if current_query and current_hits:
        annotations.append({
            'sequence_id': current_query,
            'hits': current_hits
        })

    return annotations


def extract_best_vdj(annotations):
    """Extract best V, D, J gene assignments for each sequence."""
    results = []

    for annot in annotations:
        seq_id = annot['sequence_id']
        hits = annot['hits']

        best_v = None
        best_d = None
        best_j = None

        for hit in hits:
            gene_type = hit.get('gene_type', 'unknown')

            if gene_type == 'V':
                if best_v is None or hit['bit_score'] > best_v['bit_score']:
                    best_v = hit
            elif gene_type == 'D':
                if best_d is None or hit['bit_score'] > best_d['bit_score']:
                    best_d = hit
            elif gene_type == 'J':
                if best_j is None or hit['bit_score'] > best_j['bit_score']:
                    best_j = hit

        result = {
            'sequence_id': seq_id,
            'v_gene': best_v['subject_id'] if best_v else 'NA',
            'v_identity': best_v['identity'] if best_v else 'NA',
            'd_gene': best_d['subject_id'] if best_d else 'NA',
            'd_identity': best_d['identity'] if best_d else 'NA',
            'j_gene': best_j['subject_id'] if best_j else 'NA',
            'j_identity': best_j['identity'] if best_j else 'NA',
        }

        results.append(result)

    return results


def read_consensus_sequences(consensus_dir):
    """Read consensus sequences from amplicon_sorter output directory."""
    sequences = {}

    consensus_path = Path(consensus_dir)

    # Find fasta files in the directory
    fasta_files = list(consensus_path.glob('*.fasta')) + list(consensus_path.glob('*.fa'))

    for fasta_file in fasta_files:
        with open(fasta_file, 'r') as f:
            current_id = None
            current_seq = []

            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id and current_seq:
                        sequences[current_id] = ''.join(current_seq)
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)

            if current_id and current_seq:
                sequences[current_id] = ''.join(current_seq)

    return sequences


def extract_cdr3_simple(sequence, v_end=None, j_start=None):
    """
    Simple CDR3 extraction based on conserved motifs.
    For bovine, look for conserved Cys (TGT/TGC) before CDR3 and Trp (TGG) or Phe (TTC/TTT) after.
    """
    if not sequence:
        return None

    sequence = sequence.upper()

    # Look for conserved Cys codon (typically marks start of CDR3)
    cys_pattern = r'TG[TC]'
    # Look for conserved Trp-Gly (TGGGG) or Phe-Gly (TT[TC]GG) that marks J region
    j_pattern = r'(TGGGG|TT[TC]GG)'

    cys_matches = list(re.finditer(cys_pattern, sequence))
    j_matches = list(re.finditer(j_pattern, sequence))

    if not cys_matches or not j_matches:
        return None

    # Find the last Cys before the J motif
    j_pos = j_matches[0].start()

    cdr3_start = None
    for match in reversed(cys_matches):
        if match.start() < j_pos:
            cdr3_start = match.start()
            break

    if cdr3_start is None:
        return None

    # CDR3 extends from after Cys to include the Trp/Phe
    cdr3_end = j_pos + 3  # Include the first codon of J motif

    cdr3_seq = sequence[cdr3_start:cdr3_end]

    # Sanity check - CDR3 should be reasonable length
    # Bovine can have ultralong CDR3 (up to 200+ nt), but also normal ones
    if len(cdr3_seq) < 27 or len(cdr3_seq) > 250:
        return None

    return cdr3_seq


def main():
    parser = argparse.ArgumentParser(description='Parse IgBLAST output for bovine IgG')
    parser.add_argument('--input', required=True, help='IgBLAST output file')
    parser.add_argument('--consensus_dir', required=True, help='Directory with consensus sequences')
    parser.add_argument('--chain', required=True, choices=['heavy', 'light'], help='Chain type')
    parser.add_argument('--output_tsv', required=True, help='Output TSV file with annotations')
    parser.add_argument('--output_cdr3', required=True, help='Output FASTA file with CDR3 sequences')

    args = parser.parse_args()

    # Parse IgBLAST output
    annotations = parse_igblast_output(args.input)
    vdj_results = extract_best_vdj(annotations)

    # Read consensus sequences
    sequences = read_consensus_sequences(args.consensus_dir)

    # Write annotation TSV
    with open(args.output_tsv, 'w') as f:
        header = ['sequence_id', 'chain', 'v_gene', 'v_identity', 'd_gene', 'd_identity', 'j_gene', 'j_identity', 'cdr3_length']
        f.write('\t'.join(header) + '\n')

        for result in vdj_results:
            seq_id = result['sequence_id']
            seq = sequences.get(seq_id, '')
            cdr3 = extract_cdr3_simple(seq) if seq else None
            cdr3_len = len(cdr3) if cdr3 else 'NA'

            row = [
                seq_id,
                args.chain,
                result['v_gene'],
                str(result['v_identity']),
                result['d_gene'],
                str(result['d_identity']),
                result['j_gene'],
                str(result['j_identity']),
                str(cdr3_len)
            ]
            f.write('\t'.join(row) + '\n')

    # Write CDR3 FASTA
    with open(args.output_cdr3, 'w') as f:
        for result in vdj_results:
            seq_id = result['sequence_id']
            seq = sequences.get(seq_id, '')
            cdr3 = extract_cdr3_simple(seq) if seq else None

            if cdr3:
                f.write(f'>{seq_id}_CDR3 v_gene={result["v_gene"]} j_gene={result["j_gene"]}\n')
                f.write(f'{cdr3}\n')

    # If no IgBLAST results, still output sequences with basic info
    if not vdj_results and sequences:
        with open(args.output_tsv, 'w') as f:
            header = ['sequence_id', 'chain', 'v_gene', 'v_identity', 'd_gene', 'd_identity', 'j_gene', 'j_identity', 'cdr3_length']
            f.write('\t'.join(header) + '\n')

            for seq_id, seq in sequences.items():
                cdr3 = extract_cdr3_simple(seq)
                cdr3_len = len(cdr3) if cdr3 else 'NA'

                row = [seq_id, args.chain, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', str(cdr3_len)]
                f.write('\t'.join(row) + '\n')

        with open(args.output_cdr3, 'w') as f:
            for seq_id, seq in sequences.items():
                cdr3 = extract_cdr3_simple(seq)
                if cdr3:
                    f.write(f'>{seq_id}_CDR3\n')
                    f.write(f'{cdr3}\n')

    print(f"Parsed {len(vdj_results)} sequences with annotations")
    print(f"Found {len(sequences)} consensus sequences")


if __name__ == '__main__':
    main()
