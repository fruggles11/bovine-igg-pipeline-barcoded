#!/usr/bin/env python3
"""
Pull candidate bovine ultra-long CDRH3 clones out of vdj_summary.tsv.

A clone is flagged if either:
  - its CDR3 (junction_aa_length) is >= --min_cdr3_aa, or
  - it uses the V/D/J germline combination documented as the genetic
    signature of the bovine ultra-long CDRH3 "stalk and knob" architecture
    (default: IGHV1-7 + IGHD8-2 + IGHJ2-4), all three required.

v_call/d_call/j_call fields may hold multiple comma-separated, allele-suffixed
calls when IgBLAST couldn't resolve a single best hit (e.g.
"IGHV1-20*01,IGHV1-27*01"); matching is done at the gene level (allele suffix
stripped) against the set of all called genes for that field.
"""

import argparse
import csv
import os


def base_genes(call_field):
    """Split a comma-separated, allele-suffixed call field into base gene names."""
    if not call_field or call_field == "NA":
        return set()
    genes = set()
    for allele in call_field.split(','):
        allele = allele.strip()
        gene = allele.split('*')[0]
        if gene:
            genes.add(gene)
    return genes


def main():
    parser = argparse.ArgumentParser(description="Extract candidate ultra-long CDRH3 clones from vdj_summary.tsv")
    parser.add_argument('--vdj_summary', required=True)
    parser.add_argument('--fasta_dir', default='.', help="Directory containing *_majority.fasta files")
    parser.add_argument('--min_cdr3_aa', type=int, default=40)
    parser.add_argument('--v_gene', default='IGHV1-7')
    parser.add_argument('--d_gene', default='IGHD8-2')
    parser.add_argument('--j_gene', default='IGHJ2-4')
    parser.add_argument('--output_tsv', required=True)
    parser.add_argument('--output_fasta', required=True)
    args = parser.parse_args()

    matches = []
    with open(args.vdj_summary) as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames
        for row in reader:
            length_str = row.get('junction_aa_length', 'NA')
            try:
                length = int(length_str)
            except (ValueError, TypeError):
                length = None

            is_long = length is not None and length >= args.min_cdr3_aa

            v_genes = base_genes(row.get('v_call', ''))
            d_genes = base_genes(row.get('d_call', ''))
            j_genes = base_genes(row.get('j_call', ''))
            is_vdj_combo = (
                args.v_gene in v_genes
                and args.d_gene in d_genes
                and args.j_gene in j_genes
            )

            if is_long or is_vdj_combo:
                reasons = []
                if is_long:
                    reasons.append(f"cdr3_aa>={args.min_cdr3_aa}")
                if is_vdj_combo:
                    reasons.append(f"{args.v_gene}+{args.d_gene}+{args.j_gene}")
                row['match_reason'] = ",".join(reasons)
                matches.append(row)

    out_fields = fieldnames + ['match_reason']
    with open(args.output_tsv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=out_fields, delimiter='\t')
        writer.writeheader()
        for row in matches:
            writer.writerow(row)

    with open(args.output_fasta, 'w') as out_f:
        for row in matches:
            barcode = row['barcode']
            chain = row['chain']
            fasta_path = os.path.join(args.fasta_dir, f"{barcode}_{chain}_majority.fasta")
            if not os.path.exists(fasta_path):
                continue
            with open(fasta_path) as in_f:
                content = in_f.read().strip()
            if not content:
                continue
            lines = content.splitlines()
            seq = ''.join(lines[1:])
            out_f.write(
                f">{barcode}_{chain} v_call={row.get('v_call', 'NA')} "
                f"d_call={row.get('d_call', 'NA')} j_call={row.get('j_call', 'NA')} "
                f"cdr3_aa_len={row.get('junction_aa_length', 'NA')} reason={row['match_reason']}\n"
            )
            out_f.write(f"{seq}\n")

    print(f"Matched {len(matches)} clones -> {args.output_tsv}, {args.output_fasta}")


if __name__ == '__main__':
    main()
