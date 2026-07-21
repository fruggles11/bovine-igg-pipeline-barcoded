#!/usr/bin/env python3
"""
Combine per-barcode AIRR TSVs into a single one-row-per-sequence VDJ summary,
matching the vdj_summary.tsv shape used by bovine-igg-pipeline's SUMMARIZE_VDJ
(barcode, chain, v_call, d_call, j_call, productive, v_identity, cdr3_aa,
junction_aa_length), plus junction_source to show whether the CDR3 call came
from IgBLAST directly or the J-anchor fallback.
"""

import argparse
import csv
import re

SAMPLE_SUFFIX_RE = re.compile(r'^(.*)_(heavy|light)$')


def sample_to_barcode_chain(sample_id):
    m = SAMPLE_SUFFIX_RE.match(sample_id)
    if m:
        return m.group(1), m.group(2)
    return sample_id, 'NA'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', nargs='+', required=True,
                         help='Per-sample AIRR TSV files (named <sample_id>_productive.tsv)')
    parser.add_argument('--output', required=True, help='Combined output TSV')
    args = parser.parse_args()

    cols = ['barcode', 'chain', 'v_call', 'd_call', 'j_call',
            'productive', 'v_identity', 'cdr3_aa', 'junction_aa_length', 'junction_source']

    rows_out = []
    for path in sorted(args.inputs):
        sample_id = path.split('/')[-1]
        sample_id = re.sub(r'_productive\.tsv$', '', sample_id)
        sample_id = re.sub(r'_db-pass\.tsv$', '', sample_id)
        barcode, chain = sample_to_barcode_chain(sample_id)

        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                jaa = row.get('junction_aa', '') or ''
                rows_out.append({
                    'barcode':            barcode,
                    'chain':              chain,
                    'v_call':             row.get('v_call', '') or 'NA',
                    'd_call':             row.get('d_call', '') or 'NA',
                    'j_call':             row.get('j_call', '') or 'NA',
                    'productive':         row.get('productive', '') or 'NA',
                    'v_identity':         row.get('v_identity', '') or 'NA',
                    'cdr3_aa':            jaa if jaa else 'NA',
                    'junction_aa_length': str(len(jaa)) if jaa else 'NA',
                    'junction_source':    row.get('junction_source', '') or 'NA',
                })

    with open(args.output, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows_out)

    print(f'Wrote {len(rows_out)} rows to {args.output}')


if __name__ == '__main__':
    main()
