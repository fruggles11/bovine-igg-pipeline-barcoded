# Bovine IgG Repertoire Analysis Pipeline

A Nextflow pipeline for analyzing bovine immunoglobulin heavy chain sequences from Oxford Nanopore amplicon data. Supports two input modes:

- **PCR barcode mode** — single pooled ONT run demultiplexed by plate + well PCR barcodes (up to 11 plates × 96 wells = 1,056 single cells per run)
- **ONT barcode mode** — standard MinKNOW-demultiplexed barcode directories (original behavior)

## Requirements

- [Nextflow](https://www.nextflow.io/) (v21.04+)
- [Docker](https://www.docker.com/) or [Singularity/Apptainer](https://apptainer.org/)
- [minibar](https://github.com/calacademy-research/minibar) (included in the Docker image)

## Building the Docker Image

The pipeline uses a custom Docker image (`fruggles11/bovine-igg-pipeline-barcoded:latest`) that includes all required tools. Build it on the machine where you run the pipeline (e.g. the Mac Studio):

```bash
git clone https://github.com/fruggles11/bovine-igg-pipeline-barcoded.git
cd bovine-igg-pipeline-barcoded
docker build -t fruggles11/bovine-igg-pipeline-barcoded:latest .
```

To make the image available on Docker Hub so you don't need to rebuild on every machine:

```bash
docker push fruggles11/bovine-igg-pipeline-barcoded:latest
```

You only need to rebuild the image if the `Dockerfile` or `environment.yml` changes.

## Quick Start

### PCR barcode mode (single-cell)

Sequence all cells in a single undemultiplexed ONT run and demultiplex by the plate + well PCR barcodes:

```bash
nextflow run fruggles11/bovine-igg-pipeline-barcoded \
  --pooled_fastq /path/to/your_run.fastq.gz \
  --keep_primers true
```

Use `--keep_primers true` to retain the inner gene-specific primer sequences on reads after demultiplexing. Omit it (or set `false`) to have adapter trimming run instead.

### ONT barcode mode (original)

Reads already demultiplexed by MinKNOW into barcode directories:

```bash
nextflow run fruggles11/bovine-igg-pipeline-barcoded \
  --fastq_dir /path/to/fastq_pass
```

## Input Data

### PCR barcode mode

A single gzipped (or uncompressed) FASTQ file from an ONT run with MinKNOW barcoding turned **off**. All 960–1,056 cells are pooled in one file; the pipeline demultiplexes them using the inline plate and well PCR barcodes.

The barcode index (`resources/pcr_barcodes_minibar.tsv`) is pre-populated with the full 11-plate × 96-well set. Each forward primer carries a 24 bp plate barcode and each reverse primer carries a 24 bp well barcode (minimum Levenshtein distance of 8 between any two barcodes).

### ONT barcode mode

Basecalled FASTQ files organized into barcode subdirectories:

```
fastq_pass/
├── barcode01/
│   ├── file1.fastq.gz
│   └── file2.fastq.gz
├── barcode02/
│   └── ...
└── barcode03/
    └── ...
```

All directories matching `barcode*/` are auto-detected. Reads are classified as heavy or light chain by primer sequence matching.

## Pipeline Steps

### PCR barcode mode

1. **Demultiplex** — minibar splits the pooled FASTQ into per-cell files by plate + well barcode, trimming barcodes from reads
2. **Quality Filter** — filter by length and quality score
3. **Trim Adapters** — detect and remove sequencing adapters
4. **Cluster Reads** — cluster reads into consensus sequences with amplicon_sorter
5. **Annotate** (optional) — V/D/J gene assignment with IgBLAST
6. **Report** — summary statistics

### ONT barcode mode

1. **Merge Reads** — concatenate reads within each barcode directory
2. **Classify by Primer** — sort reads into heavy/light chain
3. **Quality Filter** — filter by length and quality score
4. **Trim Adapters** — detect and remove sequencing adapters
5. **Cluster Reads** — cluster reads into consensus sequences with amplicon_sorter
6. **Annotate** (optional) — V/D/J gene assignment with IgBLAST
7. **Report** — summary statistics

## Parameters

### PCR barcode mode

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--pooled_fastq` | `""` | Path to pooled undemultiplexed FASTQ; set this to activate PCR barcode mode |
| `--barcode_index` | `resources/pcr_barcodes_minibar.tsv` | minibar index file (SampleID / Barcode1 / Barcode2) |
| `--barcode_error` | `2` | Edit distance allowed per barcode during demultiplexing |

### ONT barcode mode

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fastq_dir` | `""` | Path to directory containing barcode subdirectories |
| `--primer_table` | `resources/bovine_primers.csv` | CSV with primer sequences for chain classification |
| `--primer_mismatch` | `2` | Allowed mismatches when matching primers (max 3) |

### Shared parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_len` | `400` | Minimum amplicon length (bp) |
| `--max_len` | `1000` | Maximum amplicon length (bp) |
| `--min_qual` | `10` | Minimum Phred quality score |
| `--min_reads` | `100` | Minimum reads required per cell/barcode |
| `--similar_genes` | `85` | amplicon_sorter initial grouping threshold (%) |
| `--similar_species` | `90` | amplicon_sorter sequence addition threshold (%) |
| `--similar_consensus` | `95` | amplicon_sorter consensus merging threshold (%) |
| `--length_diff_consensus` | `10` | amplicon_sorter length variance allowance (%) |
| `--keep_primers` | `false` | Skip adapter trimming and retain inner primer sequences on reads |
| `--skip_annotation` | `false` | Skip repertoire analysis (V(D)J annotation, clonal assignment, diversity analysis) |
| `--clone_threshold` | `0.15` | Junction hamming distance threshold for clone definition |
| `--results` | `./results` | Output directory |

## Output

### PCR barcode mode

```
results/
├── 0_demuxed_reads/
│   ├── PB01_A01.fastq.gz    # One file per cell (plate_well)
│   ├── PB01_B01.fastq.gz
│   └── ...
├── 3_filtered_reads/
├── 4_consensus_sequences/
│   └── PB01_A01/
│       └── PB01_A01_heavy_consensus/
├── 5_majority_consensus/    # one dominant-cluster FASTA per cell, for e.g. Geneious import
└── 6_repertoire_analysis/   # if --skip_annotation false -- see below
```

### ONT barcode mode

```
results/
├── 1_merged_reads/
├── 2_classified_reads/
│   └── barcode01/
│       ├── barcode01_heavy.fastq.gz
│       ├── barcode01_light.fastq.gz
│       └── barcode01_unmatched.fastq.gz
├── 3_filtered_reads/
├── 4_consensus_sequences/
├── 5_majority_consensus/    # one dominant-cluster FASTA per cell, for e.g. Geneious import
└── 6_repertoire_analysis/   # if --skip_annotation false -- see below
```

### `6_repertoire_analysis/` in detail

```
6_repertoire_analysis/
├── igblast_db/               # bovine V/D/J BLAST databases built from --germline_dir
├── igblast/                  # raw per-cell IgBLAST output
├── airr/                     # per-cell AIRR-format TSVs (Change-O MakeDb.py)
├── airr_junction_filled/     # same, with missing CDR3/junction calls backfilled
├── filtered/                 # productivity-filtered AIRR TSVs
├── clones/                   # per-cell clone assignments (Change-O DefineClones.py)
├── reports/
│   └── vdj_summary.tsv       # combined V/D/J calls, one row per cell, with a
│                              # junction_source column (igblast vs anchor -- see below)
├── diversity/
│   ├── stats/                # CSV files: basic stats, gene usage, clone sizes, diversity indices
│   └── plots/                # PDF/PNG: CDR3 length, gene usage, clone size, diversity/rarefaction curves
└── repertoire_report.html
```

## Setting Up Repertoire Analysis (Optional)

Requires bovine germline gene FASTA files from IMGT — the same files used for the pipeline's other stages. Download and save them in `resources/germlines/` (or wherever `--germline_dir` points):

1. Go to [IMGT/GENE-DB](https://www.imgt.org/genedb/)
2. For each segment, select Species: *Bos taurus*, Functionality: functional, and export as FASTA
3. Save as: `bovine_IGHV.fasta`, `bovine_IGHD.fasta`, `bovine_IGHJ.fasta`, `bovine_IGKV.fasta`, `bovine_IGKJ.fasta`, `bovine_IGLV.fasta`, `bovine_IGLJ.fasta`

Run with `--skip_annotation true` if germline files are not available.

This stage uses the [Immcantation](https://immcantation.readthedocs.io/) framework (IgBLAST + Change-O + Alakazam) for V(D)J annotation, clonal assignment, and diversity analysis, via the `immcantation/suite:4.5.0` container (pulled automatically). Bovine has no official IgBLAST auxiliary data, so IgBLAST's own junction/CDR3 calling silently misses some sequences (notably bovine's hallmark ultra-long CDR3H3, which can exceed its internal detection window); an anchor-based fallback backfills those gaps by scanning for the conserved J-region motif against the real IGHJ germline. Each row in `vdj_summary.tsv` is tagged with `junction_source` (`igblast` or `anchor`) so you always know which method produced a given call.

This logic is maintained in [bovine-repertoire-analysis](https://github.com/fruggles11/bovine-repertoire-analysis), which also remains usable standalone — e.g. to reprocess consensus output from an older run, or to run its ultra-long CDR3H3 filter separately.

## Bovine IgG Considerations

- **Ultralong CDR3H**: Cattle can have CDR3H regions up to 70+ amino acids (vs ~15 for humans)
- **Limited V gene usage**: Cattle primarily use IGHV1-7 for heavy chains
- **Clustering thresholds**: Tuned for ONT error rates (~5–10%) and bovine IgH length variability

## License

MIT License
