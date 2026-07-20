#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {

	// Germline gene files (if available)
	ch_germlines = Channel
		.fromPath( "${params.germline_dir}/*.fasta", checkIfExists: false )
		.collect()
		.ifEmpty( [] )

	if ( params.pooled_fastq ) {

		// PCR barcode mode: demultiplex a single pooled FASTQ by plate + well barcodes.
		// All resulting reads are IgH-targeted, so chain classification is skipped.
		DEMULTIPLEX(
			Channel.fromPath( params.pooled_fastq ),
			Channel.fromPath( params.barcode_index )
		)

		ch_classified = DEMULTIPLEX.out
			.flatten()
			.filter { it.name ==~ /PB\d+_[A-H]\d+\.fastq\.gz/ }
			.map { fastq ->
				def cell_id = fastq.name.replaceAll( '\\.fastq\\.gz$', '' )
				tuple( cell_id, "heavy", fastq )
			}
			.filter { cell_id, chain, fastq ->
				fastq.countFastq() >= params.min_reads
			}

	} else {

		// ONT barcode directory mode: reads already demultiplexed by MinKNOW.
		ch_barcodes = Channel
			.fromPath( "${params.fastq_dir}/barcode*/", type: 'dir' )
			.map { dir -> tuple( dir.getName(), dir ) }

		// Load classification primers (separated by chain type)
		ch_class_primers = Channel
			.fromPath( params.primer_table )
			.splitCsv( header: true )
			.filter { row -> row.differentiating == "true" }

		heavy_primers = ch_class_primers
			.filter { row -> row.chain == 'heavy' }
			.map { row -> row.primer_seq }
			.collect()

		light_primers = ch_class_primers
			.filter { row -> row.chain == 'light' }
			.map { row -> row.primer_seq }
			.collect()

		// Stage 1: Read Processing - Merge reads per barcode
		MERGE_READS(
			ch_barcodes
		)

		// Stage 2: Classify reads by primer sequence
		CLASSIFY_BY_PRIMER(
			MERGE_READS.out
				.combine( heavy_primers )
				.combine( light_primers )
		)

		// Combine heavy and light outputs, filter empty files
		ch_classified = CLASSIFY_BY_PRIMER.out.heavy
			.mix( CLASSIFY_BY_PRIMER.out.light )
			.filter { barcode_id, chain, fastq ->
				fastq.countFastq() >= params.min_reads
			}

	}

	// Stage 3: Quality filter
	QUALITY_FILTER(
		ch_classified
	)

	if ( params.keep_primers ) {
		// Skip adapter trimming and pass reads directly to clustering
		ch_for_clustering = QUALITY_FILTER.out
	} else {
		// Detect and trim sequencing adapters
		FIND_ADAPTERS(
			QUALITY_FILTER.out
		)
		TRIM_PRIMERS(
			QUALITY_FILTER.out
				.join( FIND_ADAPTERS.out, by: [0, 1] )
		)
		ch_for_clustering = TRIM_PRIMERS.out
	}

	// Stage 4: Clustering & Consensus
	CONVERT_TO_FASTA(
		ch_for_clustering
	)

	CLUSTER_READS(
		CONVERT_TO_FASTA.out
	)

	// Stage 5: Annotation (conditional on germline files being available)
	if ( !params.skip_annotation ) {
		BUILD_IGBLAST_DB(
			ch_germlines
		)

		ANNOTATE_IGBLAST(
			BUILD_IGBLAST_DB.out,
			CLUSTER_READS.out
		)

		PARSE_ANNOTATIONS(
			ANNOTATE_IGBLAST.out
		)

		// Stage 6: Reporting
		// Extract only the annotations.tsv path (third element) before collecting
		COLLECT_STATS(
			PARSE_ANNOTATIONS.out.map { barcode_id, chain, annotations_tsv, cdr3_fasta -> annotations_tsv }.collect()
		)
	} else {
		// Skip annotation, just collect consensus stats
		// Extract only the path (third element) from the tuple before collecting
		COLLECT_CONSENSUS_STATS(
			CLUSTER_READS.out.map { barcode_id, chain, consensus_dir -> consensus_dir }.collect()
		)
	}

	if ( !params.skip_annotation ) {
		GENERATE_REPORT(
			COLLECT_STATS.out
		)
	} else {
		GENERATE_REPORT(
			COLLECT_CONSENSUS_STATS.out
		)
	}

}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
params.errorMode = params.debugmode == true ? 'terminate' : 'ignore'

params.demuxed_reads  = params.results + "/0_demuxed_reads"
params.merged_reads   = params.results + "/1_merged_reads"
params.classified_reads = params.results + "/2_classified_reads"
params.filtered_reads = params.results + "/3_filtered_reads"
params.consensus_seqs = params.results + "/4_consensus_sequences"
params.annotations = params.results + "/5_annotations"
params.reports = params.results + "/6_reports"
// --------------------------------------------------------------- //



// PROCESS SPECIFICATION
// --------------------------------------------------------------- //

process DEMULTIPLEX {

	publishDir "${params.demuxed_reads}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 4

	input:
	path pooled_fastq
	path barcode_index

	output:
	path "*.fastq.gz"

	script:
	"""
	minibar.py ${barcode_index} ${pooled_fastq} \
		-e ${params.barcode_error} \
		-F -P "" -T

	# Drop zero-length reads (Nextflow's FastqSplitter/.countFastq() treats
	# them as malformed records due to Groovy's empty-string falsiness),
	# then compress per-cell files; drop now-empty files (empty wells)
	shopt -s nullglob
	for f in *.fastq; do
		seqkit seq -m 1 "\$f" -o "\${f}.tmp" && mv "\${f}.tmp" "\$f"
		if [ -s "\$f" ]; then
			gzip "\$f"
		else
			rm -f "\$f"
		fi
	done
	"""

}

process MERGE_READS {

	tag { "${barcode_id}" }
	publishDir path: { "${params.merged_reads}/${barcode_id}" }, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 4

	input:
	tuple val(barcode_id), path(read_dir)

	output:
	tuple val(barcode_id), path("${barcode_id}_merged.fastq.gz")

	script:
	"""
	seqkit scat -j ${task.cpus} -f `realpath ${read_dir}` -o ${barcode_id}_merged.fastq.gz
	"""

}

process CLASSIFY_BY_PRIMER {

	tag { "${barcode_id}" }
	publishDir path: { "${params.classified_reads}/${barcode_id}" }, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 4

	input:
	tuple val(barcode_id), path(merged_reads), val(heavy_primer_list), val(light_primer_list)

	output:
	tuple val(barcode_id), val("heavy"), path("${barcode_id}_heavy.fastq.gz"), emit: heavy
	tuple val(barcode_id), val("light"), path("${barcode_id}_light.fastq.gz"), emit: light
	tuple val(barcode_id), val("unmatched"), path("${barcode_id}_unmatched.fastq.gz"), emit: unmatched

	script:
	def heavy_list = heavy_primer_list instanceof List ? heavy_primer_list : [heavy_primer_list]
	def light_list = light_primer_list instanceof List ? light_primer_list : [light_primer_list]
	def heavy_fasta = heavy_list.withIndex().collect { seq, i -> ">heavy_${i}\n${seq}" }.join('\n')
	def light_fasta = light_list.withIndex().collect { seq, i -> ">light_${i}\n${seq}" }.join('\n')
	"""
	# Create primer reference files in FASTA format
	cat << 'HEAVY_EOF' > heavy_primers.fasta
${heavy_fasta}
HEAVY_EOF

	cat << 'LIGHT_EOF' > light_primers.fasta
${light_fasta}
LIGHT_EOF

	# Step 1: Extract heavy chain reads (matching heavy primers)
	bbduk.sh in=`realpath ${merged_reads}` \
		outm=${barcode_id}_heavy.fastq.gz \
		outu=temp_not_heavy.fastq.gz \
		ref=heavy_primers.fasta \
		k=15 hdist=${params.primer_mismatch} \
		qin=33 threads=${task.cpus}

	# Step 2: From remaining reads, extract light chain reads
	bbduk.sh in=temp_not_heavy.fastq.gz \
		outm=${barcode_id}_light.fastq.gz \
		outu=${barcode_id}_unmatched.fastq.gz \
		ref=light_primers.fasta \
		k=15 hdist=${params.primer_mismatch} \
		qin=33 threads=${task.cpus}

	# Cleanup temp file
	rm -f temp_not_heavy.fastq.gz
	"""

}

process QUALITY_FILTER {

	tag { "${barcode_id}_${chain}" }
	publishDir path: { "${params.filtered_reads}/${barcode_id}" }, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 4

	input:
	tuple val(barcode_id), val(chain), path(reads)

	output:
	tuple val(barcode_id), val(chain), path("${barcode_id}_${chain}_filtered.fastq.gz")

	script:
	"""
	seqkit seq \
	--min-len ${params.min_len} \
	--max-len ${params.max_len} \
	--min-qual ${params.min_qual} \
	--validate-seq \
	--threads ${task.cpus} \
	${reads} \
	-o ${barcode_id}_${chain}_filtered.fastq.gz
	"""

}

process FIND_ADAPTERS {

	tag { "${barcode_id}_${chain}" }

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	tuple val(barcode_id), val(chain), path(reads)

	output:
	tuple val(barcode_id), val(chain), path("${barcode_id}_${chain}_adapters.fasta")

	script:
	"""
	bbmerge.sh in=`realpath ${reads}` outa="${barcode_id}_${chain}_adapters.fasta" ow qin=33
	"""

}

process TRIM_PRIMERS {

	tag { "${barcode_id}_${chain}" }
	publishDir path: { "${params.filtered_reads}/${barcode_id}" }, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 4

	input:
	tuple val(barcode_id), val(chain), path(reads), path(adapters)

	output:
	tuple val(barcode_id), val(chain), path("${barcode_id}_${chain}_trimmed.fastq.gz")

	script:
	"""
	# Check if adapters file has real sequences (not just N's)
	if grep -v '^>' ${adapters} | grep -q '[ACGT]'; then
		bbduk.sh in=`realpath ${reads}` out=${barcode_id}_${chain}_trimmed.fastq.gz \
		ref=`realpath ${adapters}` \
		ktrim=r k=19 mink=11 hdist=2 \
		minlength=${params.min_len} maxlength=${params.max_len} \
		qin=33 threads=${task.cpus}
	else
		# No real adapters found, just copy input to output
		cp `realpath ${reads}` ${barcode_id}_${chain}_trimmed.fastq.gz
	fi
	"""

}

process CONVERT_TO_FASTA {

	tag { "${barcode_id}_${chain}" }

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	tuple val(barcode_id), val(chain), path(reads)

	output:
	tuple val(barcode_id), val(chain), path("${barcode_id}_${chain}.fasta")

	script:
	"""
	seqkit fq2fa `realpath ${reads}` > ${barcode_id}_${chain}.fasta
	"""

}

process CLUSTER_READS {

	tag { "${barcode_id}_${chain}" }
	publishDir path: { "${params.consensus_seqs}/${barcode_id}" }, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 4

	input:
	tuple val(barcode_id), val(chain), path(fasta)

	output:
	tuple val(barcode_id), val(chain), path("${barcode_id}_${chain}_consensus")

	script:
	def amb_flag = params.use_ambiguous ? '-amb' : ''
	"""
	amplicon_sorter.py \
	-i ${fasta} \
	-o ${barcode_id}_${chain}_consensus \
	-min ${params.min_len} -max ${params.max_len} \
	-sg ${params.similar_genes} \
	-ss ${params.similar_species} \
	-sc ${params.similar_consensus} \
	-ldc ${params.length_diff_consensus} \
	${amb_flag} \
	-ar -maxr ${params.max_reads} -ra -np ${task.cpus}
	"""

}

process BUILD_IGBLAST_DB {

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	path fastas

	output:
	path "bovine_ig_db*"

	script:
	"""
	# Combine all germline fastas
	cat *.fasta > combined_germlines.fasta

	# Format for IgBLAST (remove gaps, standardize headers)
	sed 's/\\./-/g' combined_germlines.fasta | \
	awk '/^>/{print; next}{gsub(/\\./, ""); print}' > bovine_ig_db

	# Build BLAST database
	makeblastdb -parse_seqids -dbtype nucl -in bovine_ig_db
	"""

}

process ANNOTATE_IGBLAST {

	tag { "${barcode_id}_${chain}" }
	publishDir path: { "${params.annotations}/${barcode_id}" }, mode: 'copy', overwrite: true

	errorStrategy 'ignore'

	input:
	path db_files
	tuple val(barcode_id), val(chain), path(consensus_dir)

	output:
	tuple val(barcode_id), val(chain), path("${barcode_id}_${chain}_igblast.tsv"), path(consensus_dir)

	script:
	"""
	# Find consensus fasta files
	consensus_fasta=\$(find ${consensus_dir} -name "*.fasta" -o -name "*consensus*.fa" | head -1)

	if [ -z "\$consensus_fasta" ]; then
		# If no fasta found, create empty output
		echo "No consensus sequences found" > ${barcode_id}_${chain}_igblast.tsv
	else
		# Run IgBLAST with custom bovine database
		igblastn \
		-germline_db_V bovine_ig_db \
		-germline_db_J bovine_ig_db \
		-germline_db_D bovine_ig_db \
		-auxiliary_data optional_file/human_gl.aux \
		-query "\$consensus_fasta" \
		-outfmt "7 std qseq sseq" \
		-out ${barcode_id}_${chain}_igblast.tsv \
		|| echo "IgBLAST completed with warnings" > ${barcode_id}_${chain}_igblast.tsv
	fi
	"""

}

process PARSE_ANNOTATIONS {

	tag { "${barcode_id}_${chain}" }
	publishDir path: { "${params.annotations}/${barcode_id}" }, mode: 'copy', overwrite: true

	errorStrategy 'ignore'

	input:
	tuple val(barcode_id), val(chain), path(igblast_out), path(consensus_dir)

	output:
	tuple val(barcode_id), val(chain), path("${barcode_id}_${chain}_annotations.tsv"), path("${barcode_id}_${chain}_cdr3.fasta")

	script:
	"""
	parse_igblast.py \
	--input ${igblast_out} \
	--consensus_dir ${consensus_dir} \
	--chain ${chain} \
	--output_tsv ${barcode_id}_${chain}_annotations.tsv \
	--output_cdr3 ${barcode_id}_${chain}_cdr3.fasta
	"""

}

process COLLECT_STATS {

	publishDir params.reports, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	path annotations

	output:
	path "summary_stats.tsv"

	script:
	"""
	echo -e "barcode\tchain\tnum_sequences\tnum_unique_v\tnum_unique_j\tavg_cdr3_len" > summary_stats.tsv

	for f in *_annotations.tsv; do
		# Extract barcode and chain from filename (format: barcode_chain_annotations.tsv)
		basename=\$(basename "\$f" _annotations.tsv)
		barcode=\$(echo "\$basename" | rev | cut -d'_' -f2- | rev)
		chain=\$(echo "\$basename" | rev | cut -d'_' -f1 | rev)
		if [ -s "\$f" ]; then
			num_seqs=\$(tail -n +2 "\$f" | wc -l)
			echo -e "\${barcode}\t\${chain}\t\${num_seqs}\tNA\tNA\tNA" >> summary_stats.tsv
		fi
	done
	"""

}

process GENERATE_REPORT {

	publishDir params.reports, mode: 'copy', overwrite: true

	errorStrategy 'ignore'

	input:
	path stats

	output:
	path "*.pdf", optional: true
	path "report_summary.txt"

	script:
	"""
	echo "Bovine IgG Repertoire Analysis Report" > report_summary.txt
	echo "=====================================" >> report_summary.txt
	echo "" >> report_summary.txt
	cat ${stats} >> report_summary.txt
	echo "" >> report_summary.txt
	echo "Analysis completed: \$(date)" >> report_summary.txt
	"""

}

process COLLECT_CONSENSUS_STATS {

	publishDir params.reports, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	path consensus_dirs

	output:
	path "summary_stats.tsv"

	script:
	"""
	echo -e "barcode\tchain\tnum_consensus_sequences\ttotal_reads" > summary_stats.tsv

	for dir in */; do
		# Extract barcode and chain from directory name (format: barcode_chain_consensus)
		dirname=\$(basename "\$dir" _consensus)
		barcode=\$(echo "\$dirname" | rev | cut -d'_' -f2- | rev)
		chain=\$(echo "\$dirname" | rev | cut -d'_' -f1 | rev)
		if [ -d "\$dir" ]; then
			# Count consensus sequences
			num_seqs=\$(find "\$dir" -name "*.fasta" -exec grep -c "^>" {} + 2>/dev/null | awk -F: '{sum+=\$2} END {print sum}' || echo 0)
			echo -e "\${barcode}\t\${chain}\t\${num_seqs}\tNA" >> summary_stats.tsv
		fi
	done
	"""

}

// --------------------------------------------------------------- //
