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

	// Stage 5: Extract majority consensus sequence per barcode/chain (for
	// easy import into tools like Geneious)
	EXTRACT_MAJORITY_CONSENSUS(
		CLUSTER_READS.out
	)

	// Stage 6: Repertoire analysis (V(D)J annotation, clonal assignment,
	// and diversity analysis via the Immcantation framework -- IgBLAST +
	// Change-O + Alakazam). Runs on the majority consensus sequence per
	// cell. Requires germline files; skip with --skip_annotation true.
	if ( !params.skip_annotation ) {

		if ( !params.germline_dir ) {
			error """
			=====================================================================
			ERROR: Bovine germline database is required for repertoire analysis
			=====================================================================
			Download germline FASTA files from IMGT/GENE-DB and pass their
			directory via --germline_dir, or skip this stage entirely with
			--skip_annotation true.
			=====================================================================
			""".stripIndent()
		}

		// Skip cells EXTRACT_MAJORITY_CONSENSUS couldn't build a confident
		// consensus for (empty placeholder file) -- MakeDb.py crashes hard
		// on an empty input instead of skipping it.
		ch_repertoire_input = EXTRACT_MAJORITY_CONSENSUS.out
			.filter { barcode_id, chain, fasta -> fasta.size() > 0 }
			.map { barcode_id, chain, fasta -> tuple( "${barcode_id}_${chain}", fasta ) }

		BUILD_IGBLAST_DB(
			ch_germlines
		)

		IGBLAST_ANNOTATION(
			ch_repertoire_input,
			BUILD_IGBLAST_DB.out.db,
			BUILD_IGBLAST_DB.out.aux
		)

		MAKEDB(
			IGBLAST_ANNOTATION.out,
			ch_germlines
		)

		// Backfill CDR3/junction calls IgBLAST's heuristic missed (e.g.
		// bovine ultra-long CDR3H3, which exceed its internal detection
		// window)
		INFER_JUNCTION(
			MAKEDB.out,
			ch_germlines,
			file( "${projectDir}/bin/infer_missing_junction.py" )
		)

		FILTER_PRODUCTIVE(
			INFER_JUNCTION.out
		)

		DEFINE_CLONES(
			FILTER_PRODUCTIVE.out
		)

		// Combine per-cell AIRR TSVs into one VDJ summary TSV
		SUMMARIZE_VDJ(
			FILTER_PRODUCTIVE.out
				.map { sample_id, tsv -> tsv }
				.collect(),
			file( "${projectDir}/bin/summarize_vdj.py" )
		)

		// Pull out candidate ultra-long CDRH3 clones: CDR3 >= 40 aa, or the
		// IGHV1-7 + IGHD8-2 + IGHJ2-4 germline combination documented as
		// bovine's ultra-long "stalk and knob" genetic signature
		EXTRACT_ULTRALONG_CLONES(
			SUMMARIZE_VDJ.out,
			EXTRACT_MAJORITY_CONSENSUS.out
				.filter { barcode_id, chain, fasta -> fasta.size() > 0 }
				.map { barcode_id, chain, fasta -> fasta }
				.collect(),
			file( "${projectDir}/bin/extract_ultralong_clones.py" )
		)

		DIVERSITY_ANALYSIS(
			DEFINE_CLONES.out.collect()
		)

		GENERATE_REPORT(
			DIVERSITY_ANALYSIS.out.plots,
			DIVERSITY_ANALYSIS.out.stats
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
params.majority_consensus = params.results + "/5_majority_consensus"
params.repertoire_results = params.results + "/6_repertoire_analysis"
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

	// amplicon_sorter's -np multiprocessing gives no measurable speedup for
	// this pipeline's per-cell read counts (verified: -np 1 vs -np 4 on a
	// real 300-read cell took identical wall time, ~3s of actual CPU work
	// against ~70s wall time either way -- the runtime is dominated by
	// fixed per-iteration overhead, not parallelizable computation).
	// Reserving 4 cpus per task for no benefit was capping Nextflow to
	// ~2 concurrent CLUSTER_READS tasks instead of ~10.
	cpus 1

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

process EXTRACT_MAJORITY_CONSENSUS {

	tag { "${barcode_id}_${chain}" }
	publishDir "${params.majority_consensus}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	tuple val(barcode_id), val(chain), path(consensus_dir)

	output:
	tuple val(barcode_id), val(chain), path("${barcode_id}_${chain}_majority.fasta")

	script:
	"""
	# Find the top-level *_consensussequences.fasta (exclude group-specific *_N_consensussequences.fasta)
	consensus_file=\$(find -L ${consensus_dir} -name "*_consensussequences.fasta" \
		| grep -vE '_[0-9]+_consensussequences\\.fasta\$' | head -1)

	if [ -z "\$consensus_file" ]; then
		touch ${barcode_id}_${chain}_majority.fasta
		exit 0
	fi

	# Read counts are encoded in the header as (N) e.g. >consensus_barcode01_heavy_0_0(581)
	# Find the header with the highest read count
	grep "^>" "\$consensus_file" > headers.txt
	best_header=""
	best_count=0
	while IFS= read -r header; do
		count=\$(echo "\$header" | grep -oE '\\([0-9]+\\)' | tr -d '()')
		count=\${count:-0}
		if [ "\$count" -gt "\$best_count" ]; then
			best_count="\$count"
			best_header="\$header"
		fi
	done < headers.txt

	if [ -z "\$best_header" ]; then
		touch ${barcode_id}_${chain}_majority.fasta
		exit 0
	fi

	# Extract the sequence for that header and write with an informative name
	awk -v target="\$best_header" \
		-v new_header=">${barcode_id} chain=${chain} majority_consensus reads=\$best_count" '
		\$0 == target { found=1; next }
		/^>/ { if (found) exit }
		found { seq = seq \$0 }
		END { if (seq != "") { print new_header; print seq } }
	' "\$consensus_file" > ${barcode_id}_${chain}_majority.fasta
	"""

}

// The following processes are ported from bovine-repertoire-analysis
// (https://github.com/fruggles11/bovine-repertoire-analysis), which remains
// available standalone for reprocessing older consensus output or running
// the ultra-long CDR3H3 filter. They use the Immcantation suite's own
// container rather than this pipeline's image, and that container only
// ships an amd64 build, hence the explicit --platform on each one.
process BUILD_IGBLAST_DB {

	container 'immcantation/suite:4.5.0'
	containerOptions '--platform linux/amd64'

	publishDir "${params.repertoire_results}/igblast_db", mode: 'copy'

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	path germlines

	output:
	path "database/*", emit: db
	path "internal_data/*", emit: aux

	script:
	"""
	mkdir -p database internal_data

	# Separate V, D, J genes based on filename patterns (case-insensitive,
	# searched recursively with -L to follow Nextflow's symlinked path inputs)
	find -L . -type f \\( -iname "*IGHV*" -o -iname "*IGKV*" -o -iname "*IGLV*" \\) | while read f; do
		cat "\$f" >> database/bovine_V.fasta 2>/dev/null
	done
	find -L . -type f -iname "*IGHD*" | while read f; do
		cat "\$f" >> database/bovine_D.fasta 2>/dev/null
	done
	find -L . -type f \\( -iname "*IGHJ*" -o -iname "*IGKJ*" -o -iname "*IGLJ*" \\) | while read f; do
		cat "\$f" >> database/bovine_J.fasta 2>/dev/null
	done

	touch database/bovine_V.fasta database/bovine_D.fasta database/bovine_J.fasta

	cd database
	for f in bovine_*.fasta; do
		if [[ -s "\$f" ]]; then
			makeblastdb -parse_seqids -dbtype nucl -in "\$f"
		fi
	done
	cd ..

	# Copy internal data from the container's IGDATA if available
	if [[ -n "\${IGDATA:-}" ]] && [[ -d "\${IGDATA}/internal_data" ]]; then
		cp -r "\${IGDATA}/internal_data/"* internal_data/ 2>/dev/null || true
	fi

	# Bovine has no real IgBLAST auxiliary data (frame anchor positions per
	# J-gene name) -- this placeholder is why INFER_JUNCTION exists
	# downstream to backfill the junction calls IgBLAST can't make without it.
	mkdir -p internal_data/bovine
	cat > internal_data/bovine/bovine_gl.aux << 'AUXFILE'
# Bovine germline auxiliary data
# Frame information for bovine IG genes
AUXFILE
	"""

}

process IGBLAST_ANNOTATION {

	container 'immcantation/suite:4.5.0'
	containerOptions '--platform linux/amd64'

	tag { "${sample_id}" }
	publishDir "${params.repertoire_results}/igblast", mode: 'copy'

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 4

	input:
	tuple val(sample_id), path(fasta)
	path db
	path aux

	output:
	tuple val(sample_id), path("${sample_id}_igblast.fmt7"), path(fasta)

	script:
	"""
	mkdir -p igblast_data/database
	for f in ${db}; do
		cp "\$f" igblast_data/database/
	done

	if [[ -n "\${IGDATA:-}" ]] && [[ -d "\${IGDATA}/internal_data" ]]; then
		cp -r "\${IGDATA}/internal_data" igblast_data/
	else
		for path in /usr/local/share/igblast/internal_data /usr/share/igblast/internal_data; do
			if [[ -d "\$path" ]]; then
				cp -r "\$path" igblast_data/
				break
			fi
		done
	fi

	export IGDATA="\$(pwd)/igblast_data"

	# -organism human is required for IgBLAST's internal_data validation;
	# the actual annotation uses our bovine databases via -germline_db_*
	igblastn \
		-query ${fasta} \
		-out ${sample_id}_igblast.fmt7 \
		-num_threads ${task.cpus} \
		-ig_seqtype Ig \
		-organism human \
		-germline_db_V "\${IGDATA}/database/bovine_V.fasta" \
		-germline_db_D "\${IGDATA}/database/bovine_D.fasta" \
		-germline_db_J "\${IGDATA}/database/bovine_J.fasta" \
		-outfmt "7 std qseq sseq btop" \
		-domain_system imgt
	"""

}

process MAKEDB {

	container 'immcantation/suite:4.5.0'
	containerOptions '--platform linux/amd64'

	tag { "${sample_id}" }
	publishDir "${params.repertoire_results}/airr", mode: 'copy'

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	tuple val(sample_id), path(igblast_out), path(fasta)
	path germlines

	output:
	tuple val(sample_id), path("${sample_id}_db-pass.tsv")

	script:
	"""
	find -L . -type f \\( -iname "*IGHV*" -o -iname "*IGKV*" -o -iname "*IGLV*" \\) | while read f; do
		cat "\$f" >> combined_V.fasta 2>/dev/null
	done
	find -L . -type f -iname "*IGHD*" | while read f; do
		cat "\$f" >> combined_D.fasta 2>/dev/null
	done
	find -L . -type f \\( -iname "*IGHJ*" -o -iname "*IGKJ*" -o -iname "*IGLJ*" \\) | while read f; do
		cat "\$f" >> combined_J.fasta 2>/dev/null
	done
	touch combined_V.fasta combined_D.fasta combined_J.fasta

	# --infer-junction extracts junction/CDR3 from alignments directly,
	# since bovine has no real IgBLAST auxiliary data (see BUILD_IGBLAST_DB)
	MakeDb.py igblast \
		-i ${igblast_out} \
		-s ${fasta} \
		-r combined_V.fasta combined_D.fasta combined_J.fasta \
		--extended \
		--partial \
		--infer-junction \
		--format airr \
		-o ${sample_id}_db.tsv

	if [[ -f "${sample_id}_db.tsv" ]] && [[ ! -f "${sample_id}_db-pass.tsv" ]]; then
		mv ${sample_id}_db.tsv ${sample_id}_db-pass.tsv
	elif [[ -f "${sample_id}_db_db-pass.tsv" ]]; then
		mv ${sample_id}_db_db-pass.tsv ${sample_id}_db-pass.tsv
	fi

	if [[ ! -f "${sample_id}_db-pass.tsv" ]]; then
		echo -e "sequence_id\\tv_call\\td_call\\tj_call\\tsequence\\tproductive" > ${sample_id}_db-pass.tsv
	fi
	"""

}

process INFER_JUNCTION {

	container 'immcantation/suite:4.5.0'
	containerOptions '--platform linux/amd64'

	tag { "${sample_id}" }
	publishDir "${params.repertoire_results}/airr_junction_filled", mode: 'copy'

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	tuple val(sample_id), path(airr_tsv)
	path germlines
	path infer_script

	output:
	tuple val(sample_id), path("${sample_id}_db-pass.tsv")

	script:
	"""
	# Heavy-chain J germline only -- the anchor fallback targets CDR3H3.
	# IgBLAST's junction heuristic silently leaves junction_aa blank for
	# some bovine heavy chain sequences (e.g. ultra-long CDR3H3, which
	# exceed its internal detection window) even though V/D/J were called.
	find -L . -type f -iname "*IGHJ*" | while read f; do
		cat "\$f" >> ighj_germline.fasta 2>/dev/null
	done
	touch ighj_germline.fasta

	python3 ${infer_script} \
		--input ${airr_tsv} \
		--germline ighj_germline.fasta \
		--output ${sample_id}_db-pass.tsv
	"""

}

process FILTER_PRODUCTIVE {

	container 'immcantation/suite:4.5.0'
	containerOptions '--platform linux/amd64'

	tag { "${sample_id}" }
	publishDir "${params.repertoire_results}/filtered", mode: 'copy'

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	tuple val(sample_id), path(airr_tsv)

	output:
	tuple val(sample_id), path("${sample_id}_productive.tsv")

	script:
	if (params.skip_productive_filter)
		"""
		# Productivity filtering skipped (--skip_productive_filter true,
		# the default -- our germlines aren't IMGT-gapped, so productivity
		# calls aren't reliable)
		cp ${airr_tsv} ${sample_id}_productive.tsv
		"""
	else
		"""
		line_count=\$(wc -l < ${airr_tsv})

		if [[ \$line_count -le 1 ]]; then
			cp ${airr_tsv} ${sample_id}_productive.tsv
		else
			ParseDb.py select \
				-d ${airr_tsv} \
				-f productive \
				-u T TRUE True true \
				-o ${sample_id}_productive.tsv || cp ${airr_tsv} ${sample_id}_productive.tsv
		fi

		if [[ ! -f "${sample_id}_productive.tsv" ]]; then
			cp ${airr_tsv} ${sample_id}_productive.tsv
		fi
		"""

}

process SUMMARIZE_VDJ {

	container 'immcantation/suite:4.5.0'
	containerOptions '--platform linux/amd64'

	publishDir "${params.repertoire_results}/reports", mode: 'copy'

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	path airr_tsvs
	path summarize_script

	output:
	path "vdj_summary.tsv"

	script:
	"""
	python3 ${summarize_script} \
		--inputs ${airr_tsvs} \
		--output vdj_summary.tsv
	"""

}

process EXTRACT_ULTRALONG_CLONES {

	// Plain stdlib Python -- no need for the Immcantation container/emulation
	publishDir "${params.repertoire_results}/ultralong_cdrh3", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	path vdj_summary
	path majority_fastas
	path extract_script

	output:
	path "ultralong_candidates.tsv"
	path "ultralong_candidates.fasta"

	script:
	"""
	python3 ${extract_script} \
		--vdj_summary ${vdj_summary} \
		--fasta_dir . \
		--min_cdr3_aa ${params.min_ultralong_cdr3_aa} \
		--v_gene ${params.ultralong_v_gene} \
		--d_gene ${params.ultralong_d_gene} \
		--j_gene ${params.ultralong_j_gene} \
		--output_tsv ultralong_candidates.tsv \
		--output_fasta ultralong_candidates.fasta
	"""

}

process DEFINE_CLONES {

	container 'immcantation/suite:4.5.0'
	containerOptions '--platform linux/amd64'

	tag { "${sample_id}" }
	publishDir "${params.repertoire_results}/clones", mode: 'copy'

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	tuple val(sample_id), path(airr_tsv)

	output:
	path "${sample_id}_clones.tsv"

	script:
	"""
	line_count=\$(wc -l < ${airr_tsv})
	has_junction=\$(head -1 ${airr_tsv} | grep -c "junction" || echo "0")

	if [[ \$line_count -le 1 ]] || [[ \$has_junction -eq 0 ]]; then
		cp ${airr_tsv} ${sample_id}_clones.tsv
	else
		DefineClones.py -d ${airr_tsv} \
			--act set \
			--model ham \
			--norm len \
			--dist ${params.clone_threshold} \
			--format airr \
			-o ${sample_id}_clones.tsv || cp ${airr_tsv} ${sample_id}_clones.tsv
	fi

	if [[ ! -f "${sample_id}_clones.tsv" ]]; then
		cp ${airr_tsv} ${sample_id}_clones.tsv
	fi
	"""

}

process DIVERSITY_ANALYSIS {

	container 'immcantation/suite:4.5.0'
	containerOptions '--platform linux/amd64'

	publishDir "${params.repertoire_results}/diversity", mode: 'copy'

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	path clone_files

	output:
	path "plots/*", emit: plots
	path "stats/*", emit: stats

	script:
	"""
	mkdir -p plots stats

	Rscript - <<'RSCRIPT'

	library(alakazam)
	library(shazam)
	library(ggplot2)
	library(dplyr)
	library(airr)

	files <- list.files(pattern = "_clones.tsv\$", full.names = TRUE)

	if (length(files) == 0) {
		stop("No clone files found")
	}

	db_list <- lapply(files, function(f) {
		df <- read_rearrangement(f)
		sample_name <- gsub("_clones.tsv", "", basename(f))
		sample_name <- gsub("_[0-9]+\$", "", sample_name)
		sample_name <- gsub("_nogroup\$", "", sample_name)
		df\$sample_id <- sample_name
		return(df)
	})
	db <- bind_rows(db_list)

	# Auto-detect the dominant chain type per barcode and filter out the
	# minority chain as contamination
	db <- db %>%
		mutate(
			barcode = gsub("_(heavy|light)\$", "", sample_id),
			chain_type = ifelse(grepl("_heavy\$", sample_id), "heavy",
						 ifelse(grepl("_light\$", sample_id), "light", "unknown"))
		)

	chain_counts <- db %>%
		group_by(barcode, chain_type) %>%
		summarise(n = n(), .groups = "drop") %>%
		filter(chain_type != "unknown")

	dominant_chains <- chain_counts %>%
		group_by(barcode) %>%
		slice_max(n, n = 1, with_ties = FALSE) %>%
		select(barcode, dominant_chain = chain_type)

	valid_samples <- dominant_chains %>%
		mutate(sample_id = paste0(barcode, "_", dominant_chain)) %>%
		pull(sample_id)

	contamination <- setdiff(unique(db\$sample_id), valid_samples)
	if (length(contamination) > 0) {
		message("Auto-detected contamination (minority chain types): ", paste(contamination, collapse = ", "))
		db <- db %>% filter(sample_id %in% valid_samples)
	}
	message("Remaining samples after filtering: ", paste(unique(db\$sample_id), collapse = ", "))

	db <- db %>% select(-barcode, -chain_type)

	has_clone_id <- "clone_id" %in% colnames(db)
	has_junction_aa <- "junction_aa" %in% colnames(db)

	if (!has_clone_id) {
		message("Warning: clone_id column not found in data. Clone-based analyses will be skipped.")
	}

	stats <- db %>%
		group_by(sample_id) %>%
		summarise(
			total_sequences = n(),
			unique_clones = if (has_clone_id) n_distinct(clone_id, na.rm = TRUE) else NA_integer_,
			productive = sum(productive == TRUE | productive == "T", na.rm = TRUE),
			mean_cdr3_length = if (has_junction_aa) mean(nchar(as.character(junction_aa)), na.rm = TRUE) else NA_real_,
			median_cdr3_length = if (has_junction_aa) median(nchar(as.character(junction_aa)), na.rm = TRUE) else NA_real_
		)
	write.csv(stats, "stats/basic_stats.csv", row.names = FALSE)

	if (has_junction_aa) {
		db\$cdr3_length <- nchar(as.character(db\$junction_aa))

		p1 <- ggplot(db, aes(x = cdr3_length, fill = sample_id)) +
			geom_histogram(binwidth = 1, position = "dodge", alpha = 0.7) +
			labs(title = "CDR3 Length Distribution",
				 x = "CDR3 Length (amino acids)",
				 y = "Count") +
			theme_minimal() +
			theme(legend.position = "bottom")
		ggsave("plots/cdr3_length_distribution.pdf", p1, width = 10, height = 6)
		ggsave("plots/cdr3_length_distribution.png", p1, width = 10, height = 6, dpi = 150)
	} else {
		message("Skipping CDR3 length distribution due to missing junction_aa column")
	}

	v_usage <- db %>%
		filter(!is.na(v_call)) %>%
		mutate(v_gene = gsub("\\\\*.*", "", v_call)) %>%
		group_by(sample_id, v_gene) %>%
		summarise(count = n(), .groups = "drop") %>%
		group_by(sample_id) %>%
		mutate(freq = count / sum(count))

	write.csv(v_usage, "stats/v_gene_usage.csv", row.names = FALSE)

	p2 <- ggplot(v_usage, aes(x = reorder(v_gene, -freq), y = freq, fill = sample_id)) +
		geom_bar(stat = "identity", position = "dodge") +
		labs(title = "V Gene Usage",
			 x = "V Gene",
			 y = "Frequency") +
		theme_minimal() +
		theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
			  legend.position = "bottom")
	ggsave("plots/v_gene_usage.pdf", p2, width = 14, height = 6)
	ggsave("plots/v_gene_usage.png", p2, width = 14, height = 6, dpi = 150)

	d_usage <- db %>%
		filter(grepl("_heavy\$", sample_id)) %>%
		filter(!is.na(d_call) & d_call != "") %>%
		mutate(d_gene = gsub("\\\\*.*", "", d_call)) %>%
		group_by(sample_id, d_gene) %>%
		summarise(count = n(), .groups = "drop") %>%
		group_by(sample_id) %>%
		mutate(freq = count / sum(count))

	if (nrow(d_usage) > 0) {
		write.csv(d_usage, "stats/d_gene_usage.csv", row.names = FALSE)

		p_d <- ggplot(d_usage, aes(x = reorder(d_gene, -freq), y = freq, fill = sample_id)) +
			geom_bar(stat = "identity", position = "dodge") +
			labs(title = "D Gene Usage (Heavy Chain)",
				 x = "D Gene",
				 y = "Frequency") +
			theme_minimal() +
			theme(axis.text.x = element_text(angle = 45, hjust = 1),
				  legend.position = "bottom")
		ggsave("plots/d_gene_usage.pdf", p_d, width = 10, height = 6)
		ggsave("plots/d_gene_usage.png", p_d, width = 10, height = 6, dpi = 150)
	} else {
		message("No D gene calls found (D genes only present in heavy chain)")
	}

	j_usage <- db %>%
		filter(!is.na(j_call)) %>%
		mutate(j_gene = gsub("\\\\*.*", "", j_call)) %>%
		group_by(sample_id, j_gene) %>%
		summarise(count = n(), .groups = "drop") %>%
		group_by(sample_id) %>%
		mutate(freq = count / sum(count))

	write.csv(j_usage, "stats/j_gene_usage.csv", row.names = FALSE)

	p3 <- ggplot(j_usage, aes(x = reorder(j_gene, -freq), y = freq, fill = sample_id)) +
		geom_bar(stat = "identity", position = "dodge") +
		labs(title = "J Gene Usage",
			 x = "J Gene",
			 y = "Frequency") +
		theme_minimal() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1),
			  legend.position = "bottom")
	ggsave("plots/j_gene_usage.pdf", p3, width = 10, height = 6)
	ggsave("plots/j_gene_usage.png", p3, width = 10, height = 6, dpi = 150)

	if (has_clone_id) {
		clone_sizes <- db %>%
			filter(!is.na(clone_id)) %>%
			group_by(sample_id, clone_id) %>%
			summarise(clone_size = n(), .groups = "drop")

		write.csv(clone_sizes, "stats/clone_sizes.csv", row.names = FALSE)

		if (nrow(clone_sizes) > 0) {
			p4 <- ggplot(clone_sizes, aes(x = clone_size, fill = sample_id)) +
				geom_histogram(binwidth = 1, position = "dodge", alpha = 0.7) +
				scale_x_log10() +
				labs(title = "Clone Size Distribution",
					 x = "Clone Size (log10)",
					 y = "Count") +
				theme_minimal() +
				theme(legend.position = "bottom")
			ggsave("plots/clone_size_distribution.pdf", p4, width = 10, height = 6)
			ggsave("plots/clone_size_distribution.png", p4, width = 10, height = 6, dpi = 150)
		}

		if (nrow(db) > 0 && any(!is.na(db\$clone_id))) {

			div_curve <- tryCatch({
				alphaDiversity(db, group = "sample_id", clone = "clone_id",
							  min_q = 0, max_q = 4, step_q = 0.1,
							  ci = 0.95, nboot = 100)
			}, error = function(e) {
				message("Could not calculate diversity curve: ", e\$message)
				NULL
			})

			if (!is.null(div_curve)) {
				p5 <- plot(div_curve, legend_title = "Sample") +
					labs(title = "Repertoire Diversity (Hill Numbers)") +
					theme_minimal()
				ggsave("plots/diversity_curve.pdf", p5, width = 10, height = 6)
				ggsave("plots/diversity_curve.png", p5, width = 10, height = 6, dpi = 150)

				write.csv(div_curve@diversity, "stats/diversity_values.csv", row.names = FALSE)
			}

			rarefaction <- tryCatch({
				estimateAbundance(db, group = "sample_id", clone = "clone_id",
								ci = 0.95, nboot = 100)
			}, error = function(e) {
				message("Could not calculate rarefaction: ", e\$message)
				NULL
			})

			if (!is.null(rarefaction)) {
				tryCatch({
					p6 <- plot(rarefaction, legend_title = "Sample") +
						labs(title = "Clonal Abundance Rarefaction") +
						theme_minimal()
					ggsave("plots/rarefaction_curve.pdf", p6, width = 10, height = 6)
					ggsave("plots/rarefaction_curve.png", p6, width = 10, height = 6, dpi = 150)
				}, error = function(e) {
					message("Could not plot rarefaction curve: ", e\$message)
				})
			}
		}

		diversity_summary <- db %>%
			filter(!is.na(clone_id)) %>%
			group_by(sample_id) %>%
			summarise(
				total_sequences = n(),
				unique_clones = n_distinct(clone_id),
				simpson_index = 1 - sum((table(clone_id)/n())^2),
				shannon_index = -sum((table(clone_id)/n()) * log(table(clone_id)/n())),
				chao1 = n_distinct(clone_id) + (sum(table(clone_id) == 1)^2) / (2 * max(1, sum(table(clone_id) == 2))),
				.groups = "drop"
			)

		write.csv(diversity_summary, "stats/diversity_summary.csv", row.names = FALSE)
	} else {
		message("Skipping clone-based analyses due to missing clone_id column")
		write.csv(data.frame(), "stats/clone_sizes.csv", row.names = FALSE)
		write.csv(data.frame(), "stats/diversity_summary.csv", row.names = FALSE)
	}

	message("Diversity analysis complete!")

	RSCRIPT
	"""

}

process GENERATE_REPORT {

	container 'immcantation/suite:4.5.0'
	containerOptions '--platform linux/amd64'

	publishDir "${params.repertoire_results}", mode: 'copy'

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	path plots
	path stats

	output:
	path "repertoire_report.html"

	script:
	"""
	cat > repertoire_report.html << 'HTML'
<!DOCTYPE html>
<html>
<head>
	<title>Bovine IgG Repertoire Analysis Report</title>
	<style>
		body { font-family: Arial, sans-serif; margin: 40px; }
		h1 { color: #2c3e50; }
		h2 { color: #34495e; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
		.plot-container { margin: 20px 0; text-align: center; }
		.plot-container img { max-width: 100%; border: 1px solid #ddd; }
		.section { margin: 30px 0; }
	</style>
</head>
<body>
	<h1>Bovine IgG Repertoire Analysis Report</h1>

	<div class="section">
		<h2>Summary Statistics</h2>
		<p>See stats/ directory for detailed CSV files.</p>
	</div>

	<div class="section">
		<h2>CDR3 Length Distribution</h2>
		<div class="plot-container">
			<img src="diversity/plots/cdr3_length_distribution.png" alt="CDR3 Length Distribution">
		</div>
	</div>

	<div class="section">
		<h2>V Gene Usage</h2>
		<div class="plot-container">
			<img src="diversity/plots/v_gene_usage.png" alt="V Gene Usage">
		</div>
	</div>

	<div class="section">
		<h2>J Gene Usage</h2>
		<div class="plot-container">
			<img src="diversity/plots/j_gene_usage.png" alt="J Gene Usage">
		</div>
	</div>

	<div class="section">
		<h2>Clone Size Distribution</h2>
		<div class="plot-container">
			<img src="diversity/plots/clone_size_distribution.png" alt="Clone Size Distribution">
		</div>
	</div>

	<div class="section">
		<h2>Diversity Analysis</h2>
		<div class="plot-container">
			<img src="diversity/plots/diversity_curve.png" alt="Diversity Curve">
		</div>
		<div class="plot-container">
			<img src="diversity/plots/rarefaction_curve.png" alt="Rarefaction Curve">
		</div>
	</div>

	<div class="section">
		<h2>Output Files</h2>
		<ul>
			<li><strong>igblast/</strong> - Raw IgBLAST output</li>
			<li><strong>airr/</strong>, <strong>airr_junction_filled/</strong> - AIRR-formatted sequence annotations</li>
			<li><strong>filtered/</strong> - Productivity-filtered sequences</li>
			<li><strong>clones/</strong> - Clone assignments</li>
			<li><strong>reports/vdj_summary.tsv</strong> - Combined V/D/J call summary</li>
			<li><strong>diversity/stats/</strong> - CSV files with diversity metrics</li>
			<li><strong>diversity/plots/</strong> - PDF and PNG plots</li>
		</ul>
	</div>

</body>
</html>
HTML
	"""

}

// --------------------------------------------------------------- //
