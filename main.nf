/*
    Define pipeline parameters
    These specify required inputs and configurable analysis settings
*/

params {
    // Required inputs
    reference: Path                            // FASTA containing the ddRAD target loci.
    samplesheet: String                        // Path to the sample sheet used to locate input reads.
    slurm_account: String                      // SLURM account for submitted jobs.

    // Parallel parameters
    nchunks: Integer = 200                     // Number of groups of whole loci called in parallel.
    lower_percentile: Integer = 1              // Exclude loci below this percentile.
    upper_percentile: Integer = 99             // Exclude loci above this percentile.

    // FreeBayes: genotype model and population priors
    ploidy: Integer = 2                        // Number of chromosome copies per sample.
    theta: Double = 0.001d                     // Expected population diversity used by the population prior.
    pooled_discrete: Boolean = false           // Model samples as discrete pools (poolseq style) rather than individuals.
    hwe_priors_off: Boolean = true             // Disable the Hardy-Weinberg equilibrium prior.

    // FreeBayes: read and base filters
    min_mapping_quality: Integer = 20          // Ignore reads with mapping quality below 20.
    min_base_quality: Integer = 20             // Ignore allele observations with base quality below 20.
    min_supporting_allele_qsum: Integer = 0    // Minimum summed quality supporting an allele; 0 adds no threshold.

    // FreeBayes: site and candidate-allele thresholds
    min_coverage: Integer = 10                 // Minimum site coverage required for FreeBayes to process a site.
    limit_coverage: Integer = 100000           // Coverage-limit setting; verify its behaviour in your installed FreeBayes.
    min_alternate_fraction: Double = 0.05d     // Minimum alternate-allele fraction within a sample to evaluate a site.
    min_alternate_count: Integer = 2           // Minimum alternate-supporting observations within a sample.
    use_best_n_alleles: Integer = 4            // Evaluate at most the four best-supported SNP alleles.
  
    // FreeBayes: haplotype observations
    no_partial_observations: Boolean = false   // If true, exclude reads that do not span the detection window.

}

/*
    Define Nextflow processes
    Each process wraps a tool or a discrete processing step
*/

// Nextflow process that indexes reference fasta for alignment
process INDEX_REFERENCE {
    tag "${reference}"
    conda 'bioconda::minibwa=0.7 bioconda::samtools=1.24'
    cpus 2
    memory { 4.GB * task.attempt }
    time '1h'

    input:
    path reference, stageAs: 'reference.fa'

    output:
    tuple path("reference.fa"), path("*.{fai,sa,l2b,mbw,dict}"), emit: fasta_indexed

    script:
    """
    set -euo pipefail

    minibwa index -t ${task.cpus} reference.fa
    samtools faidx reference.fa
    """
}

// Nextflow process that trims adapters and aligns each library to reference fasta
process MINIBWA_ALIGN {
    tag "${library_id}"
    conda 'bioconda::minibwa=0.7 bioconda::samtools=1.24 bioconda::fastp=1.3.7'
    cpus 8
    memory { 16.GB * task.attempt }
    time { 8.h * task.attempt }

    input:
    tuple val(sample_name), val(library_id), path(read1), path(read2)
    tuple path(ref_fasta), path(index_files)

    output:
    tuple val(sample_name),
          path("${library_id}.bam"),
          path("${library_id}.bam.bai"),
          emit: bam

    tuple val(sample_name),
          path("${library_id}.fastp.json"),
          path("${library_id}.fastp.html"),
          emit: fastp_qc

    script:
    // All three commands run concurrently in the pipe.
    def trim_threads     = Math.max(1, task.cpus.intdiv(4))
    def sort_threads     = task.cpus >= 7 ? 1 : 0
    def overhead_threads = task.cpus >= 4 ? 3 : 1
    def aln_threads      = Math.max(1, task.cpus - trim_threads - sort_threads - overhead_threads )

    // Read group, used for library merging downstream
    def read_group = "@RG\\tID:${library_id}\\tSM:${sample_name}\\tPL:ILLUMINA"

    """
    set -euo pipefail

    fastp \\
        --in1 "${read1}" \\
        --in2 "${read2}" \\
        --stdout \\
        --thread ${trim_threads} \\
        --json "${library_id}.fastp.json" \\
        --html "${library_id}.fastp.html" \\
        --report_title "${library_id}" |
    minibwa map \\
        -t ${aln_threads} \\
        -R '${read_group}' \\
        "${ref_fasta}" - |
    samtools sort \\
        --threads ${sort_threads} \\
        -m 2G \\
        -O BAM \\
        -o "${library_id}.bam"

    samtools index "${library_id}.bam"
    """
}

// Nextflow process that merges per-library bams, if there are multiple for each sample
process MERGE_SAMPLE_BAMS {
    tag "${sample_name}"
    conda 'bioconda::samtools=1.24'
    cpus 4
    memory { 8.GB * task.attempt }
    time { 4.h * task.attempt }

    input:
    tuple val(sample_name), path(bams)

    output:
    tuple val(sample_name),
          path("${sample_name}.bam"),
          path("${sample_name}.bam.bai"),
          emit: bam

    script:
    """
    set -euo pipefail

    # Only merge if there are multiple 
    if [[ ${bams.size()} -eq 1 ]]; then
        cp "${bams[0]}" "${sample_name}.bam"
    else
        samtools merge \\
            -@ ${task.cpus} \\
            -o "${sample_name}.bam" \\
            ${bams.join(' ')}
    fi

    samtools index "${sample_name}.bam"
    """
}

// Nextflow process that calculates number of reads mapping to each contig
// then divides into nchunks with approx even workload
process MAKE_REFERENCE_CHUNKS {
    conda 'bioconda::samtools=1.24'
    tag 'cohort'
    cpus 1
    memory 2.GB
    time '1h'

    input:
    tuple path(ref_genome), path(index_files)
    path bam_files
    val nchunks
    val lower_percentile
    val upper_percentile

    output:
    path 'chunk_*.bed',     emit: chunk_beds
    path 'read_counts.tsv', emit: read_counts

    script:
    """
    set -euo pipefail
    export LC_ALL=C

    fai="${ref_genome}.fai"
    total_loci=\$(wc -l < "\$fai")

    if (( ${lower_percentile} < 0 || ${upper_percentile} > 100 || ${lower_percentile} >= ${upper_percentile} )); then
        echo "Require 0 <= lower_percentile < upper_percentile <= 100" >&2
        exit 1
    fi

    # Calculate read counts by locus from bam index: locus, mapped read segments.
    for bam in *.bam; do
        samtools idxstats "\$bam"
    done |
        awk '\$1 != "*" { reads[\$1] += \$3 }
             END { for (locus in reads) print locus, reads[locus] }' \\
        > read_counts.tsv

    if (( \$(wc -l < read_counts.tsv) != total_loci )); then
        echo "Locus counts do not match the reference index" >&2
        exit 1
    fi

    # Keep loci between the independently specified percentile ranks.
    lower_rank=\$(( total_loci * ${lower_percentile} / 100 ))
    upper_rank=\$(( (total_loci * ${upper_percentile} + 99) / 100 ))

    sort -k2,2n -k1,1 read_counts.tsv |
        awk -v lower="\$lower_rank" -v upper="\$upper_rank" \\
            'NR > lower && NR <= upper { print \$1 }' \\
        > retained.ids

    # Restore reference order for contiguous chunks and VCF concatenation.
    awk 'NR == FNR { keep[\$1] = 1; next }
         \$1 in keep { print }' \\
        retained.ids "\$fai" > retained.fai

    ncontigs=\$(wc -l < retained.fai)

    if (( ${nchunks} < 1 || ${nchunks} > ncontigs )); then
        echo "nchunks must be between 1 and \$ncontigs retained loci" >&2
        exit 1
    fi

    total_weight=\$(
        awk 'NR == FNR { reads[\$1] = \$2; next }
             { total += reads[\$1] + 1 }
             END { printf "%.0f", total }' \\
            read_counts.tsv retained.fai
    )

    # Group retained whole loci into contiguous, read-balanced chunks.
    awk -v n=${nchunks} \\
        -v ncontigs="\$ncontigs" \\
        -v total="\$total_weight" '
        BEGIN {
            OFS = "\\t"
            chunk = 0
            target = total / n
        }

        NR == FNR {
            reads[\$1] = \$2
            next
        }

        {
            weight = reads[\$1] + 1

            if (chunk < n - 1 && count > 0 &&
                (chunk_weight >= target ||
                 ncontigs - FNR + 1 == n - chunk - 1)) {
                assigned += chunk_weight
                chunk++
                chunk_weight = 0
                count = 0
                target = (total - assigned) / (n - chunk)
            }

            file = sprintf("chunk_%05d.bed", chunk)
            print \$1, 0, \$2 >> file
            chunk_weight += weight
            count++
        }
    ' read_counts.tsv retained.fai
    """
}


// Nextflow process that runs freebayes on a subset of intervals
process FREEBAYES {
    tag "${chunk_id}"
    conda 'bioconda::freebayes=1.3.10 bioconda::bcftools=1.24'
    cpus 2
    memory { 24.GB * task.attempt }
    time { 6.h * task.attempt }

    input:
    tuple val(chunk_id), path(chunk_bed)
    tuple path(ref_genome), path(genome_index_files)
    path bam_files

    output:
    tuple val(chunk_id),
          path("${chunk_id}.vcf.gz"),
          path("${chunk_id}.vcf.gz.csi"),
          emit: chunk_vcf

    script:
    """
    set -euo pipefail

    printf '%s\\n' *.bam > bam.list

    freebayes \\
        -f ${ref_genome} \\
        -L bam.list \\
        -t ${chunk_bed} \\
        --ploidy ${params.ploidy} \\
        --theta ${params.theta} \\
        --min-alternate-fraction ${params.min_alternate_fraction} \\
        --min-alternate-count ${params.min_alternate_count} \\
        --use-best-n-alleles ${params.use_best_n_alleles} \\
        --limit-coverage ${params.limit_coverage} \\
        --min-mapping-quality ${params.min_mapping_quality} \\
        --min-base-quality ${params.min_base_quality} \\
        --min-supporting-allele-qsum ${params.min_supporting_allele_qsum} \\
        --min-coverage ${params.min_coverage} \\
        ${params.pooled_discrete ? '--pooled-discrete' : ''} \\
        ${params.hwe_priors_off ? '--hwe-priors-off' : ''} \\
        ${params.no_partial_observations ? '--no-partial-observations' : ''} \\
        | bcftools view -Oz -o "${chunk_id}.vcf.gz"

    bcftools index -c "${chunk_id}.vcf.gz"
    """
}

// Nextflow process that concatenates seperate freebayes chunks into a final vcf
process CONCAT_VCFS {
    tag 'cohort'   
    conda 'bioconda::bcftools=1.24'
    cpus 2
    memory { 16.GB * task.attempt }
    time { 2.h * task.attempt }

    input:
    tuple path(vcfs), path(csis)

    output:
    tuple path('cohort.vcf.gz'), path('cohort.vcf.gz.csi'), emit: vcf

    script:
    def vcf_list = vcfs
        .toSorted { vcf -> vcf.name }
        .collect { vcf -> vcf.name }
        .join('\n')

    """
    set -euo pipefail

    printf '%s\\n' "${vcf_list}" > vcf.list

    bcftools concat \
        --threads ${task.cpus} \
        --naive \
        -f vcf.list \
        -o cohort.vcf.gz

    bcftools index -c cohort.vcf.gz
    """
}

// Nextflow process that calculates QC metrics from each merged BAM
process BAM_STATS_RIKER {
    tag "${sample}"
    conda 'bioconda::riker=0.4.1'
    cpus 4
    memory { 8.GB * task.attempt }
    time '2h'

    input:
    tuple val(sample), path(bam), path(bam_index)
    tuple path(ref_genome), path(genome_index_files)

    output:
    tuple val(sample), path('riker/*'), emit: stats

    script:
    """
    set -euo pipefail

    mkdir riker

    riker multi \\
        --threads ${task.cpus} \\
        -i "${bam}" \\
        -r "${ref_genome}" \\
        -o "riker/${sample}" \\
        --tools alignment isize basic gcbias wgs \\
        --aln::min-mapq ${params.min_mapping_quality} \\
        --aln::max-insert-size 10000 \\
        --wgs::min-mapq ${params.min_mapping_quality} \\
        --wgs::min-bq ${params.min_base_quality} \\
        --wgs::coverage-cap 10000 \\
        --isize::min-frac 0.05 \\
        --isize::deviations 10

    """
}

// Nextflow process that calculates statistics for the final cohort VCF
process VCF_STATS {
    conda 'bioconda::bcftools=1.24'
    tag 'cohort'
    cpus 4
    memory { 8.GB * task.attempt }
    time '2h'

    input:
    tuple path(vcf), path(csi)
    tuple path(ref_genome), path(genome_index_files)

    output:
    path 'vcfstats.txt', emit: vcfstats

    script:
    """
    set -euo pipefail

    bcftools stats \\
        -F "${ref_genome}" \\
        -s - \\
        "${vcf}" > vcfstats.txt
    """
}

// Nextflow process that creates multiqc report
process MULTIQC {
    conda 'bioconda::multiqc=1.35'
    cpus 1
    memory { 4.GB * task.attempt }
    time '1h'

    input:
    path multiqc_files
    path multiqc_config

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , emit: plots

    script:
    """
    set -euo pipefail

    # Prevent any loaded HPC Python modules and user packages from contaminating the Nextflow Conda environment.
    unset PYTHONPATH
    unset PYTHONHOME
    export PYTHONNOUSERSITE=1

    multiqc . \\
        --force \\
        --config "${multiqc_config}" \\
        --filename multiqc_report.html \\
        --clean-up
    """
}

/*
    Nextflow workflow 
    This is where processes are connected and pass their outputs to downstream steps
*/

workflow {
    main:
    if (!params.reference || !params.samplesheet) {
        error 'Provide both --reference and --samplesheet'
    }

    // Parse samplesheet
    def samplesheet = file(params.samplesheet, checkIfExists: true)

    read_pairs = channel.of(samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def sample_name = row.sample_name.trim()
            def read1 = file(workflow.projectDir.resolve(row.read1.trim()), checkIfExists: true)
            def read2 = file(workflow.projectDir.resolve(row.read2.trim()), checkIfExists: true)
            def library_id = "${sample_name}_${read1.name.replaceFirst(/(?i)\.(fastq|fq)(\.gz)?$/, '')}"

            tuple(sample_name, library_id, read1, read2)
        }

    // Index reference fasta
    reference_indexed = INDEX_REFERENCE(
        file(params.reference, checkIfExists: true)
    )

    // Align read pairs to indexed reference fasta
    MINIBWA_ALIGN(read_pairs, reference_indexed)

    // Group multiple bams by sample_name
    bams_by_sample = MINIBWA_ALIGN.out.bam
        .map { sample_name, bam, _bai -> tuple(sample_name, bam) }
        .groupTuple()

    // Merge bams if there are multiple fastq files per sample, then drop sample_id
    MERGE_SAMPLE_BAMS(bams_by_sample)

    // Create bam QC stats
    BAM_STATS_RIKER(
        MERGE_SAMPLE_BAMS.out.bam,
        reference_indexed
    )

    merged_bams = MERGE_SAMPLE_BAMS.out.bam
        .flatMap { _sample, bam, bai -> [bam, bai] }
        .collect()

    // Create one BED file per chunk
    MAKE_REFERENCE_CHUNKS(
        reference_indexed,
        merged_bams,
        params.nchunks,
        params.lower_percentile,
        params.upper_percentile
    )

    // One channel item per BED file; chunk_id sorts in reference order
    chunks = MAKE_REFERENCE_CHUNKS.out.chunk_beds
        .flatten()
        .map { bed -> tuple(bed.baseName, bed) }

    // Each chunk runs against the same reference and cohort BAMs
    chunk_vcfs = FREEBAYES(chunks, reference_indexed, merged_bams)

    // Concatenate freebayes called chunks into final vcf
    vcfs_to_merge = chunk_vcfs
        .collect(flat: false)
        .map { records ->
            tuple(
                records.collect { record -> record[1] },  // VCFs
                records.collect { record -> record[2] }   // CSI indexes
            )
        }

    CONCAT_VCFS(vcfs_to_merge)
    
    // Calculate statistics on final vcf
    VCF_STATS(
        CONCAT_VCFS.out.vcf,
        reference_indexed
    )

    // Combine BAM and VCF QC files for multiqc
    multiqc_files = BAM_STATS_RIKER.out.stats
        .flatMap { _sample, files -> files }
        .mix(MINIBWA_ALIGN.out.fastp_qc.flatMap { _sample, json, html -> [json, html] })
        .mix(VCF_STATS.out.vcfstats)
        .collect()

    // Use MultiQC to make a final QC report
    MULTIQC(
        multiqc_files,
        file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
    )

    publish:
    bam              = MERGE_SAMPLE_BAMS.out.bam
    read_counts      = MAKE_REFERENCE_CHUNKS.out.read_counts
    vcf              = CONCAT_VCFS.out.vcf
    riker            = BAM_STATS_RIKER.out.stats
    vcf_stats        = VCF_STATS.out.vcfstats
    multiqc_report   = MULTIQC.out.report
    multiqc_plots    = MULTIQC.out.plots
    multiqc_data     = MULTIQC.out.data
}

/*
    Publish outputs 
    This is where published outputs of the workflow are sent to your output directory
*/

output {
    bam { path 'bam' }
    vcf { path 'vcf' }
    riker { path 'qc/bam_stats' }
    read_counts { path 'qc' }
    vcf_stats { path 'qc/vcf_stats' }
    multiqc_report { path 'qc' }   
    multiqc_plots { path 'qc' }    
    multiqc_data { path 'qc' }    
}