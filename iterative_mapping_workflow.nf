#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process bwa_index {
    label 'process_bwa_index'
    tag { "iteration_${iteration}" }

    input:
    path(fasta)
    val(iteration)

    output:
    path("*.{amb,ann,bwt,pac,sa}"), emit: index

    script:
    """
    ml bwa
    bwa index ${fasta}
    """
}

process bwa_mem {
    label 'process_bwa_mem'
    tag { "iteration_${iteration}" }
    publishDir path: { "${params.outdir}/iteration_${iteration}/raw_bams/" }, mode: 'symlink'

    input:
    each path(fq)
    path(ref_index)
    val(iteration)

    output:
    path("${fq.baseName}*_rmdup.bam"), emit: sorted_rmdup_bam

    script:
    def reference = ref_index[0]
    """
    set -euo pipefail
    ml bwa samtools

    ## Mapping
    printf "STARTING mapping............................................. ${fq.baseName}\n"
    bwa mem ${params.bwa_mem_args} \\
        -R '@RG\\tID:${params.sample}\\tSM:${params.sample}\\tPL:ILLUMINA\\tLB:${params.sample}' \\
        -t ${task.cpus} ${reference} ${fq} | \\
        samtools view -@ ${task.cpus} -bh - | \\
        samtools sort -@ ${task.cpus} -o ${fq.baseName}_sorted.bam -
    samtools index ${fq.baseName}_sorted.bam
    printf "${fq.baseName} .............................................. completed\n"

    ## Filtering
    printf "STARTING filtering............................................. ${fq.baseName}\n"
    samtools view -F 4 -q ${params.mapQ} -@ ${task.cpus} -h ${fq.baseName}_sorted.bam | \\
        awk -v minlen="${params.min_readlength}" '\$1 ~ /^@/ || length(\$10) >= minlen' | \\
        samtools view -@ ${task.cpus} -bh -o ${fq.baseName}_mq${params.mapQ}_l${params.min_readlength}_filtered.bam -
    samtools index ${fq.baseName}_mq${params.mapQ}_l${params.min_readlength}_filtered.bam
    printf "${fq.baseName} .............................................. completed\n"

    ## Remove duplicates
    printf "STARTING removing duplicates............................................. ${fq.baseName}\n"
    samtools view -@ ${task.cpus} -h ${fq.baseName}_mq${params.mapQ}_l${params.min_readlength}_filtered.bam | \\
    python3 ${params.script_rmdup} | \\
    samtools view -@ ${task.cpus} -bh -o ${fq.baseName}_mq${params.mapQ}_l${params.min_readlength}_rmdup.bam -
    samtools index ${fq.baseName}_mq${params.mapQ}_l${params.min_readlength}_rmdup.bam
    printf "${fq.baseName} .............................................. completed\n"

    printf "${fq.baseName} .............................................. DONE\n\n"
    """
}

process merge_dedup_bams_all {
    label 'process_merge_dedup_bams'
    tag { "iteration_${iteration}" }

    publishDir path: { "${params.outdir}/iteration_${iteration}/merged_bams/" }, mode: 'symlink'

    input:
    path(sorted_bams)
    val(iteration)

    output:
    path("${params.sample}_merged_rmdup.bam"),     emit: merged_rmdup_bam
    path("${params.sample}_merged_rmdup.bam.bai"), emit: merged_rmdup_bai

    script:
    """
    set -euo pipefail
    ml -q samtools/1.20

    printf "Merging bam files.........................................${params.sample}\n"
    samtools merge -@ ${task.cpus} -o ${params.sample}_merged.bam ${sorted_bams}

    printf "Removing duplicates after merging.........................${params.sample}\n"
    samtools view -@ ${task.cpus} -h ${params.sample}_merged.bam | \\
    python3 ${params.script_rmdup} | \\
    samtools view -@ ${task.cpus} -bh -o ${params.sample}_merged_rmdup.bam -
    samtools index ${params.sample}_merged_rmdup.bam
    printf "${params.sample} Merging and removing duplicates................. DONE\n\n"
    """
}

process prepare_reference {
    label 'process_prepare_reference'
    tag { "iteration_${iteration}" }

    input:
    path(ref)
    val(iteration)

    output:
    path("${ref}.fai"),             emit: fai
    path("${ref.baseName}.dict"),   emit: dict

    script:
    """
    set -euo pipefail
    ml -q samtools/1.20

    samtools faidx ${ref}
    samtools dict ${ref} -o ${ref.baseName}.dict
    """
}

// ---------------------------------------------------------------------------
// GATK3 IndelRealigner — per-chromosome, phase 2 (indel iterations) only
// Without it, reads whose indel lies near the read end are aligned with
// mismatches or soft-clips instead of a gap, so the indel's alt fraction
// stays below min_alt_frac and it is never applied. SNV calls are almost
// unaffected, so in phase 1 (or with --indel_realign false) the merged BAM
// is passed through unchanged.
// ---------------------------------------------------------------------------
process gatk_indel_realigner {
    label 'process_gatk_indel_realigner'
    tag { "${chrom}_iteration_${iteration}" }

    container 'broadinstitute/gatk3:3.8-1'

    input:
    each chrom
    path(merged_bam)
    path(merged_bai)
    path(ref)
    path(fai)
    path(dict)
    val(iteration)
    val(snv_only)

    output:
    tuple val(chrom), path("${params.sample}_${chrom}.bam"), path("${params.sample}_${chrom}.bam.bai"), emit: bam

    script:
    def out = "${params.sample}_${chrom}"
    def avail_mem = task.memory ? task.memory.toGiga() - 1 : 4
    if (snv_only || !params.indel_realign)
        """
        ln -s ${merged_bam} ${out}.bam
        ln -s ${merged_bai} ${out}.bam.bai
        """
    else
        """
        java -Xmx${avail_mem}g -jar /usr/GenomeAnalysisTK.jar \\
            -T RealignerTargetCreator \\
            -nt ${task.cpus} \\
            -R ${ref} \\
            -I ${merged_bam} \\
            -L ${chrom} \\
            -o ${out}.intervals

        java -Xmx${avail_mem}g -jar /usr/GenomeAnalysisTK.jar \\
            -T IndelRealigner \\
            -R ${ref} \\
            -I ${merged_bam} \\
            -L ${chrom} \\
            -targetIntervals ${out}.intervals \\
            -o ${out}.bam

        mv ${out}.bai ${out}.bam.bai
        """
}

// ---------------------------------------------------------------------------
// bcftools mpileup/call — per-chromosome pileup calling
// BAQ is disabled (-B): it down-weights bases in clusters of mismatches,
// which in a highly diverged sample removes real alt bases. Indels are
// skipped in phase 1 (-I); in phase 2 the indel model is set by
// params.indel_model. Multiallelic sites are split and indels
// left-aligned so each record carries a single ALT allele with its own AD.
// ---------------------------------------------------------------------------
process bcftools_call {
    label 'process_bcftools_call'
    tag { "${chrom}_iteration_${iteration}" }

    input:
    tuple val(chrom), path(bam), path(bai)
    path(ref)
    path(fai)
    val(iteration)
    val(snv_only)

    output:
    tuple val(chrom), path("${params.sample}_${chrom}_raw.bcf"), emit: raw_bcf

    script:
    def prefix = "${params.sample}_${chrom}"
    def indel_opt = snv_only ? '-I' : (params.indel_model == 'default' ? '' : "--${params.indel_model}")
    """
    set -euo pipefail
    ml -q bcftools

    bcftools mpileup -B ${indel_opt} \\
        -q ${params.mapQ} -Q ${params.baseQ} \\
        -a FORMAT/AD,FORMAT/DP \\
        --ignore-RG \\
        -f ${ref} -r ${chrom} -Ou \\
        ${bam} | \\
    bcftools call -m -v -Ou | \\
    bcftools norm -m- -f ${ref} -Ob --threads ${task.cpus} -o ${prefix}_raw.bcf
    """
}

// ---------------------------------------------------------------------------
// Select sites to apply to the consensus by read counts, not genotype:
// alt reads >= min_alt_reads, alt fraction >= min_alt_frac, depth < max_depth.
// GT is then forced to 1/1, otherwise bcftools consensus writes IUPAC
// codes for 0/1 records.
// ---------------------------------------------------------------------------
process filter_bcftools_variants {
    label 'process_filter_variants'
    tag { "${chrom}_iteration_${iteration}" }

    publishDir path: { "${params.outdir}/iteration_${iteration}/variants/" }, mode: 'copy'

    input:
    tuple val(chrom), path(raw_bcf)
    val(iteration)
    val(snv_only)

    output:
    tuple val(chrom), \
        path("${params.sample}_${chrom}_alt${params.min_alt_reads}_af${params.min_alt_frac}_DPmax${params.max_depth}${snv_only ? '_no_indels' : ''}_homalt.bcf"), \
        path("${params.sample}_${chrom}_alt${params.min_alt_reads}_af${params.min_alt_frac}_DPmax${params.max_depth}${snv_only ? '_no_indels' : ''}_homalt.bcf.csi"), \
        emit: filtered_bcf

    script:
    def out = "${params.sample}_${chrom}_alt${params.min_alt_reads}_af${params.min_alt_frac}_DPmax${params.max_depth}${snv_only ? '_no_indels' : ''}_homalt.bcf"
    def var_type = snv_only ? '-v snps' : ''
    """
    set -euo pipefail
    ml -q bcftools

    bcftools view ${var_type} \\
        -i "FMT/AD[0:1]>=${params.min_alt_reads} & FMT/AD[0:1]/FMT/DP>=${params.min_alt_frac} & FMT/DP<${params.max_depth}" \\
        -Ou ${raw_bcf} | \\
    bcftools +setGT -Ob -o ${out} -- -t a -n c:1/1
    bcftools index ${out}
    """
}

// ---------------------------------------------------------------------------
// DeepVariant — per-chromosome variant calling (runs in container)
// Outputs raw VCF; filtering is done separately since bcftools is not
// available inside the DeepVariant container.
// ---------------------------------------------------------------------------
process deepvariant_call {
    label 'process_deepvariant'
    tag { "${chrom}_iteration_${iteration}" }

    container 'google/deepvariant:1.10.0'

    input:
    each chrom
    path(merged_bam)
    path(merged_bai)
    path(ref)
    path(fai)
    val(iteration)

    output:
    tuple val(chrom), path("${params.sample}_${chrom}_raw.vcf.gz"), emit: raw_vcf

    script:
    """
    /opt/deepvariant/bin/run_deepvariant \\
        --model_type=WGS \\
        --ref=${ref} \\
        --reads=${merged_bam} \\
        --regions=${chrom} \\
        --output_vcf=${params.sample}_${chrom}_raw.vcf.gz \\
        --postprocess_variants_extra_args="qual_filter=0" \\
        --num_shards=${task.cpus} \\
        --disable_small_model=true \\
        --vcf_stats_report=true \\
        --logging_dir=logs
    """
}

process filter_dv_variants {
    label 'process_filter_variants'
    tag { "${chrom}_iteration_${iteration}" }

    publishDir path: { "${params.outdir}/iteration_${iteration}/variants/" }, mode: 'copy'

    input:
    tuple val(chrom), path(raw_vcf)
    val(iteration)
    val(snv_only)

    output:
    tuple val(chrom), \
        path("${params.sample}_${chrom}_DP${params.min_depth}-${params.max_depth}${snv_only ? '_no_indels' : ''}_homalt.bcf"), \
        path("${params.sample}_${chrom}_DP${params.min_depth}-${params.max_depth}${snv_only ? '_no_indels' : ''}_homalt.bcf.csi"), \
        emit: filtered_bcf

    script:
    def prefix = "${params.sample}_${chrom}"
    def dp = "DP${params.min_depth}-${params.max_depth}"
    def indel_tag = snv_only ? '_no_indels' : ''
    def indel_filter = snv_only ? \
        "bcftools view -v snps -Ob --threads ${task.cpus} -o ${prefix}_${dp}_no_indels.bcf ${prefix}_${dp}.bcf" : \
        ""
    def homalt_input = snv_only ? "${prefix}_${dp}_no_indels.bcf" : "${prefix}_${dp}.bcf"
    """
    set -euo pipefail
    ml -q bcftools

    ## Filter: depth + quality
    bcftools filter -i "FMT/DP>=${params.min_depth} & FMT/DP<${params.max_depth} & QUAL>=${params.variant_quality}" \\
        -Ob --threads ${task.cpus} \\
        -o ${prefix}_${dp}.bcf \\
        ${raw_vcf}

    ## Conditionally remove indels (phase 1: SNV-only)
    ${indel_filter}

    ## Keep homozygous-alt only
    bcftools filter -i 'GT="1/1"' -Ob --threads ${task.cpus} \\
        -o ${prefix}_${dp}${indel_tag}_homalt.bcf \\
        ${homalt_input}
    bcftools index ${prefix}_${dp}${indel_tag}_homalt.bcf
    """
}


process merge_chrom_bcfs {
    label 'process_merge_chrom_bcfs'
    tag { "iteration_${iteration}" }

    publishDir path: { "${params.outdir}/iteration_${iteration}/variants/" }, mode: 'copy'

    input:
    path(bcfs)
    path(bcf_indexes)
    path(fai)
    val(iteration)

    output:
    path("${params.sample}_merged_homalt.bcf"),     emit: merged_bcf
    path("${params.sample}_merged_homalt.bcf.csi"), emit: merged_bcf_index

    script:
    """
    set -euo pipefail
    ml -q bcftools

    # Build BCF file list in reference chromosome order (from FAI)
    cut -f1 ${fai} | while read -r chrom; do
        ls ${params.sample}_\${chrom}_*homalt.bcf 2>/dev/null || true
    done > bcf_list.txt

    bcftools concat -a -Ob -o ${params.sample}_merged_homalt.bcf --file-list bcf_list.txt
    bcftools index ${params.sample}_merged_homalt.bcf
    """
}

process make_consensus {
    label 'process_make_consensus'
    tag { "iteration_${iteration}" }

    publishDir path: { "${params.outdir}/iteration_${iteration}/reference/" }, mode: 'copy'

    input:
    tuple path(filtered_bcf), path(filtered_bcf_index)
    path(ref)
    val(iteration)

    output:
    path("${params.sample}_consensus_itr${iteration}.fasta"),   emit: consensus_fasta
    path("${params.sample}_consensus_itr${iteration}.log"),     emit: consensus_log

    script:
    """
    set -euo pipefail
    ml -q bcftools
    bcftools consensus -f ${ref} -o ${params.sample}_consensus_itr${iteration}.fasta ${filtered_bcf} 2> ${params.sample}_consensus_itr${iteration}.log
    """
}

workflow iterative_mapping {
    take:
    ch_ref
    fastq_files
    ch_iteration
    ch_snv_only       // val channel: true = SNV-only (phase 1), false = SNPs + indels (phase 2)

    main:
    // index the genome
    bwa_index(ch_ref, ch_iteration)
    ch_ref_index = ch_ref.combine(bwa_index.out.index)

    // map all fastq files in parallel with bwa mem
    bwa_mem(fastq_files, ch_ref_index, ch_iteration)
    ch_filtered_bams = bwa_mem.out.sorted_rmdup_bam.collect()

    // merge and remove duplicates for all libraries
    merge_dedup_bams_all(ch_filtered_bams, ch_iteration)

    // prepare reference (fai + dict) — dict is needed for GATK3 indel realignment
    prepare_reference(ch_ref, ch_iteration)

    // convert to value channels so they can be reused across all per-chrom process invocations
    ch_fai  = prepare_reference.out.fai.first()
    ch_dict = prepare_reference.out.dict.first()

    // chromosomes to call, one per line; blank lines and '#' comments are skipped
    def chrom_list = file(params.chromosomes).readLines()
        .collect { line -> line.trim() }
        .findAll { line -> line && !line.startsWith('#') }

    // per-chromosome variant calling — caller selected via params.caller
    if (params.caller == 'bcftools') {
        gatk_indel_realigner(
            chrom_list,
            merge_dedup_bams_all.out.merged_rmdup_bam,
            merge_dedup_bams_all.out.merged_rmdup_bai,
            ch_ref,
            ch_fai,
            ch_dict,
            ch_iteration,
            ch_snv_only
        )
        bcftools_call(
            gatk_indel_realigner.out.bam,
            ch_ref,
            ch_fai,
            ch_iteration,
            ch_snv_only
        )
        filter_bcftools_variants(
            bcftools_call.out.raw_bcf,
            ch_iteration,
            ch_snv_only
        )
        ch_chrom_bcfs = filter_bcftools_variants.out.filtered_bcf
    } else if (params.caller == 'dv') {
        deepvariant_call(
            chrom_list,
            merge_dedup_bams_all.out.merged_rmdup_bam,
            merge_dedup_bams_all.out.merged_rmdup_bai,
            ch_ref,
            prepare_reference.out.fai,
            ch_iteration
        )
        filter_dv_variants(
            deepvariant_call.out.raw_vcf,
            ch_iteration,
            ch_snv_only
        )
        ch_chrom_bcfs = filter_dv_variants.out.filtered_bcf
    } else {
        error("Unknown caller '${params.caller}' -- use 'bcftools' or 'dv'")
    }

    // gather: split the (chrom, bcf, bcf_index) tuples so BCFs and indexes can be collected separately
    ch_chrom_bcfs.multiMap { _chrom, bcf, bcf_idx ->
        bcfs: bcf
        idxs: bcf_idx
    }.set { ch_chrom_split }

    // merge per-chromosome BCFs into one, in reference chromosome order
    merge_chrom_bcfs(
        ch_chrom_split.bcfs.collect(),
        ch_chrom_split.idxs.collect(),
        ch_fai,
        ch_iteration
    )

    // make consensus sequence from merged BCF
    ch_filtered_bcf = merge_chrom_bcfs.out.merged_bcf
        .combine(merge_chrom_bcfs.out.merged_bcf_index)
    make_consensus(ch_filtered_bcf, ch_ref, ch_iteration)

    emit:
    make_consensus.out.consensus_fasta
}
