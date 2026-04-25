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

    // Define the script to execute
    script:
    """
    ml bwa
    bwa index ${fasta}
    """
}

process bwa_aln_samse {
    label 'process_bwa_aln_samse'
    tag { "iteration_${iteration}" }
    publishDir "${params.outdir}/iteration_${iteration}/raw_bams/", mode: 'symlink'

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
    bwa aln -l 16500 -n 0.01 -o 2 -t ${task.cpus} ${reference} ${fq} | \\
        bwa samse ${reference} - ${fq} | \\
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

    publishDir "${params.outdir}/iteration_${iteration}/merged_bams/", mode: 'symlink'

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
    samtools view -@ ${task.cpus} -bh -o ${params.sample}_merged_rmdup_norg.bam -
    printf "${params.sample} Merging and removing duplicates................. DONE\n\n"

    printf "Adding read group header..................................${params.sample}\n"
    samtools addreplacerg \\
        -r "ID:${params.sample}\tSM:${params.sample}\tPL:ILLUMINA\tLB:${params.sample}" \\
        -@ ${task.cpus} \\
        -o ${params.sample}_merged_rmdup.bam \\
        ${params.sample}_merged_rmdup_norg.bam
    samtools index ${params.sample}_merged_rmdup.bam
    printf "${params.sample} Read group added................................... DONE\n\n"
    """
}

process prepare_reference_for_gatk {
    label 'process_prepare_reference_gatk'
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

process gatk_indel_realigner {
    label 'process_gatk_indel_realigner'
    tag { "iteration_${iteration}" }

    publishDir "${params.outdir}/iteration_${iteration}/realigned_bams/", mode: 'symlink'

    container 'broadinstitute/gatk3:3.8-1'

    input:
    path(bam)
    path(bai)
    path(ref)
    path(fai)
    path(dict)
    val(iteration)

    output:
    path("${params.sample}_merged_rmdup.realigned.bam"), emit: realigned_bam

    script:
    def avail_mem = task.memory ? (task.memory.toGiga()).toInteger() : 4
    """
    java -Xmx${avail_mem}g -jar /usr/GenomeAnalysisTK.jar \\
        -T RealignerTargetCreator \\
        -R ${ref} \\
        -I ${bam} \\
        -o ${params.sample}.intervals

    java -Xmx${avail_mem}g -jar /usr/GenomeAnalysisTK.jar \\
        -T IndelRealigner \\
        -R ${ref} \\
        -I ${bam} \\
        -targetIntervals ${params.sample}.intervals \\
        -o ${params.sample}_merged_rmdup.realigned.bam
    """
}

process call_fixed_variants {
    label 'process_call_fixed_variants'
    tag { "iteration_${iteration}" }

    publishDir "${params.outdir}/iteration_${iteration}/variants/", mode: 'symlink'

    input:
    path(merged_bam)
    path(ref)
    val(iteration)

    output:
    path("${params.sample}*_homalt.bcf"), emit: filtered_bcf
    path("${params.sample}*_homalt.bcf.csi"), emit: filtered_bcf_index

    script:

    """
    set -euo pipefail
    ml -q bcftools

    bcftools mpileup -q ${params.mapQ} -Q ${params.baseQ} -B -f ${ref} ${merged_bam} --ignore-RG --threads ${task.cpus} -Ou | \\
        bcftools call -mv -Ob --threads ${task.cpus} -o ${params.sample}.bcf
    bcftools sort -Ob -o ${params.sample}_sorted.bcf ${params.sample}.bcf
    bcftools filter -i "DP>=${params.min_depth} & DP<${params.max_depth} & QUAL>=${params.variant_quality}" -Ob --threads ${task.cpus} \\
        -o ${params.sample}_DP${params.min_depth}-${params.max_depth}.bcf ${params.sample}_sorted.bcf
    bcftools filter -g ${params.snp_gap_indels} -Ob --threads ${task.cpus} -o ${params.sample}_DP${params.min_depth}-${params.max_depth}_g${params.snp_gap_indels}.bcf \\
        ${params.sample}_DP${params.min_depth}-${params.max_depth}.bcf
    bcftools filter -i 'GT="1/1"' -O b --threads ${task.cpus} -o ${params.sample}_DP${params.min_depth}-${params.max_depth}_g${params.snp_gap_indels}_homalt.bcf \\
        ${params.sample}_DP${params.min_depth}-${params.max_depth}_g${params.snp_gap_indels}.bcf
    bcftools index ${params.sample}_DP${params.min_depth}-${params.max_depth}_g${params.snp_gap_indels}_homalt.bcf
    """
}

process call_fixed_snvs {
    label 'process_call_fixed_variants'
    tag { "iteration_${iteration}" }

    publishDir "${params.outdir}/iteration_${iteration}/variants/", mode: 'symlink'

    input:
    path(merged_bam)
    path(ref)
    val(iteration)

    output:
    path("${params.sample}*_homalt.bcf"), emit: filtered_bcf
    path("${params.sample}*_homalt.bcf.csi"), emit: filtered_bcf_index

    script:
    """
    set -euo pipefail
    ml -q bcftools

    bcftools mpileup -q ${params.mapQ} -Q ${params.baseQ} -B -f ${ref} ${merged_bam} --ignore-RG --threads ${task.cpus} -Ou | \\
        bcftools call -mv -Ob --threads ${task.cpus} -o ${params.sample}.bcf
    bcftools sort -Ob -o ${params.sample}_sorted.bcf ${params.sample}.bcf
    bcftools filter -i "DP>=${params.min_depth} & DP<${params.max_depth} & QUAL>=${params.variant_quality}" -Ob --threads ${task.cpus} \\
        -o ${params.sample}_DP${params.min_depth}-${params.max_depth}.bcf ${params.sample}_sorted.bcf
    bcftools filter -g ${params.snp_gap_indels} -Ob --threads ${task.cpus} -o ${params.sample}_DP${params.min_depth}-${params.max_depth}_g${params.snp_gap_indels}.bcf \\
        ${params.sample}_DP${params.min_depth}-${params.max_depth}.bcf
    bcftools filter -i 'INDEL=0' -O b --threads ${task.cpus} -o ${params.sample}_DP${params.min_depth}-${params.max_depth}_g${params.snp_gap_indels}_no_indels.bcf \\
        ${params.sample}_DP${params.min_depth}-${params.max_depth}_g${params.snp_gap_indels}.bcf
    bcftools filter -i 'GT="1/1"' -O b --threads ${task.cpus} -o ${params.sample}_DP${params.min_depth}-${params.max_depth}_g${params.snp_gap_indels}_no_indels_homalt.bcf \\
        ${params.sample}_DP${params.min_depth}-${params.max_depth}_g${params.snp_gap_indels}_no_indels.bcf
    bcftools index ${params.sample}_DP${params.min_depth}-${params.max_depth}_g${params.snp_gap_indels}_no_indels_homalt.bcf
    """
}


process make_consensus {
    label 'process_make_consensus'
    tag { "iteration_${iteration}" }

    publishDir "${params.outdir}/iteration_${iteration}/reference/", mode: 'symlink'

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

    main:
    // index the genome
    bwa_index(ch_ref, ch_iteration)
    ch_ref_index = ch_ref.combine(bwa_index.out.index)

    // run mapping on all fastq files in parallel
    bwa_aln_samse(fastq_files, ch_ref_index, ch_iteration)
    ch_filtered_bams = bwa_aln_samse.out.sorted_rmdup_bam.collect()

    // merge and remove duplicates for all libraries
    merge_dedup_bams_all(ch_filtered_bams, ch_iteration)

    // prepare reference for GATK (fai + dict)
    prepare_reference_for_gatk(ch_ref, ch_iteration)

    // indel realignment on merged BAM before SNP calling
    gatk_indel_realigner(
        merge_dedup_bams_all.out.merged_rmdup_bam,
        merge_dedup_bams_all.out.merged_rmdup_bai,
        ch_ref,
        prepare_reference_for_gatk.out.fai,
        prepare_reference_for_gatk.out.dict,
        ch_iteration
    )
    ch_realigned_bam = gatk_indel_realigner.out.realigned_bam

    // call variants on the realigned BAM
    call_fixed_variants(ch_realigned_bam, ch_ref, ch_iteration)
    ch_filtered_bcf = call_fixed_variants.out.filtered_bcf.combine(call_fixed_variants.out.filtered_bcf_index)

    // make consensus sequence
    make_consensus(ch_filtered_bcf, ch_ref, ch_iteration)
    ch_consensus_fasta = make_consensus.out.consensus_fasta

    emit:
    consensus_fasta = ch_consensus_fasta
}
