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
    path("${fq.baseName}_sorted.bam"), emit: sorted_bam

    script:
    def reference = ref_index[0]
    """
    set -euo pipefail
    ml bwa samtools

    bwa aln -l 16500 -n 0.01 -o 2 -t ${task.cpus} ${reference} ${fq} | \
        bwa samse ${reference} - ${fq} | \
        samtools view -F 4 -q 1 -@ ${task.cpus} -bh - | \
        samtools sort -@ ${task.cpus} -o ${fq.baseName}_sorted.bam -
    samtools index ${fq.baseName}_sorted.bam
    """
}

process merge_dedup_bams {
    label 'process_merge_dedup_bams'
    tag  { "iteration_${iteration}" }

    publishDir "${params.outdir}/iteration_${iteration}/merged_bams/", mode: 'symlink'

    input:
    path(sorted_bams)
    val(iteration)

    output:
    path("${params.sample}_merged_rmdup.bam"), emit: merged_rmdup_bam

    script:
    """
    set -euo pipefail
    ml -q samtools/1.20

    ############################################################################################
    ### P13117_D1
    printf "Merging for library P13117_D1...\n"
    samtools merge -@ ${task.cpus} -o P13117_D1_merged.bam \
        P13117_P13117_1001*bam \
        P13117_P13117_1002*bam \
        P13117_P13117_1003*bam \
        P13117_P13117_1004*bam \
        P13117_P13117_1005*bam \
        P13117_P13117_1006*bam \
        P13117_P13117_1007*bam \
        P13117_P13117_1008*bam \
        P13117_P13117_1009*bam \
        P13117_P13117_1010*bam


    printf "Removing duplicates for library P13117_D1...\n"
    samtools view -@ ${task.cpus} -h P13117_D1_merged.bam | \
    python3 ${params.script_rmdup} | \
    samtools view -@ ${task.cpus} -bh -o P13117_D1_merged_rmdup.bam -
    samtools index P13117_D1_merged_rmdup.bam
    printf "Merging and removing duplicates for library P13117_D1... DONE\n\n"

    ############################################################################################
    ### P15057_D1
    printf "Merging for library P15057_D1...\n"
    samtools merge -@ ${task.cpus} -o P15057_D1_merged.bam \
        P15057_P15057_1169*bam \
        P15057_P15057_1175*bam \
        P15057_P15057_1181*bam \
        P15057_P15057_1187*bam \
        P15057_P15057_1193*bam \
        P15057_P15057_1199*bam \
        P15057_P15057_1205*bam \
        P15057_P15057_1211*bam \
        P15057_P15057_1217*bam \
        P15057_P15057_1223*bam \
        P15057_P15057_1227*bam \
        P15057_P15057_1229*bam


    printf "Removing duplicates for library P15057_D1...\n"
    samtools view -@ ${task.cpus} -h P15057_D1_merged.bam | \
    python3 ${params.script_rmdup} | \
    samtools view -@ ${task.cpus} -bh -o P15057_D1_merged_rmdup.bam -
    samtools index P15057_D1_merged_rmdup.bam
    printf "Merging and removing duplicates for library P15057_D1... DONE\n\n"

    ############################################################################################
    ### P15057_D2
    printf "Merging for library P15057_D2...\n"
    samtools merge -@ ${task.cpus} -o P15057_D2_merged.bam \
        P15057_P15057_1110*bam \
        P15057_P15057_1116*bam \
        P15057_P15057_1122*bam \
        P15057_P15057_1128*bam \
        P15057_P15057_1134*bam \
        P15057_P15057_1140*bam \
        P15057_P15057_1146*bam \
        P15057_P15057_1152*bam \
        P15057_P15057_1158*bam \
        P15057_P15057_1163*bam \
        P15057_P15057_1165*bam \
        P15057_P15057_1167*bam


    printf "Removing duplicates for library P15057_D2...\n"
    samtools view -@ ${task.cpus} -h P15057_D2_merged.bam | \
    python3 ${params.script_rmdup} | \
    samtools view -@ ${task.cpus} -bh -o P15057_D2_merged_rmdup.bam -
    samtools index P15057_D2_merged_rmdup.bam
    printf "Merging and removing duplicates for library P15057_D2... DONE\n\n"

    ############################################################################################
    ### P29310_D1
    printf "Merging for library P29310_D1...\n"
    samtools merge -@ ${task.cpus} -o P29310_D1_merged.bam \
        P29310_P29310_1001*bam \
        P29310_P29310_1002*bam \
        P29310_P29310_1003*bam \
        P29310_P29310_1004*bam \
        P29310_P29310_1005*bam \
        P29310_P29310_1006*bam \
        P29310_P29310_1007*bam \
        P29310_P29310_1008*bam

    printf "Removing duplicates for library P29310_D1...\n"
    samtools view -@ ${task.cpus} -h P29310_D1_merged.bam | \
    python3 ${params.script_rmdup} | \
    samtools view -@ ${task.cpus} -bh -o P29310_D1_merged_rmdup.bam -
    samtools index P29310_D1_merged_rmdup.bam
    printf "Merging and removing duplicates for library P29310_D1... DONE\n\n"

    ############################################################################################
    ### P29310_D2
    printf "Merging for library P29310_D2...\n"
    samtools merge -@ ${task.cpus} -o P29310_D2_merged.bam \
        P29310_P29310_1009*bam \
        P29310_P29310_1010*bam \
        P29310_P29310_1011*bam \
        P29310_P29310_1012*bam \
        P29310_P29310_1013*bam \
        P29310_P29310_1014*bam \
        P29310_P29310_1015*bam \
        P29310_P29310_1016*bam

    printf "Removing duplicates for library P29310_D2...\n"
    samtools view -@ ${task.cpus} -h P29310_D2_merged.bam | \
    python3 ${params.script_rmdup} | \
    samtools view -@ ${task.cpus} -bh -o P29310_D2_merged_rmdup.bam -
    samtools index P29310_D2_merged_rmdup.bam
    printf "Merging and removing duplicates for library P29310_D2... DONE\n\n"

    ############################################################################################
    ### Merging all libraries
    printf "Merging for all libraries...\n"
    samtools merge -@ ${task.cpus} -o ${params.sample}_merged_rmdup.bam \
        P13117_D1_merged_rmdup.bam \
        P15057_D1_merged_rmdup.bam \
        P15057_D2_merged_rmdup.bam \
        P29310_D1_merged_rmdup.bam \
        P29310_D2_merged_rmdup.bam
    samtools index ${params.sample}_merged_rmdup.bam
    printf "Merging and removing duplicates for all libraries... DONE\n\n"
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
    path("${params.sample}_consensus.fasta"),   emit: consensus_fasta
    path("${params.sample}_consensus.log"),     emit: consensus_log

    script:
    """
    set -euo pipefail
    ml -q bcftools
    bcftools consensus -f ${ref} -o ${params.sample}_consensus.fasta ${filtered_bcf} 2> ${params.sample}_consensus.log
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
    ch_filtered_bams = bwa_aln_samse.out.sorted_bam.collect()

    // merge and remove duplicates for test dataset
    merge_dedup_bams(ch_filtered_bams, ch_iteration)
    ch_dedup_bams = merge_dedup_bams.out.merged_rmdup_bam

    // call variants on the merged BAM
    call_fixed_variants(ch_dedup_bams, ch_ref, ch_iteration)
    ch_filtered_bcf = call_fixed_variants.out.filtered_bcf.combine(call_fixed_variants.out.filtered_bcf_index)

    // make consensus sequence
    make_consensus(ch_filtered_bcf, ch_ref, ch_iteration)
    ch_consensus_fasta = make_consensus.out.consensus_fasta

    emit:
    consensus_fasta = ch_consensus_fasta
}
