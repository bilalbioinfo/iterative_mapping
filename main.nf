#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// --- params ---
params.sample = '' // sample name prefix
params.outdir = 'results'

params.fastq_dir = ''
params.reference = ''
params.script_rmdup = ''


// --- variant calling parameters ---
params.mapQ = 25
params.baseQ = 30
params.min_depth = 3
params.max_depth = 1000
params.variant_quality = 30
params.snp_gap_indels = 5


include { iterative_mapping as iterative_mapping_1 } from './iterative_mapping_workflow.nf'
include { iterative_mapping as iterative_mapping_2 } from './iterative_mapping_workflow.nf'

workflow {
    // first iteration
    ch_iter1 = Channel.value(1)
    ch_ref_1 = Channel.fromPath(params.reference)
    fastq_files_1 = Channel.fromPath("${params.fastq_dir}*fastq.gz")
    iterative_mapping_1(ch_ref_1, fastq_files_1, ch_iter1)

    // second iteration
    ch_iter2 = Channel.value(2)
    ch_ref_2 = iterative_mapping_1.out.consensus_fasta
    fastq_files_2 = Channel.fromPath("${params.fastq_dir}*fastq.gz")
    iterative_mapping_2(ch_ref_2, fastq_files_2, ch_iter2)

    // third iteration
    ch_iter3 = Channel.value(3)
    ch_ref_3 = iterative_mapping_2.out.consensus_fasta
    fastq_files_3 = Channel.fromPath("${params.fastq_dir}*fastq.gz")
    iterative_mapping_3(ch_ref_3, fastq_files_3, ch_iter3)

    // fourth iteration
    ch_iter4 = Channel.value(4)
    ch_ref_4 = iterative_mapping_3.out.consensus_fasta
    fastq_files_4 = Channel.fromPath("${params.fastq_dir}*fastq.gz")
    iterative_mapping_4(ch_ref_4, fastq_files_4, ch_iter4)

    // fifth iteration
    ch_iter5 = Channel.value(5)
    ch_ref_5 = iterative_mapping_4.out.consensus_fasta
    fastq_files_5 = Channel.fromPath("${params.fastq_dir}*fastq.gz")
    iterative_mapping_5(ch_ref_5, fastq_files_5, ch_iter5)

    // sixth iteration
    ch_iter6 = Channel.value(6)
    ch_ref_6 = iterative_mapping_5.out.consensus_fasta
    fastq_files_6 = Channel.fromPath("${params.fastq_dir}*fastq.gz")
    iterative_mapping_6(ch_ref_6, fastq_files_6, ch_iter6)
}


