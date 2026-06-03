# Iterative Mapping Pipeline

This pipeline is for building a sample-specific consensus reference from short-read sequencing data. The idea is simple: instead of mapping reads once to a distant reference genome, you map, call variants, update the reference with those variants, and repeat. Each round the reference gets closer to the actual sample, so reads map better and variants are called more accurately.

This is useful for ancient DNA or any situation where the sample is diverged enough from the reference that a single round of mapping leaves a lot on the table.

## What it does

The pipeline runs 8 mapping iterations in two phases:

**Phase 1 (iterations 1–3)** — SNVs only. It calls only substitutions and uses those to update the reference. Indels are excluded here on purpose because they shift coordinates and can cause instability early on. Getting the SNP landscape right first makes everything downstream more reliable.

**Phase 2 (iterations 4–8)** — SNPs + indels. Once the reference has converged on substitutions, indels are introduced. By this point the reference is close enough that indel calls are actually trustworthy.

Each iteration does the same thing: index the reference → map all FASTQs with BWA-MEM → filter by mapping quality and read length → remove duplicates → call variants per chromosome → merge variant calls → generate a new consensus FASTA. That consensus becomes the reference for the next iteration.

## Variant callers

We recommend DeepVariant (the default). It uses a neural network trained on real sequencing data and consistently produces more accurate variant calls than traditional methods, especially for ancient or low-coverage samples where read quality is variable. It runs in an Apptainer container so no module loading is needed.

GATK HaplotypeCaller is also supported (`--caller gatk`) but in our experience it does not perform as well as DeepVariant for this use case. It uses the cluster module system. Use it only if you have a specific reason to.

## Requirements

- Nextflow
- SLURM cluster (configured for NAISS/Dardel by default — change the partition and account in `nextflow.config`)
- BWA, SAMtools, bcftools available as modules
- Either GATK module or Apptainer with DeepVariant image
- Python 3 (for the duplicate removal script `samremovedup.py`)

## Input files

- FASTQ files: one or more `*.fastq.gz` files in a single directory
- Reference genome: a FASTA file to start from
- Chromosomes list: `chromosomes.tsv` — one chromosome/scaffold name per line, must match the sequence names in your reference FASTA

The chromosomes file controls which regions get variant calls. Edit it to match your reference. There is an example in the file for domestic cat.

## How to run

```bash
nextflow run main.nf \
  --sample SAMPLE_NAME \
  --outdir /path/to/output/ \
  --fastq_dir /path/to/fastqs/ \
  --reference /path/to/reference.fasta \
  --chromosomes chromosomes.tsv \
  --script_rmdup /path/to/samremovedup.py
```

Resume a run after failure or interruption:

```bash
nextflow run main.nf -resume [same parameters]
```

Resume is enabled by default in `nextflow.config`.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--sample` | (required) | Sample name, used in output file names |
| `--outdir` | (required) | Output directory |
| `--fastq_dir` | (required) | Directory containing input FASTQ files (`*.fastq.gz`) |
| `--reference` | (required) | Starting reference FASTA |
| `--caller` | `dv` | Variant caller: `dv` (DeepVariant) or `gatk` |
| `--phase1_iterations` | `3` | How many SNV-only iterations before switching to SNP+indel mode |
| `--mapQ` | `20` | Minimum mapping quality |
| `--baseQ` | `30` | Minimum base quality (GATK only) |
| `--min_readlength` | `40` | Minimum read length to keep after mapping |
| `--min_depth` | `3` | Minimum depth to call a variant |
| `--max_depth` | `50` | Maximum depth to call a variant (caps repetitive regions) |

## Output structure

Each iteration produces its own directory under `--outdir`:

```
outdir/
  iteration_1/
    raw_bams/        per-library sorted, filtered, deduplicated BAMs
    merged_bams/     all libraries merged into one BAM
    variants/        per-chromosome BCFs and the merged genome-wide BCF
    reference/       the consensus FASTA used as input for the next iteration
  iteration_2/
    ...
  iteration_8/
    reference/       final consensus — this is your end result
```

The final consensus FASTA from iteration 8 is your sample-specific reference. You can use that for any downstream analysis that benefits from a closer reference.

## Cluster configuration

The pipeline is configured to run on SLURM. The partition and account are set in `nextflow.config` — change them to match your allocation. Resource limits per process (CPUs, memory, time) are also defined there and can be adjusted if jobs are timing out or getting killed.

DeepVariant jobs are the most resource-heavy (32 CPUs, 64 GB RAM, 3h) because they run per chromosome. If your reference has many chromosomes, consider batching or increasing the queue size.
