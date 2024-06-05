#!/usr/bin/env nextflow

params.reads1 = "" // Input reads 1
params.reads2 = "" // Input reads 2
params.human_genome = "" // Path to human genome index for bowtie2
params.output_dir = "results" // Output directory

process seqtk_sample {
    input:
    path reads1 from params.reads1
    path reads2 from params.reads2

    output:
    path "*_sampled.fastq.gz"

    script:
    """
    seqtk sample ${reads1} 1000 > reads1_sampled.fastq.gz
    seqtk sample ${reads2} 1000 > reads2_sampled.fastq.gz
    """
}

process fastp {
    input:
    path reads from seqtk_sample.out

    output:
    path "*_fastp.fastq.gz"

    script:
    """
    fastp -i reads1_sampled.fastq.gz -I reads2_sampled.fastq.gz -o reads1_fastp.fastq.gz -O reads2_fastp.fastq.gz
    """
}

process bowtie2 {
    input:
    path reads from fastp.out
    path index from params.human_genome

    output:
    path "*.sam"

    script:
    """
    bowtie2 -x ${index} -1 reads1_fastp.fastq.gz -2 reads2_fastp.fastq.gz -S output.sam
    """
}

process markDuplicates {
    input:
    path samfile from bowtie2.out

    output:
    path "*.dedup.bam"

    script:
    """
    picard MarkDuplicates I=${samfile} O=output.dedup.bam M=marked_dup_metrics.txt
    """
}
