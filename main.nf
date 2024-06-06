nextflow.enable.dsl=2

params.random_seed = 42
params.input_glob = '*{1,2}.fastq{,.gz}'
params.downsample = 100000
params.temporary_path = '/tmp'
params.bowtie2_index = null
params.genome_fasta = null

if (params.bowtie2_index == null && params.genome_fasta == null) {
    log.error("Both params.bowtie2_index and params.genome_fasta cannot be null. Please provide at least one.")
    System.exit(1)
}


include { seqtk_sample; mergeAndMarkDuplicates } from './modules/preprocessing'
include { bowtie2_create_index, bowtie2_align } from './modules/aligners'

input_channel = channel.fromFilePairs(params.input_glob)
                       .map {n -> [library: n[0], read1: n[1][0], read2: n[1][1]]}

workflow {
    // Preprocessing
    seqtk_done = seqtk_sample(input_channel)
    fastp_done = fastp(seqtk_done.out)
    // Alignment
    bowtie2_index_done = params.bowtie2_index == null ? bowtie2_create_index() : params.bowtie2_index
    bowtie2_input_channel = fastp_done.out.concat(bowtie2_index_done.index)
    bowtie2_done = bowtie2_align(bowtie2_input_channel)
    picard_MarkDuplicates()
    // Quality control
    samtools_flagstat

    picard_CollectGcBiasMetrics

    picard_CollectAlignmentSummaryMetrics

    picard_EstimateLibraryComplexity
   
    picard_CollectInsertSizeMetrics

    fastqc

    multiqc
}