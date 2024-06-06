nextflow.enable.dsl=2

params.random_seed = 42
params.input_glob = '*{1,2}.fastq{,.gz}'
params.downsample = 100000
params.temporary_path = '/tmp'
params.bowtie2_index = null
params.genome_fasta = null

//if (params.bowtie2_index == null && params.genome_fasta == null) {
//    log.error("Both params.bowtie2_index and params.genome_fasta cannot be null. Please provide at least one.")
//    System.exit(1)
//}


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

    bowtie2_input_channel                   = fastp_done.out.combine(bowtie2_index_done.index)
    bowtie2_done                            = bowtie2_align(bowtie2_input_channel)
    markduplicates_done                     = picard_MarkDuplicates(bowtie2_done.out)
    
    // bedfile
    bedtools_done                           = bedtools_genome_coverage(markduplicates_done.out.bam)

    // Quality control
    flagstats_done                          = samtools_flagstat(markduplicates_done)
    gcbias_done                             = picard_CollectGcBiasMetrics(markduplicates_done)
    align_summary_done                      = picard_CollectAlignmentSummaryMetrics(markduplicates_done)
    etimater_complexity_done                = picard_EstimateLibraryComplexity(markduplicates_done)       
    insert_size_done                        = picard_CollectInsertSizeMetrics(markduplicates_done)
    fastqc_done                             = fastqc(markduplicates_done.out)

    // MultiQC
    multiqc_ch = flagstats_done.out
        .join(fastqc_done.out)
        .join(markduplicates_done.out)
        .join(insert_size_done.out)
        .join(gcbias_done.out)
        .join(align_summary_done.out)
    multiqc_done                            = multiqc(multiqc_ch)
}