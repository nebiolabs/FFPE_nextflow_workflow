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


include { seqtk_sample; fastp } from './modules/preprocessing'
include { bowtie2_create_index; bowtie2_align; picard_MarkDuplicates } from './modules/aligners'
include { samtools_flagstat; picard_CollectGcBiasMetrics; picard_EstimateLibraryComplexity; picard_CollectAlignmentSummaryMetrics; picard_CollectInsertSizeMetrics; fastqc; bedtools_genome_coverage; multiqc; tasmanian } from  './modules/qc'

input_channel = channel.fromFilePairs(params.input_glob)
                       .map {n -> [library: n[0], read1: n[1][0], read2: n[1][1]]}

workflow {
    // Preprocessing
    seqtk_done = seqtk_sample(input_channel)
    fastp_done = fastp(seqtk_done)
    
    // Alignment
    genome_channel = channel.fromPath(params.genome_fasta)
    bowtie2_index_done = params.bowtie2_index == null ? bowtie2_create_index(genome_channel) : channel.fromPath(params.bowtie2_index)

    bowtie2_done                            = bowtie2_align(fastp_done, bowtie2_index_done)
    markduplicates_done                     = picard_MarkDuplicates(bowtie2_done)
    
    // bedfile
    bedtools_done                           = bedtools_genome_coverage(markduplicates_done)

    // Quality control
    flagstats_done                          = samtools_flagstat(markduplicates_done)
    gcbias_done                             = picard_CollectGcBiasMetrics(markduplicates_done)
    align_summary_done                      = picard_CollectAlignmentSummaryMetrics(markduplicates_done)
    estimate_complexity_done                = picard_EstimateLibraryComplexity(markduplicates_done)       
    insert_size_done                        = picard_CollectInsertSizeMetrics(markduplicates_done)
    fastqc_done                             = fastqc(markduplicates_done)

    // MultiQC
    multiqc_ch = flagstats_done.groupTuple(by: [0, 1])
        .join(fastqc_done.groupTuple(by: [0, 1]))
        .join(markduplicates_done.groupTuple(by: [0, 1]))
        .join(insert_size_done.for_multiqc.groupTuple(by: [0, 1]))
        .join(gcbias_done.groupTuple(by: [0, 1]))
        .join(align_summary_done.groupTuple(by: [0, 1]))

    multiqc_done                            = multiqc(multiqc_ch)
    tasmanian_done                          = tasmanian(markduplicates_done)
}   
