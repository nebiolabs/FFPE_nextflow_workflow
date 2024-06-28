process samtools_flagstat{
    cpus 4
    tag {library}
    conda "bioconda::samtools=1.15.1"
    publishDir "${library}/"

    input:
    tuple val(library),
          path(bam)

    output:
    tuple val(library), 
          path("*flagstat.tsv")

    shell:
    '''
    samtools flagstat -@{task.cpus} !{bam} > !{library}.flagstat.tsv
    '''
}


process picard_CollectGcBiasMetrics{
    cpus 4
    tag {library}
    conda "bioconda::picard=3.1.1"
    publishDir "${library}/"

    input:
    tuple val(library),
          path(bam)

    output:
    tuple val(library), 
          path("*gcbias.out")

    shell:
    '''
    _JAVA_OPTIONS="-Xmx2048m -Xms256m" 
    export _JAVA_OPTIONS && \
    picard CollectGcBiasMetrics  --INPUT !{bam} --OUTPUT !{library}.gcbias.out  \
    --CHART_OUTPUT !{library}.gc_chart.pdf --SUMMARY_OUTPUT !{library}.gc_summary.out \
    --WINDOW_SIZE 100 --MINIMUM_GENOME_FRACTION 5e-05 --IS_BISULFITE_SEQUENCED true \
    --REFERENCE_SEQUENCE !{params.genome_fasta} --ASSUME_SORTED true  --VALIDATION_STRINGENCY LENIENT \
    --QUIET true --VERBOSITY ERROR
    '''
}

process picard_CollectAlignmentSummaryMetrics{
    cpus 4
    tag {library}
    conda "bioconda::picard=3.1.1"
    publishDir "${library}/"

    input:
    tuple val(library),
          path(bam)

    output:
    tuple val(library), 
          path("*alignment_summary.out")

    shell:
    '''
    _JAVA_OPTIONS="-Xmx2048m -Xms256m" 
    export _JAVA_OPTIONS && \
    picard CollectAlignmentSummaryMetrics --INPUT !{bam} --OUTPUT !{library}.alignment_summary.out \
    --MAX_INSERT_SIZE 100000 --METRIC_ACCUMULATION_LEVEL ALL_READS --IS_BISULFITE_SEQUENCED false  \
    --REFERENCE_SEQUENCE !{params.genome_fasta}  --ASSUME_SORTED true  --VALIDATION_STRINGENCY LENIENT \
    --QUIET true --VERBOSITY ERROR
    '''
}

process picard_EstimateLibraryComplexity {
    cpus 4
    tag {library}
    conda "bioconda::picard=3.1.1"
    publishDir "${library}/"

    input:
    tuple val(library),
          path(bam)

    output:
    tuple val(library), 
          path("*library_complexity.out")

    shell:
    '''
    _JAVA_OPTIONS="-Xmx2048m -Xms256m" 
    export _JAVA_OPTIONS && \
    picard EstimateLibraryComplexity  --INPUT !{bam} --OUTPUT !{library}.library_complexity.out  \
    --MIN_IDENTICAL_BASES 5 --MAX_DIFF_RATE 0.03 --MIN_MEAN_QUALITY 20 --MAX_GROUP_RATIO 500 \
    --READ_NAME_REGEX \'[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.\' --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100  \
    --VALIDATION_STRINGENCY LENIENT --QUIET true --VERBOSITY ERROR
    '''
}

process picard_CollectInsertSizeMetrics {
    cpus 4
    tag {library}
    conda "bioconda::picard=3.1.1"
    publishDir "${library}/"

    input:
    tuple val(library),
          path(bam)

    output:
    tuple val(library),
          path("*insert_size_metrics.out"),
          emit: for_multiqc

    tuple val(library),
          path("histogram.txt"),
          emit: other_data

    shell:
    '''
    _JAVA_OPTIONS="-Xmx2048m -Xms256m"
    export _JAVA_OPTIONS && \
    picard CollectInsertSizeMetrics --INPUT !{bam} --OUTPUT !{library}.insert_size_metrics.out --Histogram_FILE histogram.txt \
    --DEVIATIONS 10.0 --MINIMUM_PCT 0.05 --REFERENCE_SEQUENCE !{params.genome_fasta} --ASSUME_SORTED true \
    --METRIC_ACCUMULATION_LEVEL ALL_READS  --VALIDATION_STRINGENCY LENIENT --QUIET true --VERBOSITY ERROR
    '''
}

process fastqc { 
    cpus 4
    tag {library}
    conda "bioconda::fastqc=0.12.1"
    publishDir "${library}/"

    input:
    tuple val(library),
          path(bam)

    output:
    tuple val(library), 
          path("*output.txt")

    shell:
    '''
    _JAVA_OPTIONS="-Xmx2048m -Xms256m" 
    export _JAVA_OPTIONS && \
    mkdir fastqc_tmp_dir
    fastqc --outdir fastqc_tmp_dir --threads !{task.cpus} --quiet --extract  --kmers 7 -f bam !{bam}  \
    && cp fastqc_tmp_dir/*/fastqc_data.txt !{library}.fastqc.output.txt \
    && cp fastqc_tmp_dir/*.html !{library}.fastqc.output.html
    '''
}

process bedtools_genome_coverage { 
    cpus 4
    tag {library}
    conda "bioconda::bedtools=2.30"
    publishDir "${library}/"

    input:
    tuple val(library),
          path(bam)

    output:
    tuple val(library), 
          path("*.bed")

    shell:
    '''
    _JAVA_OPTIONS="-Xmx2048m -Xms256m" 
    export _JAVA_OPTIONS && \
    bedtools genomecov -ibam !{bam} -max 100000 > !{library}.bed
    '''
}

process multiqc { 
    cpus 4
    tag {library}
    conda "bioconda::multiqc=1.11"
    publishDir "${library}/"

    input:
    tuple val(library),
          path(flagstat),
          path(fastqc),
          path(markdup),
          path(insertsize),
          path(gcbias),
          path(alignmentsummary)

    output:
    tuple val(library), 
          path("*.multiqc.report_data")

    shell:
    '''
    multiqc . --filename !{library}.multiqc.report
    '''
}

process tasmanian {
    cpus 1
    tag {library}
    conda "bioconda::tasmanian-mismatch=1.0.7 bioconda::samtools=1.13"
    publishDir "${library}/"

    input:
    tuple val(library),
          path(bam)

    output:
    tuple val(library),
          path("*.html"),
          emit: html

    tuple val(library),
          path("*.csv"),
          emit: table

    shell:
    '''
    samtools view !{bam} | run_tasmanian -r !{params.genome_fasta} > !{library}.tasmanian.output.csv
    '''
}
