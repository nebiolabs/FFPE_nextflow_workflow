samtools_flagstat{
    cpus 8
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


picard_CollectGcBiasMetrics{
    cpus 8
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
    --REFERENCE_SEQUENCE !{params.genome} --ASSUME_SORTED true  --VALIDATION_STRINGENCY LENIENT \
    --QUIET true --VERBOSITY ERROR
    '''
}

picard_CollectAlignmentSummaryMetrics{
    cpus 8
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
    picard CollectAlignmentSummaryMetrics --INPUT !{ban} --OUTPUT !{library}.alignment_summary.out \
    --MAX_INSERT_SIZE 100000 --METRIC_ACCUMULATION_LEVEL ALL_READS --IS_BISULFITE_SEQUENCED false  \
    --REFERENCE_SEQUENCE !{params.genome}  --ASSUME_SORTED true  --VALIDATION_STRINGENCY LENIENT \
    --QUIET true --VERBOSITY ERROR
    '''
}

picard_EstimateLibraryComplexity {
    cpus 8
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

picard_CollectAlignmentSummaryMetrics { 
    cpus 8
    tag {library}
    conda "bioconda::picard=3.1.1"
    publishDir "${library}/"

    input:
    tuple val(library),
          path(bam)

    output:
    tuple val(library), 
          path("*insert_size_metrics.out")

    shell:
    '''
    _JAVA_OPTIONS="-Xmx2048m -Xms256m" 
    export _JAVA_OPTIONS && \
    picard CollectInsertSizeMetrics --INPUT !{bam} --OUTPUT !{library}.insert_size_metrics.out \
    --Histogram_FILE !{library}.insert_size_hist.out --DEVIATIONS 10.0   --MINIMUM_PCT 0.05 \
    --REFERENCE_SEQUENCE !{params.genome} --ASSUME_SORTED true --METRIC_ACCUMULATION_LEVEL ALL_READS  \
    --VALIDATION_STRINGENCY LENIENT --QUIET true --VERBOSITY ERROR
    '''
}

fastqc { 
    cpus 8
    tag {library}
    conda "bioconda::fastqc=0.12.1"
    publishDir "${library}/"

    input:
    tuple val(library),
          path(bam)

    output:
    tuple val(library), 
          path("*output.txt ")

    shell:
    '''
    _JAVA_OPTIONS="-Xmx2048m -Xms256m" 
    export _JAVA_OPTIONS && \
    mkdir fastqc_tmp_dir
    fastqc --outdir fastqc_tmp_dir --threads !{task.cpus} --quiet --extract  --kmers 7 -f bam !{bam}  \
    && cp fastqc_tmp_dir/*/fastqc_data.txt !{library}.fastqc.output.txt \
    && cp fastqc_tmp_dir/*\.html !{library}.fastqc.output.html
    '''
}

bedtools_genome_coverage { 
    cpus 8
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
multiqc { 
    cpus 8
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
          path("*.multiqc.report")

    shell:
    '''
    mkdir multiqc_WDir && \
    mkdir multiqc_WDir/samtools_0 && \
    mkdir multiqc_WDir/samtools_0/flagstat_0 && \
    ln -s !{flagstat} multiqc_WDir/samtools_0/flagstat_0/  && \
    mkdir multiqc_WDir/fastqc_1 && \
    mkdir multiqc_WDir/fastqc_1/data_0 && \
    mkdir multiqc_WDir/fastqc_1/data_0/file_0 && \
    ln -s !{fastqc} multiqc_WDir/fastqc_1/data_0/file_0/ && \
    mkdir multiqc_WDir/picard_2 && \
    mkdir multiqc_WDir/picard_2/markdups_0 && \
    ln -s !{markduplicates} multiqc_WDir/picard_2/markdups_0/  && \
    mkdir multiqc_WDir/picard_3 && \
    mkdir multiqc_WDir/picard_3/insertsize_0 && \
    ln -s !{insertsize} multiqc_WDir/picard_3/insertsize_0/  && \
    mkdir multiqc_WDir/picard_4 && \
    mkdir multiqc_WDir/picard_4/gcbias_0 && \
    ln -s !{gcbias} multiqc_WDir/picard_4/gcbias_0/  && \
    mkdir multiqc_WDir/picard_5 && \
    mkdir multiqc_WDir/picard_5/alignment_metrics_0 && \
    ln -s !{alignmentsummary} multiqc_WDir/picard_5/alignment_metrics_0/ && \
    multiqc multiqc_WDir --filename !{library}.multiqc.report
    '''
}