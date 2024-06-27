process bowtie2_create_index {
    tag {"creating index"}
    conda "bioconda:bowtie2=2.5.3"
    publishDir "reference_bowtie2_index/"
    cpus 24

    input:
        path(genome_fasta)

    output:
    // nf-core use the name of the path as index though...
    path("*.bt2")

    shell:
    '''
    index_name=$(basename !{genome_fasta} | awk -F"/" '{print $NF}'| cut -f 1 -d '.')
    bowtie2-build --threads !{task.cpus} !{genome_fasta} ${index_name}
    # touch testing.1.bt2 testing.2.bt2 testing.3.bt2
    '''
}

process bowtie2_align {
    cpus 8
    tag {library}
    conda "bioconda:bowtie2=2.5.3 bioconda:samtools=1.19.2"
    publishDir "${library}/"

    input:
    tuple val(library),
          path(read1),
          path(read2)
    path(index)

    output:
    tuple val(library), 
          path("*bam")

    shell:
    '''
    echo !{index}    

    INDEX=$(find -L !{index} -name "*.1.bt2" | sed 's/\\.1.bt2\$//' | grep -v "rev")

    bowtie2  -p !{task.cpus} -x $INDEX -1 !{read1} -2 !{read2} \
    | samtools sort --no-PG -@!{task.cpus} -T !{params.temporary_path} \
                    -O bam -o !{library}.sorted.bam
    '''
}

process picard_MarkDuplicates{
    cpus 8
    tag {library}
    conda "bioconda::picard=3.1.1"
    publishDir "${library}/"

    input:
    tuple val(library),
          path(bam)

    output:
    tuple val(library), 
          path("*md.bam")

    shell:
    '''
    _JAVA_OPTIONS="-Xmx2048m -Xms256m" 
    export _JAVA_OPTIONS && \
    picard MarkDuplicates  --INPUT !{bam} --OUTPUT !{library}.md.bam  \
    --METRICS_FILE !{library}.md.txt  --REMOVE_DUPLICATES true --ASSUME_SORTED true  \
    --DUPLICATE_SCORING_STRATEGY SUM_OF_BASE_QUALITIES  \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE "100"   --VALIDATION_STRINGENCY LENIENT \
    --TAGGING_POLICY All --QUIET true --VERBOSITY ERROR
    '''
}
