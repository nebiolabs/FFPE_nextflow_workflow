process seqtk_sample {
    cpus 4
    tag {library}
    conda "bioconda::seqtk=1.4"
    publishDir "${library}/"

    input:
    tuple val(library),
          path(read1),
          path(read2)

    output:
    tuple val(library), 
          path("*R1.ds.fastq"),
          path("*R2.ds.fastq")

    shell:
    '''
    seqtk sample -s!{params.random_seed} !{read1} !{params.downsample} > !{library}_R1.ds.fastq &
    seqtk sample -s!{params.random_seed} !{read2} !{params.downsample} > !{library}_R2.ds.fastq &
    wait
    '''
}

process fastp {
    cpus 4
    tag {library}
    conda "bioconda::fastp=0.23.2"
    publishDir "${library}/"

    input:
    tuple val(library),
          path(read1),
          path(read2)

    output:
    tuple val(library),
          path("*R1.ds.trimmed.fastq"),
          path("*R2.ds.trimmed.fastq")

    shell:
    '''
    fastp --thread !{task.cpus} \
          --report_title 'fastp report for !{library}' \
          -i !{read1} \
          -o !{library}.R1.ds.trimmed.fastq \
          -I !{read2} \
          -O !{library}.R2.ds.trimmed.fastq \
          --detect_adapter_for_pe 

    '''
}
