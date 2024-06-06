process bowtie2_create_index {
    tag {"creating index"}
    conda "bioconda:bowtie2=2.5.3"
    publishDir "${bt2_index}/"

    input:
    path(genome_fasta)

    output:
    // nf-core use the name of the path as index though...
    path("*.bt2", emit: 'index')

    shell:
    '''
    index_name=$(basename !{params.genome_fasta} | cut -f 1 -d '.')
    bowtie2-build !{params.genome_fasta} ${index_name}
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
          path(read2),
          path(index)

    output:
    tuple val(library), 
          path("*bam")

    shell:
    '''
    bowtie2  -p !{params.cpus} -x !{index} -1 read1 -2 read2 \
    | samtools sort --no-PG -@!{task.cpus} -T !{params.temporary_path} \
                    -O bam -o !{library}.sorted.bam
    '''
}

picard_MarkDuplicates{
    cpus 8
    tag {library}
    conda "bioconda::"
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

    bowtie2  -p !{params.cpus}
             -x '/cvmfs/data.galaxyproject.org/byhand/CHM13_T2T_v2.0/bowtie2_index/CHM13_T2T_v2.0/CHM13_T2T_v2.0'   
             -1 'input_f.fastq' -2 'input_r.fastq' | \
             samtools sort --no-PG -@${GALAXY_SLOTS:-2} -T "${TMPDIR:-.}" -O bam -o output.bam

    '''
}