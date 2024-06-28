gunzip small_ref.fa

# run tests and report
echo " test WITHOUT alignmnet index" > test_out.log
nextflow run ../main.nf --genome_fasta /mnt/home/aerijman/seqhimem02/FFPE_nextflow_workflow/test/small_ref.fa -with-conda -resume >> test_out.log 2>> test_out.log

echo "test WITH alignment index" >> test_out.log
nextflow run ../main.nf --genome_fasta /mnt/home/aerijman/seqhimem02/FFPE_nextflow_workflow/test/small_ref.fa \
                     --bowtie2_index /mnt/home/aerijman/seqhimem02/FFPE_nextflow_workflow/test/reference_bowtie2_index/ -with-conda -resume >> test_out.log 2>> test_out.log

rm -r reference_bowtie2_index/ r/ test_read/ work/
