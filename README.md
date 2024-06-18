This repository contains tools and data related to DNA libraries preserved with [formalin-fixed, paraffin-embedded (FFPE)](https://www.neb.com/en-us/products/next-generation-sequencing-library-preparation/ffpe-dna/ffpe-dna-repair)

The 3 [nextflow](https://www.nextflow.io/) scripts are written as modules:
 - preprocessing.nf (to trim and downsample reads as necessary)
 - aligners.nf (to align reads and mark duplicates)
 - qc.nf (to collect quality control metrics from a variety of tools including [Tasmanian](https://github.com/nebiolabs/tasmanian-mismatch))
