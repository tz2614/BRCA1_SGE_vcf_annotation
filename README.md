# BRCA1_SGE_annotator for use on C01 samples from clinical pool pipeline
  - By Tengyue Zheng
  - 24/10/2019

## Description
  Annotate vcfs containing BRCA1 variants present in exon2-5 and exon15-23 with function scores and class.
  Note: At the time of writing, because the clinical pool pipeline is planning to be moved to the trust cluster and changed to snakemake workflow, the software will be made to work in python3.5, and will be made available for integration with the new clinical pool pipeline instead.

## Getting Started
  1. Start by cloning the git repository into your working directory
  2. To clone git repository go here: https://gitlab.com/cuhbioinformatics/hiv_pipeline.git
  3. The vcfs files initially used for testing are available upon request by email: tengyue.zheng@addenbrookes.nhs.uk.

## Abbreviations

  1. work_dir - working directory where the software is installed
  2. runfolder - directory where the vcfs are stored
  3. root - path of parent directories to working directory or runfolder

## Prerequisites
  1. You need to have following software installed before you can run your script
  
  - python version 3.5
  - samtools version 1.5
  - bcftools version 1.7 

  2. for samtools and bcftools go to www.htslib.org/download and follow the general instructions to download the tar.bz2 file.

      link to the above versions are as follows:
  - samtools version 1.5 - https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2
  - bcftools version 1.7 - https://github.com/samtools/bcftools/releases/download/1.7/bcftools-1.7.tar.bz2

  3. for instructions on how to install bcftools version 1.7 
  - URL: https://cuhbioinformatics.atlassian.net/wiki/spaces/BT/pages/577306625/Software+installation+and+version+control+on+Trust+Cluster

## Main scripts:
  - BRCA1_SGE_ref.py
  - BRCA2_SGE_vcf_annotator.py
  - compare_CP2SGE.py
  - remove_INFO_field.py

## Unit testing scripts:

  - test_BRCA1_SGE_ref.py
  - test_BRCA1_SGE_vcf_annotator.py

## Reference csv file containing BRCA1_SGE annotations

  - /work_dir/data/41586_2018_461_MOESM3_ESM.csv

## User Requirements:
  Add function scores and class for BRCA1 variants from the 41586_2018_461_MOESM3_ESM.csv to C01####.vcf files. 
  Test the annotation tool on vcfs obtained from /data/CP/ on the university cluster and generate a new annotated vcf ready for analysis.

  The input file is vcfs that begin with C01####.vcf

  The output files for each vcf should include:

  1. C01####.vcf.gz
  2. C01####.vcf.gz.tbi
  3. C01####.BRCA1_annotated.vcf

  within the runfolder the vcf is located.

## Instructions:

To run the annotation software follow the instructions below:

1. check to the content of the git repository you have cloned, make sure it includes the main scripts and testing scripts above, and
   the reference files containing the BRCA1 variants with function scores and class.

2. Navigate to the runfolder where the vcf files are stored, the files should be in the format C01####.vcf

3. Load the software required on the trust cluster to run the software

```Bash
$ ml python3.5 samtools/1.5 bcftools/1.7
``` 

4. To annotate a vcf for a specific runfolder execute the following

```Bash
$  python3.5 work_dir/bin/BRCA1_SGE_vcf_annotator.py work_dir/data/41586_2018_461_MOESM3_ESM.csv /path/to/<sample>.vcf
```

5. To check that the annotation is being processed by assessing the output on command line

It should display the information about annotation and what files are being generated see below:

```Bash
$ creating header file
$ header file already exists
$ sorted reference vcf already exists
$ bgzip file /root/runfolder/BRCA1_SGE_ref.BRCA1_annotated.sorted.vcf.gz already exists
$ tabix file /root/runfolder/BRCA1_SGE_ref.BRCA1_annotated.sorted.vcf.gz.tbi already exists
$ merged_bgzip file already exists
$ merged tabix file already exists
$ bgzipping /root/work_dir/<sample>.vcf
$ indexing /root/work_dir/<sample>.vcf.gz
$ create new vcf with additional BRCA1 annotation
$ annotated vcf: /root/runfolder/<sample>.BRCA1_annotated.vcf
```

6. When the annotation is complete, check that the output files have been generated.
   If in doubt consult Tengyue Zheng or senior member of the bioinformatics team.

For a break down of the workflow of the pipeline please see workflow_to_annotate_BRCA1_variants(for_reference_file).png and workflow_to_annotate_BRCA1_variants(for_vcf_file).png in the repository
