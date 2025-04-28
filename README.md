# nf-larry

nf-larry is the Nextflow implementation of LARRY analysis workflow that utilises CBUtools package which is a private package. 




* QC metrics to output
* Some standardised plots
* Additional output as notebook or HTML like 10x cell ranger (summary report)
* 









## How to run:

The recommended way to use nextflow is to run it in a screen session. These steps can be directly used in Sanger's FARM, but you can modify each step according to the environment you're working on or the job scheduler your HPC uses:

1. Start a screen session: `screen -S nf_run1`
2. Start a small interactive job for nextflow: `bsub -G cellgeni -n1 -R"span[hosts=1]" -Is -q long -R"select[mem>2000] rusage[mem=2000]" -M2000 bash`
3. Modify one of RESUME scripts in examples folder (pre-made Nextflow run scripts)
4. Run the RESUME scripts you modified: `./RESUME-scautoqc-all`
5. You can leave your screen session and let it run in the background: `Ctrl+A, D`
6. You can go back to the screen session by: `screen -x nf_run1`

## Files:

* `main.nf` - the Nextflow pipeline that executes scAutoQC pipeline.
* `nextflow.config` - the configuration script that allows the processes to be submitted to IBM LSF on Sanger's HPC and ensures correct environment is set via singularity container (this is an absolute path). Global default parameters are also set in this file and some contain absolute paths.
* `bin/` - a folder that includes Python scripts used in each step
  * `find_larry_seqs.py`: 
* `RESUME-larry`
* `Dockerfile` - a dockerfile to reproduce the environment used to run the pipeline.

## Workflow

![](images/nf-cellsnplite-light.png#gh-light-mode-only)
![](images/nf-cellsnplite-dark.png#gh-dark-mode-only)  
This pipeline produces 10 outputs, each were detailed in their corresponding steps:
* ***[output 1]:*** A VCF file listing genotyped SNPs and aggregated AD & DP infomation (without GT) `cellSNP.base.vcf.gz`
* ***[output 2]:*** A SNP x cell sparse matrix in “Matrix Market exchange formats”, containing the allele depths of the alternative (ALT) alleles, the sum of allele depths of the reference and alternative alleles (REF+ALT) and the sum of allele depths of all the alleles other than REF and ALT.  `cellSNP.tag.(AD/DP/OTH).mtx`
* ***[output 3]:*** A VCF file listing genotyped SNPs and aggregated AD & DP infomation (without GT) as a final output `cellSNP.base.vcf.gz`
* ***[output 4]:*** A VCF file listing genotyped SNPs and aggregated AD & DP infomation per cell `cellSNP.cells.vcf.gz`
* ***[output 5]:*** A SNP x cell sparse matrix in “Matrix Market exchange formats”, containing the allele depths of the alternative (ALT) alleles, the sum of allele depths of the reference and alternative alleles (REF+ALT) and the sum of allele depths of all the alleles other than REF and ALT as the final outputs `cellSNP.tag.(AD/DP/OTH).mtx`

![](images/workflow_modes.png)  
The default version of the pipeline runs 2-step approach as shown the diagram above. The other workflow (1-step approach) is more useful for smaller datasets since it involves both joint calling and genotyping, and therefore is substantially slower than 2-step approach.
* `2step`: runs 2b-1a steps 
* `1step`: runs 2a step

<details>

<summary>Workflow: 2step/1step</summary>

```
# to run 2step approach
nextflow run cellgeni/nf-cellsnplite -r main \
  -entry 2step \            # to choose run mode
  --SAMPLEFILE /path/to/sample/file \
  --project_tag test1 \   # to specify the run to add to the end of output folder (e.g. cellsnplite-results-test1)
  --ansi-log false \
  -resume
```
</details>



### 1. `2b`  

This step genotypes the 10x sample in a pseudo-bulk manner without given SNPs.

This step requires a single input:
  * BAM file from Cell Ranger or STARsolo


This step produces:  
* * ***[output 1]:*** A VCF file listing genotyped SNPs and aggregated AD & DP infomation (without GT) `cellSNP.base.vcf.gz`
* ***[output 2]:*** A SNP x cell sparse matrix in “Matrix Market exchange formats”, containing the allele depths of the alternative (ALT) alleles, the sum of allele depths of the reference and alternative alleles (REF+ALT) and the sum of allele depths of all the alleles other than REF and ALT.  `cellSNP.tag.(AD/DP/OTH).mtx`

### 2. `1a`

This step genotypes single cells at a list of given SNPs (called heterouzygous variants in previous step 2b).


This step produces: 
* ***[output 3]:*** A VCF file listing genotyped SNPs and aggregated AD & DP infomation (without GT) as a final output `cellSNP.base.vcf.gz`
* ***[output 4]:*** A VCF file listing genotyped SNPs and aggregated AD & DP infomation per cell `cellSNP.cells.vcf.gz`
* ***[output 5]:*** A SNP x cell sparse matrix in “Matrix Market exchange formats”, containing the allele depths of the alternative (ALT) alleles, the sum of allele depths of the reference and alternative alleles (REF+ALT) and the sum of allele depths of all the alleles other than REF and ALT as the final outputs `cellSNP.tag.(AD/DP/OTH).mtx`

### 3. `2a`  

This step genotypes the single cells without given SNPs. It does joint calling and genotyping, therefore it is substantially slower than 2-step approach, and it is only useful for smaller datasets, or for small chromosomes (mitochondrial).


This step produces:  
* ***[output 3]:*** A VCF file listing genotyped SNPs and aggregated AD & DP infomation (without GT) as a final output `cellSNP.base.vcf.gz`
* ***[output 4]:*** A VCF file listing genotyped SNPs and aggregated AD & DP infomation per cell `cellSNP.cells.vcf.gz`
* ***[output 5]:*** A SNP x cell sparse matrix in “Matrix Market exchange formats”, containing the allele depths of the alternative (ALT) alleles, the sum of allele depths of the reference and alternative alleles (REF+ALT) and the sum of allele depths of all the alleles other than REF and ALT as the final outputs `cellSNP.tag.(AD/DP/OTH).mtx`


# Changelog

For a version history/changelog, please see the [CHANGELOG file](CHANGELOG.md).