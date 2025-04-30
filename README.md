# nf-larry

**nf-larry** is a Nextflow pipeline for the LARRY analysis workflow, leveraging the private `CBUtools` Python package for barcode analysis in single-cell sequencing data.

LARRY technique was introduced in Weinreb et al 2020 Science paper: https://www.science.org/doi/10.1126/science.aaw3381

## Overview

LARRY analysis identifies and tracks clonal barcodes in single-cell data. The pipeline processes raw sequencing reads, applies quality control, assigns clones, and integrates clone information with gene expression (GEX) data.

**CBUtools** provides two key data structures:
- **CBUSeries**: A Pandas Series with MultiIndex (Cell Barcode, LARRY Barcode, UMI), counting the frequency of each group.
- **CBSeries**: A Pandas Series with MultiIndex (Cell Barcode, LARRY Barcode), counting the number of UMIs per group.

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
  * `find_larry_seqs.py` - Extracts LARRY barcode sequences from sequencing reads and outputs a CBUSeries pickle file containing all matching reads.
  * `larry_qc.py` - Runs the LARRY QC pipeline: applies read count, Hamming distance, and UMI filtering to a CBUSeries, generates a PDF QC report, and outputs a filtered CBSeries pickle.
  * `assign_clones.py` - Assigns clone identities to cells based on filtered LARRY barcodes and outputs clone assignments for downstream analysis.
  * `match_gex` - Matches LARRY barcode-derived clone assignments to gene expression (GEX) data, integrating clone and transcriptome information for each cell.
* `RESUME-larry` -  a pre-made Nextflow run scripts for all steps
* `Dockerfile` - a dockerfile to reproduce the environment used to run the pipeline.

## Workflow

![](images/nf-larry-light.png#gh-light-mode-only)
![](images/nf-larry-dark.png#gh-dark-mode-only)

### Inputs

The pipeline requires the following primary inputs, typically configured in RESUME script:

*  The main directory containing the paired-end FASTQ files (R1, R2) from the LARRY library preparation for each sample.
*  The location of the corresponding Gene Expression data. This can be either a directory containing pre-processed H5AD files or STARsolo output directories.
*  A path to a JSON file that maps the sample identifiers used in the LARRY library to the corresponding sample identifiers in the GEX library.

### Outputs

The pipeline generates the following outputs, organized by sample and step:

1.  **Raw Barcode Reads (CBUSeries PKL):** Contains all sequences matching the LARRY barcode pattern extracted from FASTQ files, grouped by Cell Barcode, LARRY Barcode, and UMI (`find_larry_seqs.py`).
2.  **Filtered Barcodes (CBSeries PKL):** High-confidence LARRY barcodes remaining after applying read count, Hamming distance, and UMI filtering (`larry_qc.py`).
3.  **QC Report (PDF):** A detailed report with plots and metrics illustrating the effects of each QC filtering step (`larry_qc.py`).
4.  **Clone Assignments (PKL):** A cell-level summary dataframe mapping each cell barcode to its assigned clone identity (`assign_clones.py`).
5.  **Annotated GEX Data (H5AD):** Gene expression data integrated with the corresponding LARRY clone assignments for each cell (`match_gex.py`).
6.  **Clone Assignment Table (TSV/CSV):** A simple tabular file listing cell barcodes and their assigned clone identities (`match_gex.py`).
7.  **Clone Size Plot (PNG):** A plot visualizing the cumulative frequency distribution of clone sizes within the analyzed sample(s) (`match_gex.py`).

## Pipeline Steps

### 1. `FIND_LARRY_SEQS`

This initial step takes paired-end FASTQ files (R1, R2) from LARRY library for each sample as input.

It uses `cutadapt` to identify and extract LARRY barcode sequences from the raw sequencing reads based on a defined LARRY pattern. It then aggregates these sequences, grouping them by Cell Barcode (CBC), LARRY Barcode, and Unique Molecular Identifier (UMI), counting the reads for each unique combination.

The output is a CBUSeries PKL file for each sample (Output 1), storing a Pandas Series with a MultiIndex (CBC, Barcode, UMI) where the values represent the read counts, capturing all potential LARRY reads before filtering.

### 2. `LARRY_QC`

Using the CBUSeries PKL file generated by `FIND_LARRY_SEQS`, this step performs crucial quality control on the extracted LARRY barcodes. It applies a series of filters, including minimum read count thresholds, Hamming distance-based correction/filtering for sequencing errors, and UMI count thresholds for barcode reliability. After filtering, the data is collapsed into a CBSeries (grouping by CBC and Barcode, summing UMIs).

The process generates two outputs: a filtered CBSeries PKL file for each sample containing high-confidence LARRY barcodes (Output 2), and a comprehensive PDF QC report detailing the filtering process and results with plots and metrics (Output 3).

### 3. `ASSIGN_CLONES`

This step takes the filtered CBSeries PKL file from `LARRY_QC` as input. Based on the QC'd LARRY barcodes present in each cell (CBC), it assigns a unique clone identity to each cell. Specific logic is employed to handle cells with multiple barcodes or ambiguous assignments.

The output is a PKL file (potentially aggregated for multiple samples) containing a cell-level summary dataframe that maps each cell barcode to its assigned clone identity (Output 4).

### 4. `MATCH_GEX`

This final step requires three inputs: the cell-level summary PKL file from `ASSIGN_CLONES`, the Gene Expression (GEX) data path (H5AD files or STARsolo outputs), and a JSON file mapping LARRY sample IDs to GEX sample IDs. It integrates the LARRY-derived clone assignments with the GEX data by matching cells based on their barcodes. The clone information is added to the cell metadata of GEX data.

The outputs include annotated H5AD object(s) with integrated clone information (Output 5), a TSV/CSV file mapping cell barcodes to clone assignments (Output 6), and a PNG plot visualizing the cumulative frequency of clone sizes (Output 7).


## Changelog

See [CHANGELOG.md](CHANGELOG.md) for version history.