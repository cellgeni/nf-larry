# nf-larry

**nf-larry**: A scalable and reproducible workflow for lineage tracing with **LARRY** (Lineage And RNA Recovery) barcodes in single-cell RNA-seq.

LARRY delivers heritable DNA barcodes via lentiviral libraries; because the barcode is transcribed and captured in scRNA-seq, it enables clonal lineage tracing while measuring gene expression.  

**nf-larry** is a Nextflow pipeline for the LARRY analysis workflow, leveraging the private `CBUtools` Python package for barcode analysis in single-cell sequencing data.

LARRY technique was introduced in Weinreb et al 2020 Science paper: https://www.science.org/doi/10.1126/science.aaw3381

## Overview

LARRY analysis identifies and tracks clonal barcodes in single-cell data. The pipeline processes raw sequencing reads, applies quality control, assigns clones, and integrates clone information with gene expression (GEX) data.

**CBUtools** provides two key data structures:
- **CBUSeries**: A Pandas Series with MultiIndex (Cell Barcode, LARRY Barcode, UMI), counting the frequency of each group.
- **CBSeries**: A Pandas Series with MultiIndex (Cell Barcode, LARRY Barcode), counting the number of UMIs per group.

## Key Features

- **Scalable pipeline**: Automates LARRY analysis end-to-end.
- **Automatic FASTQ merging** when samples are split across lanes.
- **Built-in QC**: Generates a concise one-page PDF with plots and metrics.
- **Clone assignment per cell**, including multi-barcode cells.
- **Direct integration** with STARsolo or Cell Ranger outputs.
- **Reproducible setup** via Dockerfile, powered by the `CBUtools` Python library.

## How to run:

The recommended way to use nextflow is to run it in a screen session. These steps can be directly used in Sanger's FARM, but you can modify each step according to the environment you're working on or the job scheduler your HPC uses:

1. Start a screen session: `screen -S nf_run1`
2. Start a small interactive job for nextflow: `bsub -G cellgeni -n1 -R"span[hosts=1]" -Is -q long -R"select[mem>2000] rusage[mem=2000]" -M2000 bash`
3. Modify one of RESUME scripts in examples folder (pre-made Nextflow run scripts)
4. Run the RESUME scripts you modified: `./RESUME-larry-all`
5. You can leave your screen session and let it run in the background: `Ctrl+A, D`
6. You can go back to the screen session by: `screen -x nf_run1`

## Files:

* `main.nf` - the main Nextflow script that runs the nf-larry workflow.
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

* **Sample map (CSV):** The path to a CSV file mapping sample IDs from the LARRY libraries to their corresponding GEX sample IDs, along with a **group ID**.  
  - The group ID defines which LARRY samples should be treated as belonging to the same biological or technical group (e.g., technical replicates of the same culture or multiple lanes from one experiment).  
  - During the `ASSIGN_BARCODES` step, samples with the same group ID are combined before clone calling to ensure consistent assignment.  
* **FASTQ directory:** The path for the directory containing the paired-end FASTQ files (R1, R2) from the LARRY library preparation for each sample.
* **GEX outputs:** The path for the corresponding Gene Expression data. This can be either a directory Cell Ranger or STARsolo output folders.

### Outputs

The pipeline generates the following outputs, organized by sample and step:

0. **Concatenated FASTQs:** Contains all the reads in R1 in a single R1/R2 FASTQ files. (Only runs if a FASTQ file is coming from different lanes for the same sample)
1.  **Raw Barcode Reads (CBUSeries PKL):** Contains all sequences matching the LARRY barcode pattern extracted from FASTQ files, grouped by Cell Barcode, LARRY Barcode, and UMI (`find_larry_seqs.py`).
2.  **Filtered Barcodes (CBSeries PKL):** High-confidence LARRY barcodes remaining after applying read count, Hamming distance, and UMI filtering (`larry_qc.py`).
3.  **QC Report (PDF):** A detailed report with plots and metrics illustrating the effects of each QC filtering step (`larry_qc.py`).
4.  **Clone Assignments (PKL):** A cell-level summary dataframe mapping each cell barcode to its assigned clone identity (`assign_clones.py`).
5.  **Annotated GEX Data (H5AD):** Gene expression data integrated with the corresponding LARRY clone assignments for each cell (`match_gex.py`).
6.  **Clone Assignment Table (TSV/CSV):** A simple tabular file listing cell barcodes and their assigned clone identities (`match_gex.py`).
7.  **Clone Size Plot (PNG):** A plot visualizing the cumulative frequency distribution of clone sizes within the analyzed sample(s) (`match_gex.py`).

## Pipeline Steps

### 1. `CONCAT_FASTQS`

This initial step collapses lane-split FASTQs into one R1/R2 file per sample. It performs lossless concatenation and standardizes file naming for clean downstream matching.  

This step runs only for the samples if its FASTQ files are coming from multiple lanes, skips if a sample already has a single pair of R1/R2 FASTQ files.

It uses `cat` command for concatenation of FASTQ files. 

### 2. `FIND_LARRY_SEQS`

This initial step takes paired-end FASTQ files (R1, R2) from LARRY library for each sample as input.

It uses `cutadapt` to identify and extract LARRY barcode sequences from the raw sequencing reads based on a defined LARRY pattern. It then aggregates these sequences, grouping them by Cell Barcode (CBC), LARRY Barcode, and Unique Molecular Identifier (UMI), counting the reads for each unique combination. This step is parallelised if runs in a machine with multiple CPU cores.

This extracts candidate LARRY barcodes and produces a raw evidence table (CBUSeries) used for QC.  

The output is a CBUSeries PKL file for each sample (Output 1), storing a Pandas Series with a MultiIndex (CBC, Barcode, UMI) where the values represent the read counts, capturing all potential LARRY reads before filtering.

### 3. `LARRY_QC`

This step denoises barcodes, applies filtering, and generates a transparent one-page QC PDF. Using the CBUSeries PKL file generated by `FIND_LARRY_SEQS`, this step performs crucial quality control on the extracted LARRY barcodes. It applies a series of filters, including minimum read count thresholds, Hamming distance-based correction/filtering for sequencing errors, and UMI count thresholds for barcode reliability. After filtering, the data is collapsed into a CBSeries (grouping by CBC and Barcode, summing UMIs).

The process generates two outputs: a filtered CBSeries PKL file for each sample containing high-confidence LARRY barcodes (Output 2), and a comprehensive PDF QC report detailing the filtering process and results with plots and metrics (Output 3).

#### Default QC Thresholds
- Reads per CBC-barcode ≥ 8
- Hamming distance (edit-distance) ≤ 3
- UMIs per CBC-barcode ≥ 3
- Dominance filter = 0.1

### 4. `ASSIGN_BARCODES`

This step combines samples by group, assigns barcode(s) per cell from retained LARRY barcodes, and handles multi-barcode cells with an optional dominance filter.  This step takes the filtered CBSeries PKL file from `LARRY_QC` as input. Based on the QC'd LARRY barcodes present in each cell (CBC), it assigns a unique clone identity to each cell. Specific logic is employed to handle cells with multiple barcodes or ambiguous assignments.

The output is a PKL file (potentially aggregated for multiple samples) containing a cell-level summary dataframe that maps each cell barcode to its assigned LARRY clone identity (Output 4).

### 5. `MATCH_GEX`

This final step annotates GEX matrices with clone information, producing analysis-ready objects. This step requires three inputs: the cell-level summary PKL file from `ASSIGN_BARCODES`, the Gene Expression (GEX) data path (Cell Ranger or STARsolo outputs), and a sample map CSV linking LARRY sample IDs, GEX sample IDs, and group IDs. It integrates the LARRY-derived clone assignments with the GEX data by matching cells based on their barcodes. The clone information is added to the cell metadata of GEX data.

The outputs include annotated H5AD object(s) with integrated clone information (Output 5), a TSV/CSV file mapping cell barcodes to clone assignments (Output 6), and a PNG plot visualizing the cumulative frequency of clone sizes (Output 7).

## Future Work
- Add user-settable options for QC (read count, edit-distance, UMI thresholds, dominance filter).
- Migrate to nf-core to align the pipeline with nf-core standards for community support, modularity and long-term sustainability.
- Expand functionality with modules for state-fate clone matching, enabling direct linkage of clonal identity to cell-state transitions and outcomes.


## Changelog

See [CHANGELOG.md](CHANGELOG.md) for version history.
