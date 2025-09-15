# Changelog

## 25-258 (v0.2.0)

**First public release.**
This version introduces major structural updates, improved performance, and clearer input handling.

### 🚀 New

* **New step: `CONCAT_FASTQS`**
  Collapses lane-split FASTQs into single R1/R2 files per sample.

  * Skips samples that already have a single pair.
  * Standardizes naming for clean downstream matching.

### ✨ Improvements

* **Inputs & Groups**

  * Restructured input style: added clearer **group ID** handling in sample maps.
  * Whitelist CSV is now stored directly in the repository.

* **FIND\_LARRY\_SEQS**

  * Underlying `CBUtools` and Dockerfile changes enable **parallelisation**, yielding **>4× faster runtime**.

* **LARRY\_QC**

  * Added more logging between filters for transparency.
  * Fixed QC plots and text in the PDF report.
  * Whitelist directory now defaults to the repo path.

* **ASSIGN\_BARCODES** *(formerly `ASSIGN_CLONES`)*

  * Renamed for clarity and consistency with the poster/README.
  * Refactored to support the new input structure.
  * Removed `combine_samples` functionality (now handled via group IDs).

* **MATCH\_GEX**

  * Added **support for Cell Ranger outputs** (in addition to STARsolo).
  * Refactored to support the new input structure.

### 🛠 Internal

* Updated Dockerfile: streamlined how `CBUtools` is installed and used.
* Modified process resources and updated Singularity image.
* Cleaned up example files/scripts.
* Removed old git submodule.

---

## (v0.1.0) *(Internal only)*

Prototype version used internally before public release.

* Core workflow with steps: `FIND_LARRY_SEQS`, `LARRY_QC`, `ASSIGN_CLONES`, `MATCH_GEX`.
* Added early prototype workflow: `until_clones`.

