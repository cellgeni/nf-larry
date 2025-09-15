# Changelog

## 25-258 (v0.2.0)

This release introduces significant changes to the pipeline.  

* New step: CONCAT_FASTQS
#### ➤ Improvements
- Input style has been modified.

#### ➤ Major Workflow & Script Changes
**FIND_LARRY_SEQS**
- Minor change on underlying CBUtools package and the Dockerfile allows parallelisation and significantly speeds up this process (more than 4x)

**LARRY_QC**
- Added more logging before each filter
- Fixed QC plots and text in PDF report
- Changed whitelist directory to the repository directory

**ASSIGN_CLONES → ASSIGN_BARCODES**
- Changed the step name for clarity.
- Refactored to support new input structure

**MATCH_GEX**
- Added support for Cell Ranger outputs.
- Refactored to support new input structure
