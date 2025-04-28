from .larry import process_larry
from scipy.stats.contingency import crosstab
import numpy as np
import pandas as pd
from . import hamming
from . import CBUclasses
from .CBUclasses import load_barcodes
import papermill as pm
from pathlib import Path
from subprocess import run

HERE = Path(__file__).parent


def get_barcodes(files, type="larry", valid_CBC=None, save="barcodes.csv"):
    """Retrieve barcodes from fastq

    Uses one a of the recipes to read barcodes from read1 and read2 files
    corresponding to a 10x Genomics experiment

    Parameters
    ----------
    files : dict
        Dictionary with paths to fastq files, needs to have 'r1' and 'r2' keys
    type : str
        recipe used to recover barcodes
    valid_CBC : Optional[list, np.ndarray]
        List or array with valid cell barcodes
    save : string
        Path where the barcode counts will be saved

    Returns
    --------
    CBUSeries
        Returns CBUSeries (derivation of pd.Series) with CBC, Barcode, UMI read counts

    """
    if type == "larry":
        counts = process_larry(files, valid_CBC=valid_CBC)
        return counts


def get_barcodes_report(
    files=None,
    report_template="larry1",
    save_path=".",
    prefix="larry1",
    make_html=True,
    include_code=True,
    min_reads=2,
    min_hamming=2,
    min_umis=2
):
    """Retrieve barcodes and generate report

    Retrieves the barcodes using the get_barcodes function and generates the
    report using the specified template.

    Parameters
    ----------
    files : dict
        Dictionary with paths to fastq files, needs to have 'r1' and 'r2' keys
    report_template : str
        Name of the template used for report generation (stored as .ipynb in the
        templates directory)
    save_path : PATH or str
        path where the report and barcodes should be saved
    prefix : str
        prefix for the report and barcodes files
    make_html : bool
        whether the .ipynb report should be automatically converted to html
    include_code : bool
        whether source code should be included in the report

    """
    if save_path is str:
        save_path = Path(save_path)
    ipynb_path = Path(save_path, prefix).with_suffix(".ipynb")
    csv_path = Path(save_path, prefix).with_suffix(".csv")

    if report_template == "larry1":
        template = str(HERE / "templates/larry1_template.ipynb")
        print("Converting notebook to HTML...")
        pm.execute_notebook(
            template,
            ipynb_path,
            parameters=dict(
                read1=files["r1"],
                read2=files["r2"],
                min_reads=min_reads,
                min_hamming=min_hamming,
                min_umis=min_umis,
                csv_path=str(csv_path),
            ),
        )

    print("Converting notebook to HTML...")
    if include_code:
        out = run(
            f"jupyter nbconvert --to html {ipynb_path}",
            capture_output=True,
            text=True,
            shell=True,
        )
    else:
        out = run(
            f"jupyter nbconvert --no-input --to html {ipynb_path}",
            capture_output=True,
            text=True,
            shell=True,
        )
    print(out.stderr)
    print(out.stdout)
