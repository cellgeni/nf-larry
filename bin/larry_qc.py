#!/usr/bin/env python
"""
larry_qc.py – single‑file QC pipeline + PDF report for LARRY libraries
=====================================================================

Usage (as in the Nextflow process)
----------------------------------
```
python larry_qc.py pkl_file=<raw_series.pkl> sample_id=SAMPLE42
```
Run `python larry_qc.py --help` to see all command‑line options.

* Expects **a `cbu.CBUSeries`** in the input pickle – aborts otherwise.
* Performs read‑count, Hamming‑distance and UMI filtering.
* Captures key metrics **before and after every step**.
* Computes the percentage of barcodes that appear in the official LARRY
  whitelist (optional) and shows the numbers in the PDF report.
* Produces:
  * `SAMPLEID_QC_report.pdf` (figures + metrics table)
  * `SAMPLEID_CBU_bar5.pkl`   (final QC‑passed `CBSeries`)
"""

# ────────────────────────────────────────────────────────────────────────
#  Standard library
import sys, pickle, logging, textwrap
from dataclasses import dataclass
from math import ceil
from typing import List

#  Third‑party libraries
import fire
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#  Project code (mounted in /opt inside the container)
sys.path.append("/opt")
import cbutools as cbu  # gives us CBUSeries / CBSeries

# ────────────────────────────────────────────────────────────────────────
#  📊  Plotting helpers

def _hist_counts(series, groupby, title, xlabel, *, bins: int = 50):
    """Return a Figure with a log-scaled histogram of grouped counts."""
    # use constrained layout to keep axes/labels visible
    fig, ax = plt.subplots(figsize=(4, 4), constrained_layout=True)

    grouped = series.groupby(groupby).sum()
    ax.hist(np.log10(grouped), bins=bins)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Frequency", labelpad=6)
    ax.set_yscale("log")

    vmax = grouped.max()
    ticks = ceil(np.log10(vmax)) + 1
    ax.set_xticks(range(ticks))
    ax.set_xticklabels([rf"$10^{i}$" if i else "1" for i in range(ticks)])
    ax.tick_params(axis="x", labelsize=8)
    ax.tick_params(axis="y", labelsize=8)
    ax.margins(x=0.02)
    return fig


def _barcode_elbow(series, title):
    """Return a Figure with the classic elbow plot (#barcodes per cell)."""
    fig, ax = plt.subplots(figsize=(4, 4), constrained_layout=True)
    ordered = series.sort_values(ascending=False)
    ax.plot(range(len(ordered)), ordered, marker=".")
    ax.set_title(title)
    ax.set_xlabel("Cells (sorted)")
    ax.set_ylabel("# barcodes per cell", labelpad=6)
    ax.tick_params(axis="x", labelsize=8)
    ax.tick_params(axis="y", labelsize=8)
    return fig

# ────────────────────────────────────────────────────────────────────────
#  📈  Metrics data‑class
@dataclass
class StepMetrics:
    name: str
    total_reads: int
    total_umis: int
    total_cbc: int
    total_barcode: int
    mean_bc_per_cell: float
    pct_whitelist: float | None = None   # filled later (optional)

# ────────────────────────────────────────────────────────────────────────
#  🖨️  Simple PDF helper
class PdfReport:
    """Collect matplotlib figures + text into one PDF."""

    def __init__(self, path: str):
        self._pdf = PdfPages(path)

    # —— high‑level helpers ———————————————————————————————
    def add_text(self, title: str, body: str):
        # A4 portrait with constrained layout so margins/labels are clean
        fig = plt.figure(figsize=(8.27, 11.69), constrained_layout=True)
        fig.text(0.5, 0.92, title, ha="center", va="top", size=18, weight="bold")

        # Render each line with spacing; preserves explicit \n in body
        y = 0.82
        for line in body.splitlines():
            if not line.strip():
                y -= 0.02  # small gap for blank lines
                continue
            fig.text(0.08, y, line, va="top", size=12)
            y -= 0.05

        self._pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    def add_fig(self, fig):
        self._pdf.savefig(fig)
        plt.close(fig)

    def add_metrics_table(self, metrics: List[StepMetrics]):
        fig, ax = plt.subplots(figsize=(8.27, 4))
        ax.axis("off")
        colnames = [
            "Step", "Reads", "UMIs", "Cells",
            "Barcodes", "Mean BC/cell", "% in whitelist",
        ]
        rows = [[
            m.name,
            f"{m.total_reads:,}",
            f"{m.total_umis:,}",
            m.total_cbc,
            m.total_barcode,
            f"{m.mean_bc_per_cell:.2f}",
            f"{m.pct_whitelist:.1f}%" if m.pct_whitelist is not None else "–",
        ] for m in metrics]
        table = ax.table(cellText=rows, colLabels=colnames, loc="center")
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1, 1.6)
        self._pdf.savefig(fig)
        plt.close(fig)

    def close(self):
        self._pdf.close()

# ────────────────────────────────────────────────────────────────────────
#  🔬  QC core – runs filter cascade & returns everything we need

def _run_filters(bar_m, *, min_reads: int, min_hamming: int, min_umis: int):
    """Apply all QC filters and collect figures + metrics."""
    figs: List[plt.Figure] = []
    metrics: List[StepMetrics] = []
    objs: List = []  # keep references for whitelist calc

    def _collect(name: str, obj):
        metrics.append(
            StepMetrics(
                name=name,
                total_reads=int(obj.sum()),
                total_umis=len(obj),
                total_cbc=len(obj.index.get_level_values("CBC").unique()),
                total_barcode=len(obj.index.get_level_values("Barcode").unique()),
                mean_bc_per_cell=obj.groupby("CBC").size().mean(),
            )
        )
        objs.append(obj)

    # —— step 0 : raw ——
    logging.info("Collecting raw metrics")
    _collect("raw", bar_m)
    figs.append(_hist_counts(bar_m, ["CBC", "Barcode"],
                             "Raw reads per CBC–Barcode", "Read counts (log10)"))

    # —— step 1 : minimum reads ——
    logging.info("Filtering by minimum reads")
    bar2 = bar_m.filter_by_reads(["CBC", "Barcode"], min_reads)
    _collect(f"reads ≥ {min_reads}", bar2)
    figs.append(_hist_counts(bar2, ["CBC", "Barcode"],
                             f"After reads ≥ {min_reads}", "Read counts (log10)"))
    figs.append(_barcode_elbow(bar2.groupby("CBC").size(),
                               "# barcodes per cell after reads filter"))

    # —— step 2 : hamming distance ——
    logging.info("Filtering by Hamming distance")
    prev_figs = set(plt.get_fignums())
    bar3 = bar2.filter_by_hamming(min_distance=min_hamming)
    if hasattr(bar3, "_hamming_fig"):
        figs.append(bar3._hamming_fig)
    _collect(f"hamming ≥ {min_hamming}", bar3)
    figs.append(_barcode_elbow(bar3.groupby("CBC").size(),
                               f"# barcodes per cell after hamming ≥ {min_hamming}"))

    # —— step 3 : collapse to UMI —— : collapse to UMI ——
    logging.info("Collecting collapse to UMI metrics")
    bar4 = bar3.count_UMI()  # CBSeries
    _collect("collapse to UMI", bar4)
    figs.append(_hist_counts(bar4, ["CBC", "Barcode"],
                             "UMIs per CBC–Barcode", "UMI counts (log10)"))

    # —— step 4 : minimum UMIs ——
    logging.info("Filtering by minimum UMIs")
    bar5 = bar4.filter_by_UMI(["CBC", "Barcode"], min_umis)
    _collect(f"UMIs ≥ {min_umis}", bar5)
    figs.append(_hist_counts(bar5, ["CBC", "Barcode"],
                             f"Final – UMIs ≥ {min_umis}", "UMI counts (log10)"))
    figs.append(_barcode_elbow(bar5.groupby("CBC").size(),
                               "# barcodes per cell – FINAL"))

    return bar5, metrics, figs, objs

# ────────────────────────────────────────────────────────────────────────
#  🚀  Public CLI

def larry_qc(
    pkl_file: str,
    sample_id: str,
    *,
    min_reads: int = 8,
    min_hamming: int = 3,
    min_umis: int = 3,
    whitelist_csv: str = "/opt/cbutools/larry_whitelist.csv",
    check_whitelist: bool = False,
    make_pdf: bool = False,
):
    """Run QC & generate PDF report for one library."""

    # —— set up logging ——
    logging.basicConfig(
        filename=f"{sample_id}.log",
        filemode="a",
        level=logging.INFO,
        format="%(asctime)s  %(levelname)s: %(message)s",
    )
    logging.info("Starting LARRY_QC")

    # —— load raw data ——
    with open(pkl_file, "rb") as handle:
        bar_m = pickle.load(handle)

    # —— run filter cascade ——
    bar5, metrics, figs, objs = _run_filters(bar_m,
                                             min_reads=min_reads,
                                             min_hamming=min_hamming,
                                             min_umis=min_umis)

    # —— whitelist percentages ——
    if check_whitelist:
        whitelist = pd.read_csv(whitelist_csv)["barcode"]

        def pct_in_whitelist(obj):
            uniq = obj.index.get_level_values("Barcode").unique()
            return 100 * uniq.isin(whitelist).sum() / len(uniq)

        for m, obj in zip(metrics, objs):
            m.pct_whitelist = pct_in_whitelist(obj)
        logging.info("Whitelist stats calculated (added to PDF only)")

    # —— PDF report ——
    if make_pdf:
        pdf_name = f"{sample_id}_QC_report.pdf"
        rpt = PdfReport(pdf_name)

        rpt.add_text(
            "LARRY QC report",
            "\n".join([
                f"Sample: {sample_id}",
                f"min reads ≥ {min_reads}",
                f"Hamming distance ≥ {min_hamming}",
                f"min UMIs ≥ {min_umis}",
                "Percentages refer to the fraction of unique barcodes that match the official LARRY whitelist."
            ]),
        )
        rpt.add_metrics_table(metrics)
        for fig in figs:
            rpt.add_fig(fig)
        rpt.close()
        logging.info(f"PDF report written → {pdf_name}")

    # —— save final object ——
    out_pkl = f"{sample_id}_CBU_bar5.pkl"
    with open(out_pkl, "wb") as handle:
        pickle.dump(bar5, handle)
    logging.info(f"Final CBSeries pickled → {out_pkl}")

    return out_pkl


if __name__ == "__main__":
    fire.Fire(larry_qc)
