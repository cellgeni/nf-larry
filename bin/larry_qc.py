#!/usr/bin/env python
"""
larry_qc.py - single-file QC pipeline + HTML report for LARRY libraries
=======================================================================

Usage (as in the Nextflow process)
----------------------------------
```
python larry_qc.py pkl_file=<raw_series.pkl> sample_id=SAMPLE42
```
Run `python larry_qc.py --help` to see all command‑line options.

* Expects **a `cbu.CBUSeries`** in the input pickle - aborts otherwise.
* Performs read‑count, Hamming‑distance and UMI filtering.
* Captures key metrics **before and after every step**.
* Computes the percentage of barcodes that appear in the official LARRY
    whitelist and shows the numbers in the report.
* Produces:
    * `SAMPLEID_QC_report.html` (interactive HTML report)
    * `SAMPLEID_CBU_bar5.pkl`   (final QC‑passed `CBSeries`)
"""

# ────────────────────────────────────────────────────────────────────────
#  Standard library
import sys, pickle, logging, base64
from dataclasses import dataclass, field
from math import ceil
from typing import List, Dict, Any, Optional
from datetime import datetime
from html import escape
from io import BytesIO

#  Third‑party libraries
import fire
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import cbutools as cbu  # gives us CBUSeries / CBSeries
from cbutools import hamming as cbu_hamming


# ────────────────────────────────────────────────────────────────────────
#  Metrics & data classes

@dataclass
class StepMetrics:
    name: str
    total_reads: int
    total_umis: int
    total_cbc: int
    total_barcode: int
    mean_bc_per_cell: float
    pct_whitelist: float | None = None


@dataclass
class BarcodeStats:
    """Statistics about barcodes after final filtering."""
    barcode_cell_counts: pd.Series = field(default_factory=pd.Series)
    barcode_umi_counts: pd.Series = field(default_factory=pd.Series)
    total_unique_barcodes: int = 0
    total_cells: int = 0


# ────────────────────────────────────────────────────────────────────────
#  Matplotlib plotting functions

def _fig_to_base64(fig: Figure, *, dpi: int = 150) -> str:
    """Return a base64-encoded PNG rendering of a matplotlib Figure."""
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", facecolor='white')
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("ascii")


def _hist_counts(series, groupby, title, xlabel, *, bins: int = 50, color: str = "#3b82f6"):
    """Return a Figure with a log-scaled histogram of grouped counts."""
    fig, ax = plt.subplots(figsize=(5, 4), constrained_layout=True)
    
    grouped = series.groupby(groupby).sum()
    log_data = np.log10(grouped[grouped > 0])
    
    ax.hist(log_data, bins=bins, color=color, edgecolor='white', linewidth=0.5, alpha=0.85)
    ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel("Frequency", fontsize=10)
    ax.set_yscale("log")
    
    if len(grouped) > 0:
        vmax = grouped.max()
        ticks = ceil(np.log10(vmax)) + 1
        ax.set_xticks(range(ticks))
        ax.set_xticklabels([f"$10^{i}$" if i else "1" for i in range(ticks)])
    
    ax.tick_params(axis="both", labelsize=9)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    return fig


def _barcode_elbow(series, title, color: str = "#3b82f6"):
    """Return a Figure with the classic elbow plot (#barcodes per cell)."""
    fig, ax = plt.subplots(figsize=(5, 4), constrained_layout=True)
    
    ordered = series.sort_values(ascending=False).reset_index(drop=True)
    ax.plot(range(len(ordered)), ordered, marker='.', markersize=3, 
            color=color, linewidth=1.5, alpha=0.8)
    
    ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
    ax.set_xlabel("Cells (sorted by barcode count)", fontsize=10)
    ax.set_ylabel("# barcodes per cell", fontsize=10)
    ax.tick_params(axis="both", labelsize=9)
    ax.grid(axis='both', alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    return fig


def _hamming_histogram(distances: np.ndarray, min_distance: int, seq_length: int | None, skip_filtering: bool = False):
    """Build a histogram figure for pairwise Hamming distances."""
    fig, ax = plt.subplots(figsize=(5, 4), constrained_layout=True)

    if distances.size:
        max_dist = int(distances.max()) + 1 if distances.size else (seq_length or 30)
        bins = np.arange(-0.5, max_dist + 1.5, 1)
        
        # Color bars below threshold differently only if filtering is applied
        n, bin_edges, patches = ax.hist(distances, bins=bins, color="#3b82f6", 
                                         edgecolor='white', linewidth=0.5, alpha=0.85)
        
        # Color bars below threshold red only if filtering is applied
        if not skip_filtering:
            for i, patch in enumerate(patches):
                if i < min_distance:
                    patch.set_facecolor('#ef4444')
    else:
        ax.text(0.5, 0.5, "Not enough barcodes\nto compute distances",
                ha="center", va="center", fontsize=11, transform=ax.transAxes, color="gray")

    ax.set_title("Pairwise Hamming Distance Distribution", fontsize=12, fontweight='bold', pad=10)
    ax.set_xlabel("Hamming distance", fontsize=10)
    ax.set_ylabel("Barcode pairs", fontsize=10)
    
    if seq_length:
        ax.set_xlim(-0.5, min(seq_length + 0.5, 30))
    
    # Only show threshold line if filtering is applied
    if not skip_filtering:
        ax.axvline(min_distance - 0.5, color="#f59e0b", linestyle="--", linewidth=2, label=f"Threshold = {min_distance}")
        ax.legend(loc='upper right', fontsize=9)
    
    ax.tick_params(axis="both", labelsize=9)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    return fig


def _barcode_cell_bar_chart(barcode_stats: BarcodeStats, top_n: int = 30):
    """Generate bar chart showing top barcodes by cell count."""
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    
    if barcode_stats.barcode_cell_counts.empty:
        ax.text(0.5, 0.5, "No barcode data available",
                ha="center", va="center", fontsize=11, transform=ax.transAxes, color="gray")
        return fig
    
    top_barcodes = barcode_stats.barcode_cell_counts.nlargest(min(top_n, len(barcode_stats.barcode_cell_counts)))
    
    # Truncate labels
    labels = [bc[:12] + "..." if len(str(bc)) > 12 else str(bc) for bc in top_barcodes.index]
    
    bars = ax.bar(range(len(top_barcodes)), top_barcodes.values, color="#10b981", 
                  edgecolor='#047857', linewidth=0.5, alpha=0.85)
    
    ax.set_title(f"Top {len(top_barcodes)} Barcodes by Cell Count", fontsize=12, fontweight='bold', pad=10)
    ax.set_xlabel("Barcode", fontsize=10)
    ax.set_ylabel("Number of cells", fontsize=10)
    ax.set_xticks(range(len(top_barcodes)))
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    ax.tick_params(axis="y", labelsize=9)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    return fig


def _barcode_umi_scatter(barcode_stats: BarcodeStats):
    """Generate scatter plot of UMI counts vs cell count for each barcode."""
    fig, ax = plt.subplots(figsize=(5, 4), constrained_layout=True)
    
    if barcode_stats.barcode_cell_counts.empty:
        ax.text(0.5, 0.5, "No barcode data available",
                ha="center", va="center", fontsize=11, transform=ax.transAxes, color="gray")
        return fig
    
    cell_counts = barcode_stats.barcode_cell_counts
    umi_counts = barcode_stats.barcode_umi_counts
    
    common_barcodes = cell_counts.index.intersection(umi_counts.index)
    x_data = cell_counts.loc[common_barcodes].values
    y_data = umi_counts.loc[common_barcodes].values
    
    scatter = ax.scatter(x_data, y_data, c=y_data, cmap='viridis', 
                         s=30, alpha=0.7, edgecolors='none')
    
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
    cbar.set_label('Total UMIs', fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    
    ax.set_title("Barcode Quality: Cells vs UMIs", fontsize=12, fontweight='bold', pad=10)
    ax.set_xlabel("Number of cells", fontsize=10)
    ax.set_ylabel("Total UMIs", fontsize=10)
    ax.tick_params(axis="both", labelsize=9)
    ax.grid(axis='both', alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    return fig


# ────────────────────────────────────────────────────────────────────────
#  Hamming filter

def _detect_hamming_threshold(distances: np.ndarray, default_threshold: int = 3) -> int:
    """
    Detect optimal Hamming distance threshold from bimodal distribution.
    Finds the valley (local minimum) between two peaks.
    
    Args:
        distances: Array of pairwise Hamming distances
        default_threshold: Fallback value if detection fails
    
    Returns:
        Detected threshold (valley between peaks) or default_threshold
    """
    if distances.size < 10:
        logging.info("Not enough data for auto threshold detection, using default")
        return default_threshold
    
    try:
        from scipy.signal import find_peaks
        from scipy.ndimage import gaussian_filter1d
        
        # Create histogram
        max_dist = int(distances.max()) + 1
        hist, bin_edges = np.histogram(distances, bins=np.arange(-0.5, max_dist + 1.5, 1))
        
        # Smooth histogram to reduce noise
        smoothed = gaussian_filter1d(hist.astype(float), sigma=1.0)
        
        # Find peaks
        peaks, properties = find_peaks(smoothed, prominence=smoothed.max() * 0.1)
        
        if len(peaks) < 2:
            logging.info("Less than 2 peaks detected, using default threshold")
            return default_threshold
        
        # Get the two highest peaks
        peak_heights = smoothed[peaks]
        top_two_indices = np.argsort(peak_heights)[-2:]
        top_two_peaks = peaks[top_two_indices]
        top_two_peaks = np.sort(top_two_peaks)  # Sort by position
        
        # Find minimum between the two peaks
        valley_region = smoothed[top_two_peaks[0]:top_two_peaks[1]+1]
        if len(valley_region) == 0:
            return default_threshold
        
        valley_idx = np.argmin(valley_region)
        threshold = top_two_peaks[0] + valley_idx
        
        # Ensure threshold is reasonable (between 1 and max_dist)
        threshold = max(1, min(threshold, max_dist - 1))
        
        logging.info(f"Auto-detected Hamming threshold: {threshold} (peaks at {top_two_peaks[0]} and {top_two_peaks[1]})")
        return int(threshold)
        
    except ImportError:
        logging.warning("scipy not available for auto threshold detection, using simpler method")
        # Fallback: simple method without scipy
        max_dist = int(distances.max()) + 1
        hist, _ = np.histogram(distances, bins=np.arange(-0.5, max_dist + 1.5, 1))
        
        # Simple smoothing with moving average
        if len(hist) < 5:
            return default_threshold
        
        smoothed = np.convolve(hist, np.ones(3)/3, mode='same')
        
        # Find local maxima (simple approach)
        peaks = []
        for i in range(1, len(smoothed) - 1):
            if smoothed[i] > smoothed[i-1] and smoothed[i] > smoothed[i+1] and smoothed[i] > hist.max() * 0.1:
                peaks.append(i)
        
        if len(peaks) < 2:
            logging.info("Less than 2 peaks detected (fallback method), using default threshold")
            return default_threshold
        
        # Find valley between first two peaks
        peak1, peak2 = sorted(peaks[:2])
        valley_idx = peak1 + np.argmin(smoothed[peak1:peak2+1])
        threshold = max(1, min(valley_idx, max_dist - 1))
        
        logging.info(f"Auto-detected Hamming threshold (fallback): {threshold}")
        return int(threshold)
        
    except Exception as e:
        logging.warning(f"Error in auto threshold detection: {e}, using default")
        return default_threshold


def _calculate_hamming_distances(bar_series):
    """Calculate Hamming distances without filtering."""
    counts = bar_series.groupby(["Barcode"]).sum()
    if counts.empty:
        logging.info("No barcodes available for Hamming distance calculation.")
        return bar_series, np.array([]), None, []

    seqs = counts.index.astype(str).to_numpy()
    logging.info("Calculating pairwise Hamming distances for %d barcodes", len(seqs))

    # Compute Hamming distance matrix
    hdist = cbu_hamming.compute_hamming_matrix(seqs)
    
    # Get sequence length from the actual sequences
    seq_length = cbu_hamming.get_sequence_length(seqs)
    
    # Extract upper triangle for histogram
    condensed = (
        hdist[np.triu_indices_from(hdist, k=1)] if len(seqs) > 1 else np.array([])
    )
    
    return bar_series, condensed, seq_length, counts.index.tolist()


def _apply_hamming_filter(bar_series, *, min_distance: int):
    """Apply Hamming-distance filtering and return filtered series + distance stats."""
    counts = bar_series.groupby(["Barcode"]).sum()
    if counts.empty:
        logging.info("No barcodes available for Hamming filter; skipping distance plot.")
        return bar_series, np.array([]), None, []

    seqs = counts.index.astype(str).to_numpy()
    logging.info("Calculating pairwise Hamming distances for %d barcodes", len(seqs))

    # Compute Hamming distance matrix
    hdist = cbu_hamming.compute_hamming_matrix(seqs)
    
    # Get sequence length from the actual sequences
    seq_length = cbu_hamming.get_sequence_length(seqs)
    
    # Extract upper triangle for histogram
    condensed = (
        hdist[np.triu_indices_from(hdist, k=1)] if len(seqs) > 1 else np.array([])
    )

    toreject: set[str] = set()
    ties: dict[str, List[str]] = {}

    for row_idx, seq in enumerate(seqs):
        row = hdist[row_idx]
        mask = row < min_distance
        mask[row_idx] = False
        if not mask.any():
            continue

        neighbours = np.where(mask)[0]
        neighbour_counts = counts.iloc[neighbours]
        seq_count = counts.iloc[row_idx]

        if np.any(seq_count < neighbour_counts):
            toreject.add(seq)
            continue

        max_neighbour = neighbour_counts.max()
        if seq_count == max_neighbour:
            tied = neighbour_counts.index[neighbour_counts == max_neighbour].tolist()
            if tied:
                ties[seq] = tied

    tokeep = counts.index[~counts.index.isin(toreject)]
    logging.info(
        "Hamming filter retained %d barcodes and rejected %d (ties: %d)",
        len(tokeep), len(toreject), len(ties),
    )
    if ties:
        logging.warning("Hamming filter ties detected for %d barcodes.", len(ties))

    filtered = bar_series.loc[:, tokeep, :]
    return filtered, condensed, seq_length, tokeep


# ────────────────────────────────────────────────────────────────────────
#  HTML Report Generator

def _write_html_report(
    path: str,
    *,
    sample_id: str,
    metrics: List[StepMetrics],
    figures: Dict[str, str],  # name -> base64 encoded PNG
    filter_params: dict,
    barcode_stats: BarcodeStats,
    barcode_table_data: List[Dict],
):
    """Render an HTML QC report with embedded matplotlib figures."""
    generated_ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    if not metrics:
        raise ValueError("No metrics were collected; HTML report cannot be rendered.")

    def safe_pct(part, whole):
        return 0.0 if not whole else 100 * float(part) / float(whole)

    raw_metrics = metrics[0]
    final_metrics = metrics[-1]
    
    # Build summary stats
    reads_retained_pct = safe_pct(final_metrics.total_reads, raw_metrics.total_reads)
    umi_retained_pct = safe_pct(final_metrics.total_umis, raw_metrics.total_umis)
    cbc_retained_pct = safe_pct(final_metrics.total_cbc, raw_metrics.total_cbc)
    barcode_retained_pct = safe_pct(final_metrics.total_barcode, raw_metrics.total_barcode)

    safe_sample = escape(sample_id)
    
    # Build metrics table rows
    metrics_rows = ""
    for idx, m in enumerate(metrics):
        reads_vs_raw = safe_pct(m.total_reads, raw_metrics.total_reads)
        whitelist_pct = f"{m.pct_whitelist:.1f}" if m.pct_whitelist is not None else "-"
        row_class = "highlight" if idx == len(metrics) - 1 else ""
        metrics_rows += f"""
            <tr class="{row_class}">
                <td class="step-name">{escape(m.name)}</td>
                <td class="num">{m.total_reads:,}</td>
                <td class="num pct">{reads_vs_raw:.1f}</td>
                <td class="num">{m.total_umis:,}</td>
                <td class="num">{m.total_cbc:,}</td>
                <td class="num">{m.total_barcode:,}</td>
                <td class="num">{m.mean_bc_per_cell:.2f}</td>
                <td class="num">{whitelist_pct}</td>
            </tr>"""

    # Build barcode table (top 100)
    barcode_rows = ""
    for i, bc in enumerate(barcode_table_data[:100], 1):
        wl_mark = "Y" if bc['in_whitelist'] else "-"
        barcode_rows += f"""
            <tr>
                <td class="num">{i}</td>
                <td class="barcode-seq">{escape(str(bc['barcode']))}</td>
                <td class="num">{bc['cells']:,}</td>
                <td class="num">{bc['umis']:,}</td>
                <td class="num">{bc['mean_umis_per_cell']:.2f}</td>
                <td class="whitelist">{wl_mark}</td>
            </tr>"""

    # Build figure HTML
    def fig_html(name: str, caption: str) -> str:
        if name in figures:
            return f'''<figure>
                <img src="data:image/png;base64,{figures[name]}" alt="{escape(caption)}">
                <figcaption>{escape(caption)}</figcaption>
            </figure>'''
        return ""

    # Whitelist info
    whitelist_summary = (
        f"{final_metrics.pct_whitelist:.1f}% of unique barcodes in whitelist"
        if final_metrics.pct_whitelist is not None
        else "Whitelist check disabled"
    )
    whitelist_info = ""  # Not used in header anymore

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>LARRY QC Report - {safe_sample}</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
:root {{
    color-scheme: light dark;
    --accent: #1f77b4;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
}}
* {{ box-sizing: border-box; }}
body {{
    margin: 0;
    background: #f7f9fc;
    color: #1a1c1f;
}}
main {{
    max-width: 1180px;
    margin: 0 auto;
    padding: 2rem 1.5rem 4rem;
}}
header {{
    border-bottom: 3px solid var(--accent);
    margin-bottom: 2rem;
    padding-bottom: 1rem;
}}
header h1 {{
    margin: 0 0 0.5rem;
    font-size: clamp(2rem, 4vw, 2.75rem);
}}
h2 {{
    font-size: 1.35rem;
    margin: 0 0 1rem;
}}
.summary {{
    display: grid;
    gap: 0.25rem;
    font-size: 1rem;
}}
.toc {{
    display: flex;
    gap: 0.75rem;
    flex-wrap: wrap;
    margin-top: 1.25rem;
}}
.toc button {{
    border: 1px solid rgba(15, 23, 42, 0.15);
    background: rgba(15, 23, 42, 0.05);
    padding: 0.5rem 0.9rem;
    border-radius: 999px;
    cursor: pointer;
    font-weight: 600;
}}
.toc button:hover {{
    border-color: var(--accent);
    color: var(--accent);
}}
section {{
    margin-bottom: 2.5rem;
}}
.summary-grid {{
    display: grid;
    gap: 1rem;
    grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
}}
.summary-card {{
    background: #fff;
    border-radius: 0.5rem;
    padding: 1.25rem 1.5rem;
    box-shadow: 0 1px 3px rgba(15, 23, 42, 0.1);
}}
.summary-card h3 {{
    margin: 0;
    font-size: 0.85rem;
    color: rgba(15, 23, 42, 0.6);
    font-weight: 500;
}}
.summary-card .metric {{
    font-size: 1.75rem;
    font-weight: 700;
    margin: 0.3rem 0 0.2rem;
    color: #1a1c1f;
}}
.summary-card .delta {{
    color: rgba(15, 23, 42, 0.55);
    margin: 0;
    font-size: 0.9rem;
}}
.summary-note {{
    margin-top: 1rem;
    color: rgba(15, 23, 42, 0.7);
}}
.filters {{
    background: rgba(31, 119, 180, 0.08);
    padding: 1.2rem 1.5rem;
    border-radius: 0.75rem;
    margin: 0 auto 2.5rem;
    max-width: 640px;
    width: 100%;
}}
.filters ul {{
    list-style: none;
    margin: 0;
    padding: 0;
}}
.filters li {{
    display: flex;
    justify-content: space-between;
    gap: 0.75rem;
    padding: 0.45rem 0;
    border-bottom: 1px solid rgba(15, 23, 42, 0.08);
}}
.filters li:last-child {{
    border-bottom: none;
}}
.filters .label {{
    font-weight: 600;
}}
.filters .value {{
    font-variant-numeric: tabular-nums;
}}
.table-container {{
    overflow-x: auto;
    border-radius: 0.75rem;
    box-shadow: 0 12px 30px rgba(15, 23, 42, 0.07);
}}
table {{
    width: 100%;
    border-collapse: collapse;
    background: #fff;
    font-size: 0.9rem;
}}
thead {{
    background: var(--accent);
    color: #f8fafc;
}}
th, td {{
    padding: 0.75rem 1rem;
    text-align: left;
}}
th {{
    font-weight: 600;
    font-size: 0.8rem;
    text-transform: uppercase;
    letter-spacing: 0.5px;
}}
tbody tr {{
    border-bottom: 1px solid rgba(15, 23, 42, 0.08);
    transition: background 0.2s;
}}
tbody tr:hover {{
    background: #f8fafc;
}}
tbody tr.highlight {{
    background: #dbeafe;
    font-weight: 600;
}}
td.num {{
    text-align: right;
    font-variant-numeric: tabular-nums;
}}
td.step-name {{
    font-weight: 500;
}}
td.barcode-seq {{
    font-family: "SF Mono", Monaco, "Cascadia Code", monospace;
    font-size: 0.85rem;
}}
td.whitelist {{
    text-align: center;
    color: var(--accent);
    font-weight: bold;
}}
.figures-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
    gap: 1.5rem;
}}
figure {{
    margin: 0;
    background: #fff;
    border-radius: 0.75rem;
    padding: 1rem 1.25rem 1.5rem;
    box-shadow: 0 12px 30px rgba(15, 23, 42, 0.08);
}}
figure img {{
    width: 100%;
    height: auto;
    border-radius: 0.5rem;
}}
figcaption {{
    margin-top: 0.75rem;
    font-size: 0.95rem;
    color: rgba(15, 23, 42, 0.7);
    text-align: center;
}}
.barcode-summary {{
    display: grid;
    gap: 1rem;
    grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
    margin-bottom: 1.5rem;
}}
.barcode-summary .summary-card {{
    background: #fff;
    border-radius: 0.5rem;
    padding: 1.25rem 1.5rem;
    box-shadow: 0 1px 3px rgba(15, 23, 42, 0.1);
}}
.barcode-summary .summary-card h3 {{
    margin: 0;
    font-size: 0.85rem;
    color: rgba(15, 23, 42, 0.6);
    font-weight: 500;
}}
.barcode-summary .summary-card .metric {{
    font-size: 1.75rem;
    font-weight: 700;
    margin: 0.3rem 0 0.2rem;
    color: #1a1c1f;
}}
.barcode-summary .summary-card .delta {{
    color: rgba(15, 23, 42, 0.55);
    margin: 0;
    font-size: 0.9rem;
}}
.card {{
    background: #fff;
    border-radius: 0.75rem;
    padding: 1.5rem;
    box-shadow: 0 12px 30px rgba(15, 23, 42, 0.07);
}}
.card h3 {{
    margin: 0 0 1rem;
    color: #1a1c1f;
}}
footer {{
    text-align: center;
    padding: 2rem;
    color: rgba(15, 23, 42, 0.6);
    font-size: 0.85rem;
    border-top: 1px solid rgba(15, 23, 42, 0.1);
    margin-top: 2rem;
}}
@media (max-width: 768px) {{
    .figures-grid {{ grid-template-columns: 1fr; }}
    .barcode-summary {{ flex-direction: column; gap: 1rem; }}
}}
@media (prefers-color-scheme: dark) {{
    body {{
        background: #0f172a;
        color: #e2e8f0;
    }}
    header {{
        border-bottom-color: rgba(125, 211, 252, 0.5);
    }}
    .filters {{
        background: rgba(56, 189, 248, 0.12);
        border: 1px solid rgba(148, 163, 184, 0.3);
    }}
    .filters li {{
        border-bottom-color: rgba(148, 163, 184, 0.2);
    }}
    table {{
        background: #16213f;
    }}
    thead {{
        background: rgba(56, 189, 248, 0.35);
    }}
    tbody tr {{
        border-bottom-color: rgba(148, 163, 184, 0.15);
    }}
    tbody tr:hover {{
        background: #334155;
    }}
    tbody tr.highlight {{
        background: #1e3a8a;
    }}
    .summary-card {{
        background: #1e293b;
        box-shadow: 0 1px 3px rgba(0, 0, 0, 0.2);
        border: 1px solid rgba(148, 163, 184, 0.15);
    }}
    .summary-card h3 {{
        color: rgba(226, 232, 240, 0.6);
    }}
    .summary-card .metric {{
        color: #e2e8f0;
    }}
    .summary-card .delta {{
        color: rgba(226, 232, 240, 0.55);
    }}
    .summary-note {{
        color: rgba(226, 232, 240, 0.7);
    }}
    figure {{
        background: #16213f;
        box-shadow: none;
        border: 1px solid rgba(148, 163, 184, 0.25);
    }}
    figcaption {{
        color: rgba(226, 232, 240, 0.7);
    }}
    .barcode-summary .summary-card {{
        background: #1e293b;
        box-shadow: 0 1px 3px rgba(0, 0, 0, 0.2);
        border: 1px solid rgba(148, 163, 184, 0.15);
    }}
    .barcode-summary .summary-card h3 {{
        color: rgba(226, 232, 240, 0.6);
    }}
    .barcode-summary .summary-card .metric {{
        color: #e2e8f0;
    }}
    .barcode-summary .summary-card .delta {{
        color: rgba(226, 232, 240, 0.55);
    }}
    .card {{
        background: #16213f;
        box-shadow: none;
        border: 1px solid rgba(148, 163, 184, 0.25);
    }}
    .card h3 {{
        color: #e2e8f0;
    }}
    .table-container {{
        box-shadow: none;
    }}
    td.whitelist {{
        color: #38bdf8;
    }}
    footer {{
        color: rgba(148, 163, 184, 0.7);
        border-top-color: rgba(148, 163, 184, 0.2);
    }}
    .toc button {{
        background: #0f1b3a;
        border-color: rgba(148, 163, 184, 0.4);
        color: rgba(226, 232, 240, 0.85);
    }}
    .toc button:hover {{
        color: #7dd3fc;
        border-color: #7dd3fc;
    }}
}}
html {{ scroll-behavior: smooth; }}
</style>
<script>
document.addEventListener("DOMContentLoaded", function() {{
    document.querySelectorAll("[data-scroll]").forEach(function(btn) {{
        btn.addEventListener("click", function() {{
            const target = document.getElementById(btn.dataset.scroll);
            if (target) {{
                target.scrollIntoView({{ behavior: "smooth", block: "start" }});
            }}
        }});
    }});
}});
</script>
</head>
<body>
<main>
<header>
    <h1>LARRY QC Report</h1>
    <div class="summary">
        <span><strong>Sample:</strong> {safe_sample}</span>
        <span><strong>Generated:</strong> {escape(generated_ts)}</span>
        {whitelist_info if final_metrics.pct_whitelist is not None else ""}
    </div>
    <nav class="toc">
        <button type="button" data-scroll="summary">Summary</button>
        <button type="button" data-scroll="pipeline">Pipeline</button>
        <button type="button" data-scroll="figures">Figures</button>
        <button type="button" data-scroll="barcodes">Barcodes</button>
    </nav>
</header>

    <section id="summary">
        <h2>Summary</h2>
        <div class="summary-grid">
            <div class="summary-card">
                <h3>Reads Retained</h3>
                <div class="metric">{final_metrics.total_reads:,}</div>
                <div class="delta">{reads_retained_pct:.1f}% of raw reads</div>
            </div>
            <div class="summary-card">
                <h3>UMIs Retained</h3>
                <div class="metric">{final_metrics.total_umis:,}</div>
                <div class="delta">{umi_retained_pct:.1f}% of raw UMIs</div>
            </div>
            <div class="summary-card">
                <h3>Cells Retained</h3>
                <div class="metric">{final_metrics.total_cbc:,}</div>
                <div class="delta">{cbc_retained_pct:.1f}% of raw cells</div>
            </div>
            <div class="summary-card">
                <h3>Unique Barcodes</h3>
                <div class="metric">{final_metrics.total_barcode:,}</div>
                <div class="delta">{barcode_retained_pct:.1f}% of raw barcodes</div>
            </div>
            <div class="summary-card">
                <h3>Mean BC per cell</h3>
                <div class="metric">{final_metrics.mean_bc_per_cell:.2f}</div>
                <div class="delta">Final step</div>
            </div>
        </div>
        <p class="summary-note">{escape(whitelist_summary)}</p>

    <h2>Filtering Parameters</h2>
    <ul>
        <li>
            <span class="label">Min reads per CBC-Barcode</span>
            <span class="value">&ge; {filter_params['min_reads']}</span>
        </li>
        <li>
            <span class="label">Min Hamming distance</span>
            <span class="value">{'SKIPPED' if filter_params['min_hamming_skipped'] else f"≥ {filter_params['min_hamming']}"}</span>
        </li>
        <li>
            <span class="label">Min UMIs per CBC-Barcode</span>
            <span class="value">&ge; {filter_params['min_umis']}</span>
        </li>
    </ul>
    </section>


    <section id="pipeline">
        <h2>QC Pipeline Steps</h2>
        <div class="table-container">
            <table>
                <thead>
                    <tr>
                        <th>Step</th>
                        <th class="num">Reads</th>
                        <th class="num">% Retained</th>
                        <th class="num">UMIs</th>
                        <th class="num">Cells</th>
                        <th class="num">Barcodes</th>
                        <th class="num">Mean BC/Cell</th>
                        <th class="num">% Whitelist</th>
                    </tr>
                </thead>
                <tbody>
                    {metrics_rows}
                </tbody>
            </table>
        </div>
    </section>

    <section id="figures">
        <h2>Figures</h2>
        <div class="figures-grid">
            {fig_html('raw_hist', 'Raw reads per CBC-Barcode')}
            {fig_html('reads_hist', f"After reads ≥ {filter_params['min_reads']}")}
            {fig_html('elbow_reads', 'Barcodes per cell (after reads filter)')}
            {fig_html('hamming', 'Pairwise Hamming distance distribution')}
            {fig_html('elbow_hamming', 'Barcodes per cell (Hamming filter SKIPPED)' if filter_params['min_hamming_skipped'] else f"Barcodes per cell (after Hamming ≥ {filter_params['min_hamming']})")}
            {fig_html('umi_hist', 'UMIs per CBC-Barcode')}
            {fig_html('final_hist', f"Final - UMIs ≥ {filter_params['min_umis']}")}
            {fig_html('elbow_final', 'Barcodes per cell (FINAL)')}
        </div>
    </section>

    <section id="barcodes">
        <h2>Barcode Analysis (Post-QC)</h2>
        
        <div class="barcode-summary">
            <div class="summary-card">
                <h3>Total cells in sample</h3>
                <div class="metric">{barcode_stats.total_cells:,}</div>
                <div class="delta">After QC filtering</div>
            </div>
            <div class="summary-card">
                <h3>Unique barcodes</h3>
                <div class="metric">{barcode_stats.total_unique_barcodes:,}</div>
                <div class="delta">Detected in sample</div>
            </div>
        </div>
        
        <div class="figures-grid" style="margin-bottom: 1.5rem;">
            {fig_html('bc_cells', 'Top barcodes by cell count')}
            {fig_html('bc_scatter', 'Barcode quality: Cells vs UMIs')}
        </div>
        
        <div class="card">
            <h3>Top Barcodes (by cell count)</h3>
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>#</th>
                            <th>Barcode</th>
                            <th>Cells</th>
                            <th>Total UMIs</th>
                            <th>Mean UMI/Cell</th>
                            <th>Whitelist</th>
                        </tr>
                    </thead>
                    <tbody>
                        {barcode_rows}
                    </tbody>
                </table>
            </div>
        </div>
    </section>

<footer>
    Generated by larry_qc.py | LARRY Barcode QC Pipeline
</footer>
</main>
</body>
</html>'''

    with open(path, "w", encoding="utf-8") as handle:
        handle.write(html)


# ────────────────────────────────────────────────────────────────────────
#  QC core - runs filter cascade & collects data

def _run_filters(bar_m, *, min_reads: int, min_hamming: int, min_umis: int, auto_hamming: bool = False, skip_hamming: bool = False):
    """Apply all QC filters and collect figures + metrics."""
    figures: Dict[str, str] = {}
    metrics: List[StepMetrics] = []
    objs: List = []

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

    # Step 0: Raw
    logging.info("Collecting raw metrics")
    _collect("raw", bar_m)
    figures["raw_hist"] = _fig_to_base64(
        _hist_counts(bar_m, ["CBC", "Barcode"], "Raw Reads per CBC-Barcode", "Read counts (log10)")
    )

    # Step 1: Minimum reads
    logging.info("Filtering by minimum reads")
    bar2 = bar_m.filter_by_reads(["CBC", "Barcode"], min_reads)
    _collect(f"reads ≥ {min_reads}", bar2)
    figures["reads_hist"] = _fig_to_base64(
        _hist_counts(bar2, ["CBC", "Barcode"], f"After Reads ≥ {min_reads}", "Read counts (log10)")
    )
    figures["elbow_reads"] = _fig_to_base64(
        _barcode_elbow(bar2.groupby("CBC").size(), "Barcodes per Cell (after reads filter)")
    )

    # Step 2: Hamming distance
    if skip_hamming:
        logging.info("Calculating Hamming distances (filtering SKIPPED)")
        bar3, hdist_values, seq_length, _ = _calculate_hamming_distances(bar2)
        actual_hamming = None  # No threshold when skipped
        figures["hamming"] = _fig_to_base64(
            _hamming_histogram(hdist_values, min_hamming, seq_length, skip_filtering=True)
        )
        _collect("hamming filter SKIPPED", bar3)
        figures["elbow_hamming"] = _fig_to_base64(
            _barcode_elbow(bar3.groupby("CBC").size(), "Barcodes per Cell (Hamming filter SKIPPED)")
        )
    else:
        logging.info("Filtering by Hamming distance")
        bar3, hdist_values, seq_length, _ = _apply_hamming_filter(bar2, min_distance=min_hamming)
        
        # Auto-detect threshold if requested
        actual_hamming = min_hamming
        if auto_hamming and hdist_values.size > 0:
            actual_hamming = _detect_hamming_threshold(hdist_values, default_threshold=min_hamming)
            logging.info(f"Re-applying Hamming filter with auto-detected threshold: {actual_hamming}")
            bar3, hdist_values, seq_length, _ = _apply_hamming_filter(bar2, min_distance=actual_hamming)
        
        figures["hamming"] = _fig_to_base64(
            _hamming_histogram(hdist_values, actual_hamming, seq_length)
        )
        _collect(f"hamming ≥ {actual_hamming}", bar3)
        figures["elbow_hamming"] = _fig_to_base64(
            _barcode_elbow(bar3.groupby("CBC").size(), f"Barcodes per Cell (after Hamming ≥ {actual_hamming})")
        )

    # Step 3: Collapse to UMI
    logging.info("Collapsing to UMI")
    bar4 = bar3.count_UMI()
    _collect("collapse to UMI", bar4)
    figures["umi_hist"] = _fig_to_base64(
        _hist_counts(bar4, ["CBC", "Barcode"], "UMIs per CBC-Barcode", "UMI counts (log10)", color="#10b981")
    )

    # Step 4: Minimum UMIs
    logging.info("Filtering by minimum UMIs")
    bar5 = bar4.filter_by_UMI(["CBC", "Barcode"], min_umis)
    _collect(f"UMIs ≥ {min_umis}", bar5)
    figures["final_hist"] = _fig_to_base64(
        _hist_counts(bar5, ["CBC", "Barcode"], f"Final - UMIs ≥ {min_umis}", "UMI counts (log10)", color="#8b5cf6")
    )
    figures["elbow_final"] = _fig_to_base64(
        _barcode_elbow(bar5.groupby("CBC").size(), "Barcodes per Cell - FINAL", color="#8b5cf6")
    )

    return bar5, metrics, figures, objs, actual_hamming, skip_hamming


def _compute_barcode_stats(final_series, whitelist: pd.Series | None = None) -> tuple[BarcodeStats, List[Dict]]:
    """Compute barcode statistics after final filtering."""
    
    # Cell counts per barcode
    barcode_cell_counts = final_series.groupby("Barcode").apply(
        lambda x: len(x.index.get_level_values("CBC").unique())
    )
    
    # UMI counts per barcode
    barcode_umi_counts = final_series.groupby("Barcode").sum()
    
    stats = BarcodeStats(
        barcode_cell_counts=barcode_cell_counts,
        barcode_umi_counts=barcode_umi_counts,
        total_unique_barcodes=len(barcode_cell_counts),
        total_cells=len(final_series.index.get_level_values("CBC").unique())
    )
    
    # Build table data
    table_data = []
    sorted_barcodes = barcode_cell_counts.sort_values(ascending=False)
    
    for bc in sorted_barcodes.index:
        cells = int(barcode_cell_counts.get(bc, 0))
        umis = int(barcode_umi_counts.get(bc, 0))
        in_wl = whitelist is not None and bc in whitelist.values
        table_data.append({
            "barcode": bc,
            "cells": cells,
            "umis": umis,
            "mean_umis_per_cell": umis / cells if cells > 0 else 0,
            "in_whitelist": in_wl
        })
    
    return stats, table_data


# ────────────────────────────────────────────────────────────────────────
#  Public CLI

def larry_qc(
    pkl_file: str,
    sample_id: str,
    *,
    min_reads: int = 8,
    min_hamming: int = 3,
    min_umis: int = 3,
    auto_hamming: bool = False,
    skip_hamming: bool = False,
    whitelist_csv: str = "/opt/cbutools/larry_whitelist.csv",
):
    """Run QC & generate an HTML report for one library.
    
    Args:
        pkl_file: Path to input pickle file containing CBUSeries
        sample_id: Sample identifier for output files
        min_reads: Minimum reads per CBC-Barcode pair
        min_hamming: Minimum Hamming distance (or default for auto mode)
        min_umis: Minimum UMIs per CBC-Barcode pair
        auto_hamming: If True, automatically detect optimal Hamming threshold from bimodal distribution
        skip_hamming: If True, calculate distances but skip Hamming filtering
        whitelist_csv: Path to LARRY whitelist CSV file
    """

    # Set up logging
    logging.basicConfig(
        filename=f"{sample_id}.log",
        filemode="a",
        level=logging.INFO,
        format="%(asctime)s  %(levelname)s: %(message)s",
    )
    logging.info("Starting LARRY_QC")
    if auto_hamming:
        logging.info("Auto Hamming threshold detection enabled")
    if skip_hamming:
        logging.info("Hamming distance filtering SKIPPED (distances will still be calculated)")

    # Load raw data
    with open(pkl_file, "rb") as handle:
        bar_m = pickle.load(handle)

    # Run filter cascade
    bar5, metrics, figures, objs, actual_hamming, skip_hamming_used = _run_filters(
        bar_m, min_reads=min_reads, min_hamming=min_hamming, min_umis=min_umis, auto_hamming=auto_hamming, skip_hamming=skip_hamming
    )

    # Load whitelist
    try:
        whitelist = pd.read_csv(whitelist_csv)["barcode"]
    except Exception as e:
        logging.warning(f"Could not load whitelist: {e}")
        whitelist = None

    # Whitelist percentages
    def pct_in_whitelist(obj):
        if whitelist is None:
            return None
        uniq = obj.index.get_level_values("Barcode").unique()
        return 100 * uniq.isin(whitelist).sum() / len(uniq) if len(uniq) > 0 else 0

    for m, obj in zip(metrics, objs):
        m.pct_whitelist = pct_in_whitelist(obj)
    logging.info("Whitelist stats calculated")

    # Compute barcode statistics
    barcode_stats, barcode_table_data = _compute_barcode_stats(bar5, whitelist)
    
    # Generate barcode figures
    figures["bc_cells"] = _fig_to_base64(_barcode_cell_bar_chart(barcode_stats))
    figures["bc_scatter"] = _fig_to_base64(_barcode_umi_scatter(barcode_stats))

    # Report generation
    filter_params = {
        "min_reads": min_reads,
        "min_hamming": actual_hamming,
        "min_hamming_skipped": skip_hamming_used,
        "min_umis": min_umis,
    }

    html_name = f"{sample_id}_QC_report.html"
    _write_html_report(
        html_name,
        sample_id=sample_id,
        metrics=metrics,
        figures=figures,
        filter_params=filter_params,
        barcode_stats=barcode_stats,
        barcode_table_data=barcode_table_data,
    )
    logging.info(f"HTML report written -> {html_name}")

    # Save final object
    out_pkl = f"{sample_id}-CBU_bar5.pkl"
    with open(out_pkl, "wb") as handle:
        pickle.dump(bar5, handle)
    logging.info(f"Final CBSeries pickled -> {out_pkl}")

    return out_pkl


if __name__ == "__main__":
    fire.Fire(larry_qc)
