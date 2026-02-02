from pathlib import Path
import pandas as pd
import fire
import scanpy as sc
import pickle
import logging


def open_from_pickle(pkl_file):
    with open(pkl_file,"rb") as handle:
        bar_m = pickle.load(handle)
    return bar_m
	
DEFAULT_CB_FILENAMES = ["cellbender_filtered.h5", "cellbender_out_filtered.h5"]


def _cb_filenames_from_param(cb_filename: str) -> list[str]:
    cb_filename = (cb_filename or "").strip()
    if not cb_filename or cb_filename in DEFAULT_CB_FILENAMES:
        return DEFAULT_CB_FILENAMES
    return [cb_filename]


def detect_gex_source(
    gex_root: Path,
    ss_out: str = "Gene",
    cb_filename: str = "cellbender_out_filtered.h5",
    samples: dict | None = None,
) -> dict:
    """
    Detect the layout of the GEX data.

    Returns a dict with keys:
    - kind: "mtx" for Cell Ranger / STARsolo matrix folders, "cellbender" for
      CellBender h5 outputs.
    - rel_path: relative path (for mtx) or filename (for h5) to use under each
      sample directory. If gex_root is a single h5 file, rel_path is that name
      and single_file is True.
    - single_file: True when gex_root itself is an h5 (single sample CellBender).
    """
    # Single CellBender h5 provided directly
    if gex_root.is_file() and gex_root.suffix in {".h5", ".hdf5"}:
        return {"kind": "cellbender", "rel_path": gex_root.name, "single_file": True}

    allowed = set(samples.values()) if samples else None
    subdirs = [d for d in sorted(gex_root.iterdir()) if d.is_dir() and (allowed is None or d.name in allowed)]
    if not subdirs:
        if allowed:
            raise FileNotFoundError(
                f"No subdirectories in {gex_root} match provided samples: {sorted(allowed)}"
            )
        raise FileNotFoundError(f"No subdirectories found in {gex_root}")
    first_subdir = subdirs[0]

    # Cell Ranger default
    cr_rel = Path("outs/filtered_feature_bc_matrix")
    # STARsolo (e.g. sample/output/<ss_out>/filtered)
    ss_rel = Path("output") / ss_out / "filtered"

    if (first_subdir / cr_rel).is_dir():
        return {"kind": "mtx", "rel_path": cr_rel, "single_file": False}
    if (first_subdir / ss_rel).is_dir():
        return {"kind": "mtx", "rel_path": ss_rel, "single_file": False}

    # CellBender output h5 inside sample directory (require explicit filename)
    cb_names = _cb_filenames_from_param(cb_filename)
    cb_candidates = []
    for name in cb_names:
        cb_name = Path(name).name
        cb_candidates.extend(list(first_subdir.glob(cb_name)))
        # Allow .hdf5 if caller gave .h5 default name
        if not cb_candidates and cb_name.endswith(".h5"):
            cb_candidates.extend(list(first_subdir.glob(cb_name.replace(".h5", ".hdf5"))))
        if cb_candidates:
            return {
                "kind": "cellbender",
                "rel_path": cb_candidates[0].name,
                "single_file": False,
            }

    if cb_names:
        raise FileNotFoundError(
            f"Expected CellBender file(s) {cb_names} (or .hdf5) under {first_subdir}"
        )
    raise FileNotFoundError(
        "Could not detect GEX source: expected Cell Ranger, STARsolo, or CellBender outputs"
    )

def match_gex(samples_larry, sample_csv, ss_out, group_id, res3_pkl, gex_path, plot_cumulative, cb_filename):
    gex_root = Path(gex_path)  # pathlib early
    res3_tabs = open_from_pickle(res3_pkl)
    samples = pd.read_csv(sample_csv, index_col=0)['sample_gex'].to_dict()

    samples_larry = samples_larry.split(',')

    samples = {s: samples[s] for s in samples_larry}

    if len(samples_larry) == 1:
        samp_gex = samples[samples_larry[0]]

    tmp = res3_tabs.drop('Barcode_tuple', axis=1)

    barcode_to_clone = {
        barcode: f"Clone_{i+1}"
        for i, barcode in enumerate(tmp["Barcode"].unique())
    }
    tmp["Clone"] = tmp["Barcode"].map(barcode_to_clone)

    gex_source = detect_gex_source(
        gex_root,
        ss_out=ss_out,
        cb_filename=cb_filename,
        samples=samples,
    )
    gex_rel_subdir = gex_source["rel_path"]
    is_cellbender = gex_source["kind"] == "cellbender"
    single_file_cb = gex_source.get("single_file", False)
    cb_names = _cb_filenames_from_param(cb_filename)

    if len(samples_larry) > 1:
        tmp['sample_larry'] = [i.split("_")[0] for i in tmp.index]
        tmp['sample_gex'] = tmp['sample_larry'].map(samples)
        tmp.index = [f"{i.split('_')[1]}-{j}" for i, j in zip(tmp.index, tmp['sample_gex'])]
        tmp['sample_larry'] = tmp.index.str.split('_').str[0]

        adatas = []
        if is_cellbender:
            # Iterate over sample directories containing CellBender h5 files
            for sample_dir in sorted(gex_root.iterdir()):
                if not sample_dir.is_dir():
                    continue
                sample_id = sample_dir.name
                if sample_id in samples.values():
                    h5_path = None
                    for name in cb_names:
                        candidate = sample_dir / Path(name).name
                        if candidate.is_file():
                            h5_path = candidate
                            break
                        if candidate.suffix == ".h5":
                            alt = candidate.with_suffix(".hdf5")
                            if alt.is_file():
                                h5_path = alt
                                break
                    if h5_path is None:
                        continue
                    adata_tmp = sc.read_10x_h5(h5_path)
                    adata_tmp.obs['barcodes'] = adata_tmp.obs_names
                    adata_tmp.obs['sanger_id'] = sample_id
                    adata_tmp.var['gene_symbols'] = adata_tmp.var_names
                    adata_tmp.var_names = adata_tmp.var['gene_ids']
                    adatas.append(adata_tmp)
            if not adatas:
                raise FileNotFoundError(
                    f"No CellBender h5 files matching {cb_names} under {gex_root}"
                )
        else:
            # Iterate over sample/*/<relative_subdir>
            for path in sorted(gex_root.glob(f"*/{gex_rel_subdir.as_posix()}")):
                sample_id = path.relative_to(gex_root).parts[0]  # first dir under root
                if sample_id in samples.values():
                    adata_tmp = sc.read_10x_mtx(path)
                    adata_tmp.obs['barcodes'] = adata_tmp.obs_names
                    adata_tmp.obs['sanger_id'] = sample_id
                    adata_tmp.var['gene_symbols'] = adata_tmp.var_names
                    adata_tmp.var_names = adata_tmp.var['gene_ids']
                    adatas.append(adata_tmp)
            if not adatas:
                raise FileNotFoundError(f"No GEX directories matched pattern */{gex_rel_subdir} under {gex_root}")

        adata = sc.concat(adatas)
        adata.obs_names = [
            f"{bc}-{sid}" for sid, bc in zip(adata.obs['sanger_id'], adata.obs['barcodes'])
        ]
        adata.var = adatas[0].var
    else:
        tmp['sample_larry'] = samples_larry[0]
        tmp['sample_gex'] = tmp['sample_larry'].map(samples)
        tmp.index = [f"{i}-{j}" for i, j in zip(tmp.index, tmp['sample_gex'])]
        if is_cellbender:
            if single_file_cb:
                h5_path = gex_root
            else:
                h5_path = None
                sample_dir = gex_root / samp_gex
                for name in cb_names:
                    candidate = sample_dir / Path(name).name
                    if candidate.is_file():
                        h5_path = candidate
                        break
                    if candidate.suffix == ".h5":
                        alt = candidate.with_suffix(".hdf5")
                        if alt.is_file():
                            h5_path = alt
                            break
            if h5_path is None or not h5_path.is_file():
                raise FileNotFoundError(
                    f"CellBender h5 not found for sample {samp_gex} using {cb_names} under {gex_root}"
                )
            adata = sc.read_10x_h5(h5_path)
        else:
            adata = sc.read_10x_mtx(gex_root / samp_gex / gex_rel_subdir)
        adata.obs_names = [f"{i}-{samp_gex}" for i in adata.obs_names]

    clones = tmp[['Clone', 'Barcode', 'Barcode_n']].reset_index(drop=True).set_index("Clone").drop_duplicates()
    clones['sing_or_mult'] = ['multi' if i > 1 else 'single' for i in clones['Barcode_n']]

    tmp1 = tmp[tmp.index.isin(adata.obs_names)].copy()

    adata.obs["Clone"] = tmp1["Clone"].reindex(adata.obs_names)

    if plot_cumulative:
        import matplotlib.pyplot as plt

        cumulative_counts = adata.obs['Clone'].value_counts().sort_values(ascending=False).cumsum()
        # Reset index to get a proper DataFrame for plotting
        cum_df = pd.concat([cumulative_counts, adata.obs['Clone'].value_counts()],axis=1)
        cum_df = cum_df.reset_index()
        cum_df.columns = ['Clone', 'Cumulative', 'Count']

        plot_df = cum_df[cum_df['Count'] > 1]

        # Plot
        plt.figure(figsize=(9, 7))
        plt.plot(plot_df['Cumulative'], marker='o', linestyle='-', linewidth=2, label='Cumulative Count')

        # Aesthetics
        plt.title('Cumulative Clone Frequency Distribution', fontsize=16)
        plt.xlabel('Clone Rank', fontsize=14)
        plt.ylabel('Cumulative Frequency', fontsize=14)
        plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.tight_layout()

        plt.savefig(f"{group_id}_cumulative_clone_frequency.png", dpi=300, bbox_inches='tight')

    adata.obs["Clone"] = adata.obs["Clone"].fillna("NA").astype(str)
    for i in clones.columns:
        adata.obs[i] = adata.obs['Clone'].map(clones[i])
        adata.obs[i] = adata.obs[i].fillna("NA").astype(str)

    adata.write(f"{group_id}_larry_gex_clones.h5ad", compression='gzip')
    clones.to_csv(f"{group_id}_larry_gex_clones.csv")
	

if __name__ == '__main__':
    fire.Fire(match_gex)