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
	
def detect_gex_source(gex_root: Path, ss_out: str = "Gene") -> Path:
    """
    Inspect the first sample directory to infer whether data come from Cell Ranger
    or STARsolo and return the relative subdirectory (under each sample) that
    contains the 10X-style matrix for sc.read_10x_mtx.
    """
    subdirs = [d for d in sorted(gex_root.iterdir()) if d.is_dir()]
    if not subdirs:
        raise FileNotFoundError(f"No subdirectories found in {gex_root}")
    first_subdir = subdirs[0]

    # Cell Ranger default
    cr_rel = Path("outs/filtered_feature_bc_matrix")
    # STARsolo (e.g. sample/output/<ss_out>/filtered)
    ss_rel = Path("output") / ss_out / "filtered"

    if (first_subdir / cr_rel).is_dir():
        return cr_rel
    if (first_subdir / ss_rel).is_dir():
        return ss_rel
    raise FileNotFoundError(
        f"Neither Cell Ranger ({cr_rel}) nor STARsolo ({ss_rel}) folder found under {first_subdir}"
    )

def match_gex(samples_larry, sample_csv, ss_out, group_id, res3_pkl, gex_path, plot_cumulative):
    gex_root = Path(gex_path)  # pathlib early
    res3_tabs = open_from_pickle(res3_pkl)
    samples = pd.read_csv(sample_csv, index_col=0)['sample_gex'].to_dict()

    samples_larry = samples_larry.split(',')

    if len(samples_larry) == 1:
        samp_gex = samples[samples_larry[0]]

    tmp = res3_tabs.drop('Barcode_tuple', axis=1)

    barcode_to_clone = {
        barcode: f"Clone_{i+1}"
        for i, barcode in enumerate(tmp["Barcode"].unique())
    }
    tmp["Clone"] = tmp["Barcode"].map(barcode_to_clone)

    gex_rel_subdir = detect_gex_source(gex_root, ss_out=ss_out)

    if len(samples_larry) > 1:
        tmp['sample_larry'] = [i.split("_")[0] for i in tmp.index]
        tmp['sample_gex'] = tmp['sample_larry'].map(samples)
        tmp.index = [f"{i.split('_')[1]}-{j}" for i, j in zip(tmp.index, tmp['sample_gex'])]
        tmp['sample_larry'] = tmp.index.str.split('_').str[0]

        adatas = []
        # Iterate over sample/*/<relative_subdir>
        for path in sorted(gex_root.glob(f"*/{gex_rel_subdir.as_posix()}")):
            sample_id = path.relative_to(gex_root).parts[0]  # first dir under root
            if sample_id in samples.values():
                adata_tmp = sc.read_10x_mtx(path)
                adata_tmp.obs['barcodes'] = adata_tmp.obs_names
                adata_tmp.obs['sanger_id'] = sample_id
                adatas.append(adata_tmp)
        if not adatas:
            raise FileNotFoundError(f"No GEX directories matched pattern */{gex_rel_subdir} under {gex_root}")
        adata = sc.concat(adatas)
        adata.obs_names = [
            f"{bc}-{sid}" for sid, bc in zip(adata.obs['sanger_id'], adata.obs['barcodes'])
        ]
    else:
        tmp['sample_larry'] = samples_larry[0]
        tmp['sample_gex'] = tmp['sample_larry'].map(samples)
        tmp.index = [f"{i}-{j}" for i, j in zip(tmp.index, tmp['sample_gex'])]
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