import pandas as pd
import fire
import json
import glob
import scanpy as sc
import pickle
from collections import defaultdict

def open_from_pickle(pkl_file):
    with open(pkl_file,"rb") as handle:
        bar_m = pickle.load(handle)
    return bar_m

def open_from_json(json_file):
    with open(json_file,"r") as f:
        data = json.load(f)
    return data

def match_gex(sampl, bar5_pkl, samples_json, gex_path, combine_samples, plot_cumulative):
	res3_tabs = open_from_pickle(bar5_pkl)
	samples = open_from_json(samples_json)

	if not len(sampl.split(',')) > 1:
		samp_gex = samples[sampl]

	tmp = res3_tabs.drop('Barcode_tuple', axis=1)

	barcode_to_clone = {
		barcode: f"Clone_{i+1}" 
		for i, barcode in enumerate(tmp["Barcode"].unique())
	}

	tmp["Clone"] = tmp["Barcode"].map(barcode_to_clone)

	

	if combine_samples == 'true':
		tmp['sample_larry'] = [i.split("_")[0] for i in tmp.index]
		tmp['sample_gex'] = tmp['sample_larry'].map(samples)
		tmp.index = [f"{j}_{i.split('_')[1]}" for i, j in zip(tmp.index, tmp['sample_gex'])]
		tmp['sample_larry'] = tmp.index.str.split('_').str[0]
		if all(i.endswith("h5ad") for i in glob.glob(f"{gex_path}/*")): # h5ads are ready
			# assuming all h5ad files include sanger_id and barcodes columns
			adata = sc.concat([sc.read(i) for i in sorted(glob.glob(f"{gex_path}/*.h5ad")) if i.split('/')[-1].startswith(tuple(samples.values()))])
		else: # h5ads are not ready
			adata = []
			for i in sorted(glob.glob(f"{gex_path}/*/output/Gene/filtered")):
				samp = i.split('/')[-4]
				if samp.startswith(tuple(samples.values())):
					adata_tmp = sc.read_10x_mtx(i)
					adata_tmp.obs['barcodes'] = adata_tmp.obs_names
					adata_tmp.obs['sanger_id'] = samp
					adata.append(adata_tmp)
		adata.obs_names = [f"{i}_{j}"for i, j in zip(adata.obs['sanger_id'], adata.obs['barcodes'])]
	else:
		tmp['sample_larry'] = sampl
		if all(i.endswith("h5ad") for i in glob.glob(f"{gex_path}/*")):
			h5ad_path = glob.glob(f"{gex_path}/{samp_gex}*.h5ad")[0]
			adata = sc.read(h5ad_path)
		else:
			adata = sc.read_10x_mtx(f"{gex_path}/{sampl}/output/Gene/filtered")

	tmp['sample_gex'] = tmp['sample_larry'].map(samples)
	
	clones = tmp[['Clone', 'Barcode', 'Barcode_n']].reset_index(drop=True).set_index("Clone").drop_duplicates()
	clones['sing_or_mult'] = ['multi' if i > 1 else 'single' for i in clones['Barcode_n']]

	tmp1 = tmp[tmp.index.isin(adata.obs_names)].copy()

	adata.obs["Clone"] = tmp1["Clone"].reindex(adata.obs_names)

	if plot_cumulative == 'true':
		import matplotlib.pyplot as plt

		cumulative_counts = adata.obs['Clone'].value_counts().sort_values(ascending=False).cumsum()
		# Reset index to get a proper DataFrame for plotting
		cum_df = pd.concat([cumulative_counts, adata.obs['Clone'].value_counts()],axis=1)
		cum_df = cum_df.reset_index()
		cum_df.columns = ['Clone', 'Cumulative', 'Count']

		plot_df = cum_df[cum_df['Count'] > 1]

		clones_single = clones[clones['sing_or_mult'] == 'single']
		clones_multi = clones[clones['sing_or_mult'] == 'multi']

		# Preprocess: build a mapping from each barcode fragment to the full row indices of clones_multi
		barcode_to_mult_rows = defaultdict(list)

		for idx, row in clones_multi.iterrows():
			for part in row['Barcode'].split('-'):
				barcode_to_mult_rows[part].append(idx)

		# Build the rows
		rows = []

		for cl_sing_idx, row_sing in clones_single.iterrows():
			barcode = row_sing['Barcode']
			matching_idxs = barcode_to_mult_rows.get(barcode, [])
			clones_multi_similar = list(clones_multi.loc[matching_idxs].index)
			
			rows.append({
				'Clone_single': cl_sing_idx,
				'Clones_multi_similar': clones_multi_similar,
				'n_clones_similar': len(clones_multi_similar)
			})

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

		if len(sampl.split(',')) > 1:
			plt.savefig("cumulative_clone_frequency.png", dpi=300, bbox_inches='tight')
		else:
			plt.savefig(f"{sampl}_cumulative_clone_frequency.png", dpi=300, bbox_inches='tight')

	adata.obs["Clone"] = adata.obs["Clone"].fillna("NA").astype(str)
	for i in clones.columns:
		adata.obs[i] = adata.obs['Clone'].map(clones[i])
		adata.obs[i] = adata.obs[i].fillna("NA").astype(str)

	if len(sampl.split(',')) > 1:
		adata.write("larry_gex_clones.h5ad", compression='gzip')
		clones.to_csv("larry_gex_clones.csv")
	else:
		adata.write(f"{sampl}_larry_gex_clones.h5ad", compression='gzip')
		clones.to_csv(f"{sampl}_larry_gex_clones.csv")
	
	


if __name__ == '__main__':
    fire.Fire(match_gex)