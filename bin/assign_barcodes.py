import sys
import cbutools as cbu
import pandas as pd
import fire
import pickle

def open_from_pickle(sample_i):
    with open(sample_i,"rb") as handle:
        bar_m = pickle.load(handle)
    return bar_m

def assign_barcodes_examine_combination(bar_5, dispr_filter=None):

    """Assign barcodes to cell

        Loops through each cell and assigns the barcode(s) to each one

        Parameters
        ----------
        dispr_filter : float
            If None all barcodes are assigned to the cell. Otherwise all
            barcodes < dispr_filter * max_count are rejected. max_count
            correspond to the counts observed for the most abundant barcode

        Returns
        -------
        pd.DataFrame
            DataFrame with barcodes assigned per cell

        """

    if len(bar_5.split(',')) > 1:
        res_tabs2 = sorted(bar_5.split(','))
        res_tabs2_new = []
        for bar5 in res_tabs2:
            samp = bar5.split("_")[0]
            s = open_from_pickle(bar5)	
            new_cbc = s.index.get_level_values(0).map(lambda x: f"{samp}_{x}")
            new_index = pd.MultiIndex.from_arrays(
                [new_cbc, s.index.get_level_values(1)],
                names=s.index.names
            )
            s.index = new_index
            res_tabs2_new.append(s)
        bar5_pkl = pd.concat(res_tabs2_new)
    else:
        bar5_pkl = open_from_pickle(bar_5)	

    df = pd.DataFrame()
    for i in bar5_pkl.index.get_level_values("CBC").unique():
        counts = bar5_pkl[i]
        if dispr_filter is not None:
            max_counts = max(counts)
            counts = counts[counts >= (dispr_filter * max_counts)]
        n_total_reads = counts.sum()
        barcodes = tuple(counts.index.tolist())
        row = pd.DataFrame(
            {"Barcode_tuple": [barcodes], "Barcode_n": len(barcodes),"UMI_count":n_total_reads},index=[i]
        )
        df = pd.concat([df, (row)])
    def catl(x):
        """Function for concatenating strings"""
        return sorted(x,reverse=False)
        
    def str_barcode(x):
        return "-".join(x)
    
    df["Barcode_list"] = df["Barcode_tuple"].apply(func=catl)
    df["Barcode"] = df["Barcode_list"].apply(func=str_barcode)
    
    with open('res3_tabs.pkl', 'wb') as f:
        pickle.dump(df, f)


if __name__ == '__main__':
    fire.Fire(assign_barcodes_examine_combination)
