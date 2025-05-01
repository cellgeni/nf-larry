import os
import sys
sys.path.append("/opt")
import cbutools as cbu
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import fire
import logging

def open_from_pickle(sample_i):
    with open(sample_i,"rb") as handle:
        bar_m = pickle.load(handle)
    return bar_m

def check_whitelist(CBUseries, whitelist, step):
    logging.info(f"~~~~~~~ {step}")
    len_uniq = len(CBUseries.index.get_level_values("Barcode").unique())
    len_uniq_white = CBUseries.index.get_level_values("Barcode").unique().isin(whitelist).sum()
    logging.info(f"Number of unique barcodes: {len_uniq}")
    logging.info(f"Number of unique barcodes in whitelist: {len_uniq_white}")
    logging.info(f"{'{:.2f}'.format(100*(len_uniq_white/len_uniq))}% of the barcodes from whitelist")
	
def qc_LARRY_69(pkl_file,sample_i,check_whitelist,plot_qc):
    logging.basicConfig(
        filename=f"{sample_i}_whitelist_check.log",
        filemode='a',
        format='%(asctime)s %(levelname)s: %(message)s',
        level=logging.INFO
    )
    
    bar_m = open_from_pickle(pkl_file)
	
    logging.basicConfig(
		filename=f"{sample_i}_whitelist_check.log",
		filemode='a',
		format='%(asctime)s %(levelname)s: %(message)s',
		level=logging.INFO
	)  
    
    bar2 = bar_m.filter_by_reads(groupby=["CBC", "Barcode"], min_counts=8) # maybe cut from 100-200, the parameter is not logged
    bar3 = bar2.filter_by_hamming(min_distance=3) # 8-10 should be fine
    bar4 = bar3.count_UMI() 
    bar5 = bar4.filter_by_UMI(groupby=["CBC", "Barcode"], min_counts=3)
   
    if plot_qc == 'true':
        bar_m.plot_hist(groupby=["CBC", "Barcode"])
        bar2.plot_hist(groupby=["CBC", "Barcode"])
        bar2.plot_barcode_no()
        bar2.summary()
        fig,ax = plt.subplots(figsize=(4,4))
        bar3.plot_barcode_no()
        bar3.summary()
        fig,ax = plt.subplots(figsize=(4,4))
        bar4.plot_hist(groupby=["Barcode", "CBC"])
        bar4.summary()
        fig,ax = plt.subplots(figsize=(4,4))
        bar5.plot_hist(groupby=["Barcode", "CBC"])
        fig,ax = plt.subplots(figsize=(4,4))
        bar5.plot_barcode_no()
        bar5.summary()

    if check_whitelist == 'true':
        whitelist = pd.read_csv('/opt/cbutools/larry_whitelist.csv')['barcode']
        check_whitelist(bar_m, whitelist, 'no filtering')
        check_whitelist(bar2, whitelist, 'filtering counts < 8')
        check_whitelist(bar3, whitelist, 'filtering by hamming distance < 3')
        check_whitelist(bar4, whitelist, 'counting UMIs')
        check_whitelist(bar5, whitelist, 'filtering UMIs < 3')

    with open(f"{sample_i}-CBU_bar5.pkl","wb") as handle:
        pickle.dump(bar5,handle)
    
if __name__ == '__main__':
    fire.Fire(qc_LARRY_69)