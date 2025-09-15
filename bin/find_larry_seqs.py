import sys
import cbutools as cbu
import pickle
import fire

def pro_sample_i_v2(fq_dir,r1_fq,r2_fq,sample_i):
    fq_dir = fq_dir
    bar_m = cbu.get_barcodes(files={'r1' : r1_fq, 'r2' : r2_fq})
    opt_file = f"{sample_i}-CBU_bar_m.pkl"
    with open(opt_file,"wb") as handle:
        pickle.dump(bar_m,handle)

if __name__ == '__main__':
    fire.Fire(pro_sample_i_v2)