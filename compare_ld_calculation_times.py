import numpy as np
import tskit
import time
import os
import argparse
import pandas as pd

#12 July 2023

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--num_indivs", type=int)
parser.add_argument("-w", "--win_size", type=int)

args = parser.parse_args()
num_indivs = args.num_indivs
win_size = args.win_size
win_kb_size = 500
print(f"input args: num of indivs={num_indivs} and window size={win_size}")
log_file=f"logs/compute_times/r2_array_method/win_kb_{win_kb_size}/ld_calculation_times_N{num_indivs}_win{win_size}.csv"

def update_log_file(log_file, num_samples, method, ld_window, ld_window_kb, time_taken):
    if not os.path.isfile(log_file):
        with open(log_file, 'w') as f:
            f.write("num_samples,method,ld_window,ld_window_kb,time_taken\n")
    with open(log_file, 'a') as f:
        f.write(f"{num_samples},{method},{ld_window},{ld_window_kb},{time_taken}\n")

def calc_ld_with_tskit(ts, ld_window, ld_window_kb, out_file=None):
    start = time.perf_counter()
    ldcalc = tskit.LdCalculator(ts)
    sites_pos = ts.sites_position

    df = pd.DataFrame()
    focal_sites=[]
    focal_positions=[]
    target_sites=[]
    target_positions=[]
    r2=[]

    for focal_site in range(ts_gene_biallelic.num_sites):
        r2_fwd = ldcalc.r2_array(focal_site, direction=tskit.FORWARD, max_mutations=ld_window)
        r2_rev = ldcalc.r2_array(focal_site, direction=tskit.REVERSE, max_mutations=ld_window)
        r2.append(np.round(r2_fwd,2))
        r2.append(np.round(r2_rev,2))
        focal_sites.append(np.repeat(focal_site,len(r2_fwd)+len(r2_rev)))
        focal_target_sites_fwd=np.arange(focal_site+1, focal_site+1+len(r2_fwd))
        target_sites.append(focal_target_sites_fwd)
        focal_target_sites_rev=np.arange(focal_site-1,focal_site-1-len(r2_rev),-1)
        target_sites.append(focal_target_sites_rev)
    target_sites=np.concatenate(target_sites)
    focal_sites=np.concatenate(focal_sites)
    focal_positions=np.take(sites_pos, focal_sites).astype(int).flatten()
    target_positions=np.take(sites_pos, target_sites).astype(int).flatten()
    df['focal_site']=focal_sites
    df['focal_position']=focal_positions
    df['target_site']=target_sites
    df['target_position']=target_positions
    df['r2']=np.concatenate(r2)
    
    if out_file is not None:
        df.to_pickle(out_file)
        #df.to_csv(out_file.replace(".pkl", ".csv"),index=False, sep="\t")
    update_log_file(log_file, num_indivs, "tskit", ld_window, ld_window_kb, 
                    time.perf_counter()-start)

def calc_ld_with_plink(bfile, out_file, ld_window, ld_window_kb):
    start = time.perf_counter()
    os.system(f"plink --bfile {bfile} --r2 --ld-window-r2 0 --ld-window {ld_window} --ld-window-kb {ld_window_kb} --threads 1 --out {out_file}")
    update_log_file(log_file, num_indivs, "plink", ld_window, ld_window_kb, time.perf_counter() - start)

ts_gene_biallelic= tskit.load(
    f"/nfs_home/users/osvk/projects/tskit-ld/data/trees/ts_gene_biallelic_N{num_indivs}.trees")
print(f"num sites in input trees: {ts_gene_biallelic.num_sites}")

plink_file = f"/nfs_home/users/osvk/projects/tskit-ld/data/bed/ts_gene_biallelic_N{num_indivs}"
calc_ld_with_tskit(ts_gene_biallelic, 
                   ld_window=win_size, ld_window_kb=win_kb_size,
                   out_file=f"data/tskit_ld/r2_array_method/win_kb_{win_kb_size}/N{num_indivs}_w{win_size}_ld.pkl",
                   )

calc_ld_with_plink(f"data/bed/ts_gene_biallelic_N{num_indivs}",
                    f"data/plink_ld/r2_array_method/win_kb_{win_kb_size}/N{num_indivs}_win{win_size}", 
                    ld_window=win_size, 
                    ld_window_kb=win_kb_size)