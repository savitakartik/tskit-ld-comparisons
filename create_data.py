import tskit
import tszip
import time
import os
import argparse
import gzip
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--num_indivs", type=int)
args = parser.parse_args()
num_indivs = args.num_indivs #num of indivs
print(f"input args: num of indivs ={num_indivs}")

def subset_trees(N, ts, trees_file, log_file, start, end, seed=123):
    np.random.seed(seed)
    start_time = time.perf_counter()
    indivs = np.random.choice(ts.individuals(), N, replace=False)
    ind_nodes = [] 
    for ind in indivs:
        ind_nodes.extend(ind.nodes)
    ts_sub= ts.simplify(samples=ind_nodes)
    ts_sub_gene = ts_sub.keep_intervals([[start,end]])
    non_biallelic_sites = np.where(np.bincount(ts_sub_gene.mutations_site, minlength=ts_sub_gene.num_sites) != 1)[0]
    ts_gene_biallelic = ts_sub_gene.delete_sites(non_biallelic_sites)
    ts_gene_biallelic.dump(trees_file)
    end_time = time.perf_counter()
    with open(log_file, "a") as log_file:
        log_file.write(f"Time to create and dump tree subset of {ts_gene_biallelic.num_samples} samples: {end_time - start_time:.2f} seconds\n")
    return ts_gene_biallelic

def ts_to_vcf(ts_to_conv, vcf_file, log_file):
    start = time.perf_counter()
    with gzip.open(vcf_file, "wt") as vcf_file:
        ts_to_conv.write_vcf(vcf_file)
    end = time.perf_counter()
    with open(log_file, "a") as log_file:
        log_file.write(f"Time to write VCF for {ts_to_conv.num_samples} samples: {end - start:.2f} seconds\n")

def vcf_to_plink(vcf_file, plink_file, log_file):
    start = time.perf_counter()
    os.system(f"plink --vcf {vcf_file} --out {plink_file} --double-id --vcf-half-call 'h'")
    end = time.perf_counter()
    with open(log_file, "a") as log_file:
        log_file.write(f"Time to convert VCF to PLINK for {vcf_file}: {end - start:.2f} seconds\n")

ts = tszip.decompress(
    "/nfs_home/users/osvk/projects/tskit-ld/data/simulated_chrom_22.ts.tsz"
)

ts_gene_biallelic = subset_trees(num_indivs, ts, 
                            f"data/trees/ts_gene_biallelic_N{num_indivs}.trees",
                            f"logs/conversions/file_conversions_N{num_indivs}.log",
                            start=49_000_000, end=50_000_000
                            )

ts_to_vcf(ts_gene_biallelic, f"data/vcf/ts_gene_biallelic_N{num_indivs}.vcf.gz", 
          log_file=f"logs/conversions/file_conversions_N{num_indivs}.log")
vcf_to_plink(f"data/vcf/ts_gene_biallelic_N{num_indivs}.vcf.gz", 
             f"data/bed/ts_gene_biallelic_N{num_indivs}",
             f"logs/conversions/file_conversions_N{num_indivs}.log")