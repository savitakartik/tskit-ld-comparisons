#! /bin/bash

#SBATCH --job-name=calc_ld
#SBATCH --partition=compute
#SBATCH --mem=60G
#SBATCH --array=17-18
#SBATCH --output=logs/slurm/calc_ld/win_500kb/win_75_100/%a.out
#SBATCH --error=logs/slurm/calc_ld/win_500kb/win_75_100/%a.err
#SBATCH --time=300:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=osvk@novonordisk.com

. /etc/profile.d/modules.sh
module load plink/1.90
conda activate ~/.conda/envs/msprime-tskit-env/

date

# do all win sizes for win_kb=500 and 1000 as well
i=${SLURM_ARRAY_TASK_ID}
i=$((i-1))
arr=(10 100 1000 10000 25000 50000 75000 100000 250000 500000)
win_sizes=(75 100)
sample_sizes=($(for i in ${arr[@]}; do for ((n=1; n<=2; n++)) do echo -n "$i ";done; done;))
#win_sizes=( ${arr2[*]} ${arr2[*]} ${arr2[*]})
sample_size=${sample_sizes[i]}
#win_size=${win_sizes[i]}
win_size=${win_sizes[i%2]}
echo -ne "Sample size: ${sample_size} and window size: ${win_size}\n"
rm logs/compute_times/r2_array_method/ld_calculation_times_N${sample_size}_win${win_size}.csv
python compare_ld_calculation_times.py -n ${sample_size} -w ${win_size}

date