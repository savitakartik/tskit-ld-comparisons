#! /bin/bash

#SBATCH --job-name=create_data
#SBATCH --partition=compute
#SBATCH --mem=60G
#SBATCH --array=9-10
#SBATCH --output=logs/slurm/create_data/%a.out
#SBATCH --error=logs/slurm/create_data/%a.err
#SBATCH --time=300:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=osvk@novonordisk.com

. /etc/profile.d/modules.sh
module load plink/1.90
conda activate ~/.conda/envs/msprime-tskit-env/

date

i=${SLURM_ARRAY_TASK_ID}
i=$((i-1))
num_indivs=(10 100 1000 10000 25000 50000 75000 100000 250000 500000)
num_indivs=${num_indivs[i]}
echo -ne "Num of indivs: ${num_indivs}\n"
python create_data.py -n ${num_indivs}

date