#!/bin/bash
#
#SBATCH --job-name="Process_New_Ablate_Data"
#SBATCH -o Process_New_Ablate_Data.out
#SBATCH -e Process_New_Ablate_Data.err
#SBATCH --mail-user=dwyerdei@buffalo.edu
#SBATCH --mail-type=end
#SBATCH --cluster=ub-hpc
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --mem=100G
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=30
#SBATCH --time=00:20:00

source ~/.bashrc
conda activate

#export data_dir=~/data/_1DSampleDiffusionFlame/flowField_mixtureFraction
#export data_dir=~/ClionProjects/ablate/ablateInputs/chemtabDiffusionFlame/_1DSampleDiffusionFlame/flowField_mixtureFraction
export data_dir=~/ClionProjects/ablate/ablateInputs/chemtabDiffusionFlame/_1DSampleChemTabDiffusionFlame/flowField_chemTab
export should_rebalance=T
post_glob="*.?????.csv.gz"
[ $should_rebalance ] && post_glob="*-rebalanced.csv.gz"

# optionally remove all previous files (e.g. if they are defective/outdated)
rm $data_dir/$post_glob

. wrangle_all_ablate_data.sh

Rscript Collate_Ablate_Data.R $data_dir/$post_glob
