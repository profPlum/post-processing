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

source ~/.bash_profile
conda activate

#export data_dir=~/data/_1DSampleDiffusionFlame/flowField_mixtureFraction
#export data_dir=~/ClionProjects/ablate/ablateInputs/chemtabDiffusionFlame/_1DSampleDiffusionFlame/flowField_mixtureFraction
#export data_dir=~/ClionProjects/ablate/ablateInputs/chemtabDiffusionFlame/_1DSampleChemTabDiffusionFlame/flowField_chemTab
#export data_dir=~/ClionProjects/ablate/ablateInputs/chemtabDiffusionFlame/_1DSampleDiffusionFlame.10504.bak/flowField_mixtureFraction

# a little automation to select most recent simulation run!
data_dir=~/ClionProjects/ablate/ablateInputs/chemtabDiffusionFlame/_1DSampleChemTabDiffusionFlame_*/flowField_chemTab
#data_dir=~/ClionProjects/ablate/ablateInputs/chemtabDiffusionFlame/_1DSampleDiffusionFlame*/flowField_mixtureFraction
export data_dir=$(ls -dt $data_dir | tr ' ' '\n' | head -n 1)
echo data_dir=$data_dir

export should_rebalance=0
#export should_rebalance=1
post_glob="*.?????.csv.gz"
((should_rebalance)) && post_glob="*-rebalanced.csv.gz"
rm $data_dir/$post_glob # optionally remove all previous files (e.g. if they are defective/outdated)

. wrangle_all_ablate_data.sh

Rscript Collate_Ablate_Data.R $data_dir/$post_glob
mkdir ablate_collated.new
mv ablate_collated.csv.gz ablate_collated.new/.
Rscript FitCPVs.R ablate_collated.new/ablate_collated.csv.gz
mv ablate_collated.new ablate_collated.$RANDOM
