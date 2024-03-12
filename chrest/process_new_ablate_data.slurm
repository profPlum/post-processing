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

# a little automation to select most recent simulation run!
#data_dir=~/ClionProjects/ablate/ablateInputs/chemtabDiffusionFlame/_1DSampleChemTabDiffusionFlame*/flowField_chemTab
#data_dir=~/ClionProjects/ablate/ablateInputs/chemtabDiffusionFlame/_1DSampleChemTabDiffusionFlame-R^2=1*/flowField_chemTab
data_dir=~/ClionProjects/ablate/ablateInputs/chemtabDiffusionFlame/_1DSampleDiffusionFlame*/flowField_mixtureFraction
#data_dir=~/CLionProjects/ablate/ablateInputs/chemTabDiffusionFlame/_tchemMethaneDiffusionFlameResults/flowField_mixtureFraction
#data_dir=~/ClionProjects/ablate/ablateInputs/chemtabDiffusionFlame/_1DSampleChemTabDiffusionFlame-Resume*/flowField_chemTab
export data_dir=$(ls -dt $data_dir | tr ' ' '\n' | head -n 1)
echo data_dir=$data_dir

export should_rebalance=0
post_glob="*.?????.csv.gz"
((should_rebalance)) && post_glob="*-rebalanced.csv.gz"
#rm $data_dir/$post_glob # optionally remove all previous files (e.g. if they are defective/outdated)

. wrangle_all_ablate_data.sh

Rscript Collate_Ablate_Data.R $data_dir/$post_glob
mkdir ablate_collated_new
mv ablate_collated.csv.gz ablate_collated_new/.
Rscript FitCPVs.R ablate_collated_new/ablate_collated.csv.gz
echo ^ NOTE: it\'s ok if the FitCPVs.R script failed above

if (( $(/bin/ls ablate_collated_new | wc -l)==1 )); then
    mv ablate_collated_new/ablate_collated.csv.gz ablate_collated_$RANDOM.csv.gz
    rm -r ablate_collated_new
else
    mv ablate_collated_new ablate_collated_$RANDOM
fi
echo done
