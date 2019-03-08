#!/bin/bash
#
#SBATCH --job-name=h3n2skygrid
#SBATCH --output=H3N2SkygridCutoff8.out
#SBATCH --error=H3N2SkygridCutoff8.err
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --cores=4
#SBATCH --mem=25600
#SBATCH --mail-type=END,FAIL

module load java
ulimit -c unlimited
beast -overwrite -threads 4 ~/repos/PILAF/inst/extdata/H3N2/NewYork_A_H3N2_20190307_filtered_aligned.xml
