#!/bin/bash
#
#SBATCH --job-name=h1n1skygrid
#SBATCH --output=H1N1SkygridCutoff8.out
#SBATCH --error=H1N1SkygridCutoff8.err
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --cores=4
#SBATCH --mem=25600
#SBATCH --mail-type=END,FAIL

module load java
ulimit -c unlimited
beast -overwrite -threads 4 ~/repos/PILAF/inst/extdata/H1N1/NewYork_A_H1N1_20190318_filtered_aligned.xml
