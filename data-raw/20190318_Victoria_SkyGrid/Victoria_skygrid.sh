#!/bin/bash
#
#SBATCH --job-name=Victoriaskygrid
#SBATCH --output=VictoriaSkygridCutoff8.out
#SBATCH --error=VictoriaSkygridCutoff8.err
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --cores=4
#SBATCH --mem=25600
#SBATCH --mail-type=END,FAIL

module load java
ulimit -c unlimited
beast -overwrite -threads 4 ~/repos/PILAF/inst/extdata/Victoria/NewYork_B_Victoria_20190318_filtered_aligned.xml
