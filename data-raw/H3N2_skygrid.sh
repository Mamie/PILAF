#!/bin/bash
#
#SBATCH --job-name=h3n2skygrid
#SBATCH --output=H3N2SkygridCutoff8.out
#SBATCH --error=H3N2SkygridCutoff8.err
#SBATCH --time=2-00:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=64000
#SBATCH --mail-type=END,FAIL

module load java
srun ~/BEASTv1.8.4/bin/beast -threads 8 ~/repos/PILAF/inst/extdata/H3N2/NewYork_A_H3N2_20190307_filtered_aligned.xml
