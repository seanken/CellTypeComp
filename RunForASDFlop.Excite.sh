#! /bin/bash

#$ -cwd

#$ -q broad
#$ -P regevlab
#$ -l h_vmem=60g
#$ -e ErrFiles/err.agg.err
#$ -o ErrFiles/log.agg.out
#$ -l h_rt=96:00:00
#$ -l os=RedHat7

source /broad/software/scripts/useuse
use R-3.5



export R_LIBS_USER=/stanley/levin_dr/ssimmons/R/x86_64-pc-linux-gnu-library/3.5_v2


Rscript Central.Test.flop.excite.R ../data/ASD/ASD.RDS condition CellType individual Inhibitory ../results/results.ASD.excite.txt
