#! /bin/bash

#$ -cwd

#$ -q broad
#$ -P regevlab
#$ -l h_vmem=60g
#$ -e ErrFiles/err.agg.scna.err
#$ -o ErrFiles/log.agg.scna.out
#$ -l h_rt=96:00:00
#$ -l os=RedHat7

source /broad/software/scripts/useuse
use Google-Cloud-SDK
use .samtools-1.8
use BEDTools
use R-3.5
use UGER
use .python-2.7.9-sqlite3-rtrees
use .git-2.12.0 
use .seqtk-1.0
use .java-jdk-1.8.0_181-x86-64
use .jags-4.3.0 
use .hdf5-1.8.16


export R_LIBS_USER=/stanley/levin_dr/ssimmons/R/x86_64-pc-linux-gnu-library/3.5_v2


Rscript Central.Test.R ../data/SCNA.RDS condition CellType orig.ident Inhibitory ../results/results.SCNA.txt
