#library(Seurat);
library(dplyr);library(tidyr)
source("RunAll.no.multi.R")

print("Loaded!")


testResults_no_multi<-function(seur,condition,celltype,basetype,individual)
{
meta=c()
if(is.data.frame(seur))
{meta=seur}
else
{
meta=seur@meta.data
}

lst=c("ko","wt")
meta["condition"]=lst[as.numeric(factor(meta[,condition]))]
meta["CellType"]=meta[,celltype]
meta["ID_loc"]=meta[,individual]

meta["CellType"]=as.character(meta[,"CellType"])

perc=table(meta[,"CellType"])/dim(meta)[1]

#meta[meta[,"CellType"] %in% names(perc)[perc<.005],"CellType"]="rare"
meta[,"CellType"]=factor(meta[,"CellType"])

#meta=permuteIndividual(meta,individual,celltype,condition,condition_new)
print(head(meta))
dat=RunAll(meta,"condition","CellType",basetype,"ID_loc")
return(dat)
}



if(!interactive())
{
args = commandArgs(trailingOnly=TRUE)
seur=readRDS(args[1])
condition=args[2]
celltype=args[3]
individual=args[4]
basetype=args[5]
print("loaded")
ret=testResults(seur,condition,celltype,basetype,individual)
savefil=args[6]
write.table(ret,savefil,row.names=F,sep="\t",quote=F)
}









