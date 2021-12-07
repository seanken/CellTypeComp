#library(Seurat);
library(dplyr);library(tidyr)
source("RunAll.R")

##
##Permutes the individuals so 
##Here, meta is the metadata from Seurat object
##ID_loc is the column containing individuals ID
##condition is the column containing condition info (case vs control, etc)
##condition_new is the column of perturbed conditions
##
permuteIndividual_power<-function(meta,ID_loc,CellType,condition,condition_new,perc)
{
meta=meta[,c(ID_loc,condition,CellType)]
colnames(meta)=c("ID_loc","condition","CellType")
tab<-meta %>% group_by(ID_loc,condition) %>% summarise() %>% as.data.frame()
tab[condition_new]=sample(tab[,"condition"])
tab["condition"]=tab[,condition_new]
meta=meta[,c("ID_loc","CellType")]
meta=left_join(meta,tab)
lst=c("Excitatory","Gone")
num=sum(meta[,"condition"]!="wt" & meta[,"CellType"]=="Excitatory")-1
meta[meta[,"condition"]!="wt" & meta[,"CellType"]=="Excitatory","CellType"]=lst[c(1,(1+rbinom(num,1,perc)))]
meta=meta[meta[,"CellType"]!="Gone",]
return(meta)
}



testResults_power<-function(seur,condition,celltype,basetype,individual,perc,numPerm=1000)
{
condition_new="cond_temp"
meta=c()
if(is.data.frame(seur))
{meta=seur}
else
{
meta=seur@meta.data
}

ret=c()
for(i in 1:numPerm)
{

if(is.data.frame(seur))
{meta=seur}
else
{
meta=seur@meta.data
}

meta=permuteIndividual_power(meta,individual,celltype,condition,condition_new,perc)
print(head(meta))
print(i)
dat=RunAll(meta,"condition","CellType",basetype,"ID_loc")
dat["iter"]=i
if(i==1){ret=dat}
else{ret=rbind(ret,dat)}
}
return(ret)
}



if(!interactive())
{
args = commandArgs(trailingOnly=TRUE)
seur=readRDS(args[1])
condition=args[2]
celltype=args[3]
individual=args[4]
basetype=args[5]
perc=.5
if(length(args)>6){perc=as.numeric(args[7])}
ret=testResults(seur,condition,celltype,basetype,individual,perc)
savefil=args[6]
write.table(ret,savefil,row.names=F,sep="\t",quote=F)
}









