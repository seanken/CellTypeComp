#library(Seurat);
library(dplyr);library(tidyr)
source("RunAll.new_overdisp.R")

##
##Permutes the individuals so 
##Here, meta is the metadata from Seurat object
##ID_loc is the column containing individuals ID
##condition is the column containing condition info (case vs control, etc)
##condition_new is the column of perturbed conditions
##
permuteIndividual<-function(meta,ID_loc,CellType,condition,condition_new)
{
meta=meta[,c(ID_loc,condition,CellType)]
colnames(meta)=c("ID_loc","condition","CellType")
numSamp=floor(dim(meta)[1]/4)
meta=meta[sample(1:dim(meta)[1],numSamp),]
#tab<-meta %>% group_by(ID_loc,condition) %>% summarise() %>% as.data.frame()
#tab[condition_new]=sample(tab[,"condition"])
#tab["condition"]=tab[,condition_new]
#meta=meta[,c("ID_loc","CellType")]
#meta=left_join(meta,tab)
meta["CellType"]=factor(as.character(meta[,"CellType"]))
lst=c("ko","wt")
meta["condition"]=lst[as.numeric(factor(meta[,"condition"]))]

meta[condition_new]=meta[,"condition"]
return(meta)
}



testResults<-function(seur,condition,celltype,basetype,individual,numPerm=100)
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

meta=permuteIndividual(meta,individual,celltype,condition,condition_new)
print(head(meta))
print(i)
dat=tryCatch({
tmp=RunAll(meta,"condition","CellType",basetype,"ID_loc")
},error = function(err){
print("Yuck!")
print(err)
return(data.frame())
})

if(dim(dat)[1]>0)
{
dat["iter"]=i
if(i==1){ret=dat}
else{ret=rbind(ret,dat)}
}
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
ret=testResults(seur,condition,celltype,basetype,individual)
savefil=args[6]
write.table(ret,savefil,row.names=F,sep="\t",quote=F)
}









