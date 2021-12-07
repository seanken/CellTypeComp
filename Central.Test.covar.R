#library(Seurat);
library(dplyr);library(tidyr)
source("RunAll_covar.R")

##
##Permutes the individuals so 
##Here, meta is the metadata from Seurat object
##ID_loc is the column containing individuals ID
##condition is the column containing condition info (case vs control, etc)
##condition_new is the column of perturbed conditions
##
permuteIndividual_covar<-function(meta,ID_loc,CellType,condition,condition_new)
{
meta=meta[,c(ID_loc,condition,CellType)]
colnames(meta)=c("ID_loc","condition","CellType")
tab<-meta %>% group_by(ID_loc,condition) %>% summarise() %>% as.data.frame()
tab[condition_new]=sample(tab[,"condition"])
tab["condition"]=tab[,condition_new]
colnames(meta)[2]="Covar"
#meta=meta[,c("ID_loc","CellType")]
meta=left_join(meta,tab)
meta["Covar"]=factor(meta[,"Covar"])
meta["CellType"]=factor(meta[,"CellType"])
meta["condition"]=factor(meta[,"condition"])
meta["ID_loc"]=factor(meta[,"ID_loc"])
return(meta)
}



testResults_covar<-function(seur,condition,celltype,basetype,individual,numPerm=1000)
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
I=1
for(i in 1:numPerm)
{

if(is.data.frame(seur))
{meta=seur}
else
{
meta=seur@meta.data
}

meta=permuteIndividual_covar(meta,individual,celltype,condition,condition_new)

print(head(meta))
print(" ")
print(" ")
print(" ")
print("New iter!")
print(head(meta))
print(i)
dat=tryCatch({
tmp=RunAll_covar(meta,"condition","CellType",basetype,"ID_loc","Covar")
#return(tmp)
},error = function(err){
print("Yuck!")
return(data.frame())
})
print("continue")
if(dim(dat)[1]>0)
{
print(head(dat))
dat["iter"]=i
if(I==1){ret=dat;I=2}
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









