library(MASS)
library(lme4)
library(dplyr)
library(tidyr)
library (DirichletReg)
library(mclogit)
source("/stanley/levin_dr/ssimmons/TestCellComp/Code_Expanded/Propel/propel.core.modified.with.covar.R")



##Fisher Exact
test_covar_fisher<-function(meta,condition,celltype)
{
meta["Cond"]=factor(meta[,condition])
celltypes=unique(meta[,celltype])

res=c()
for(i in celltypes)
{
meta["Cur"]=factor(meta[,celltype]==i)
pval=fisher.test(meta[,"Cur"],meta[,"Cond"])$p.val
res<-c(res,pval)
}
names(res)=celltypes
return(res)
}

##Test Chi Squared
test_covar_chi<-function(meta,condition,celltype)
{
meta["Cond"]=factor(meta[,condition])
res=c()
celltypes=unique(meta[,celltype])
for(i in celltypes)
{
meta["Cur"]=meta[,"CellType"]==i
tab<-meta %>% group_by(Cur,Cond) %>% summarise(Num=length(Cond)) %>% spread(Cond,Num,fill=0) %>% as.data.frame()
tab=tab[,2:3]
pval=chisq.test(tab)$p.value
res<-c(res,pval)
}
names(res)=celltypes
return(res)

}


##Logistic Regression
test_covar_logistic<-function(meta,condition,celltype,covar)
{
meta["CellType"]=meta[,celltype]
meta["Cond"]=meta[,condition]
meta["Covar"]=meta[,covar]
celltypes=unique(meta[,celltype])
res=c()
for(i in celltypes)
{
meta["Cur"]=meta[,celltype]==i

fit=glm(Cur~Cond+Covar,data=meta,family="binomial")

#print(summary(fit)$coefficients)
pval=summary(fit)$coefficients[grep("^Condwt",rownames(summary(fit)$coefficients)),4]
res=c(res,pval)
}
names(res)=celltypes
return(res)
}


##Direchlet
test_covar_direchlet<-function(meta,condition,celltype,individual,covar,method="common",basetype="")
{
meta=meta[,c(condition,celltype,individual,covar)]
colnames(meta)=c("Cond","CellType","Ind","Covar")
tab<-meta %>% group_by(CellType,Cond,Ind,Covar) %>% summarise(Count=length(Cond)) %>% as.data.frame()
tab2<-meta %>% group_by(Ind,Cond,Covar) %>% summarise(Tot=length(Cond)) %>% as.data.frame()

tab=left_join(tab,tab2)
tab["Percent"]=tab["Count"]/tab["Tot"]
tab["CellType"]=sub("^","CellType_",tab[,"CellType"])

print(levels(tab[,"CellType"]))
if(nchar(basetype)>0)
{
tab["CellType"]=factor(tab[,"CellType"])
tab["CellType"]=relevel(tab[,"CellType"],ref=paste("CellType_",basetype,sep=""))
}

tab=tab[,c("CellType","Ind","Cond","Covar","Percent")] %>% spread(CellType,Percent,fill=0)
tab$Y=DR_data(tab[,grep("^CellType",colnames(tab))])
res=DirichReg(Y~Cond+Covar,tab,model=method)
mat=summary(res)$coef.mat
nams=sub("^CellType_","",summary(res)$varnames)
mat=mat[grep("^Cond",rownames(mat)),]
res=data.frame(mat)[,4]
if(method=="common")
{
names(res)=nams
}
else{names(res)=nams[2:length(nams)]}
return(res)
}


##Direchlet Alternative
test_covar_direchlet_base<-function(meta,condition,celltype,individual,basetype,covar)
{
meta[celltype]=factor(meta[,celltype])
meta[celltype]=relevel(meta[,celltype],ref=basetype)
return(test_covar_direchlet(meta,condition,celltype,individual,covar,method="alternative",basetype=basetype))
}


##multinomial model--just calls the multinomial model with no random effects and no overdisplersion
test_covar_multi<-function(meta,condition,celltype,individual,basetype,covar)
{
meta=meta[,c(condition,celltype,individual,covar)]
colnames(meta)=c("Cond","CellType","Ind","Covar")
meta["CellType"]=factor(meta[,"CellType"])
meta[,"CellType"]=relevel(meta[,"CellType"],ref=basetype)
res=mblogit(CellType~Cond+Covar,data=meta)
mat=summary(res)$coefficients
mat=mat[grep("Cond",rownames(mat)),]
res=data.frame(mat)[,4]
names(res)=as.character(lapply(rownames(mat),function(x){strsplit(x,"~",fixed=T)[[1]][1]}))
return(res)
}


##multinomial mixed model with overdispersion
test_covar_multi_mixed_overdisp<-function(meta,condition,celltype,individual,basetype,covar,overdisp="Afroz")
{
meta=meta[,c(condition,celltype,individual,covar)]
colnames(meta)=c("Cond","CellType","Ind","Covar")
meta["CellType"]=factor(meta[,"CellType"])
meta[,"CellType"]=relevel(meta[,"CellType"],ref=basetype)
res=mblogit(CellType~Cond+Covar,random=~1|Ind,data=meta,dispersion=overdisp)
mat=summary(res)$coefficients
mat=mat[grep("Cond",rownames(mat)),]
print(summary(res))
res=data.frame(mat)[,4]
names(res)=as.character(lapply(rownames(mat),function(x){strsplit(x,"~",fixed=T)[[1]][1]}))
return(res)
}





##multinomial mixed model with overdispersion v2
test_covar_multi_mixed_overdisp_v2<-function(meta,condition,celltype,individual,basetype,covar)
{
return(test_covar_multi_mixed_overdisp(meta,condition,celltype,individual,basetype,covar,overdisp="Fletcher"))
}

##multinomial mixed model with overdispersion v3
test_covar_multi_mixed_overdisp_v3<-function(meta,condition,celltype,individual,basetype,covar)
{
return(test_covar_multi_mixed_overdisp(meta,condition,celltype,individual,basetype,covar,overdisp="Pearson"))
}

##multinomial mixed model with overdispersion v4
test_covar_multi_mixed_overdisp_v4<-function(meta,condition,celltype,individual,basetype,covar)
{
return(test_covar_multi_mixed_overdisp(meta,condition,celltype,individual,basetype,covar,overdisp="Deviance"))
}


##multinomial mixed model with no overdispersion
test_covar_multi_mixed<-function(meta,condition,celltype,individual,basetype,covar)
{
return(test_covar_multi_mixed_overdisp(meta,condition,celltype,individual,basetype,covar,overdisp=F))
}




##Poisson Regression
test_covar_nb<-function(meta,condition,celltype,individual,covar)
{
meta=meta[,c(condition,celltype,individual,covar)]
colnames(meta)=c("Cond","CellType","Ind","Covar")
tab<-meta %>% group_by(CellType,Ind,Cond,Covar) %>% summarise(Count=length(Ind)) %>% as.data.frame()
tab2<-meta %>% group_by(Ind,Cond,Covar) %>% summarise(Tot=length(Ind)) %>% as.data.frame()
tab=left_join(tab,tab2)
tab["logTot"]=log(tab[,"Tot"])
celltypes=unique(tab[,"CellType"])

res<-as.numeric(lapply(celltypes,function(x){

cur=tab[tab[,"CellType"]==x,]
fit<-glm.nb(Count~offset(logTot)+Cond+Covar,data=cur)
mat=summary(fit)$coefficients
print(mat)
return(mat[grep("^Cond",rownames(mat)),4])
}))
names(res)=celltypes
return(res)

}


##Poisson Regression
test_covar_poisson<-function(meta,condition,celltype,individual,covar)
{
meta=meta[,c(condition,celltype,individual,covar)]
colnames(meta)=c("Cond","CellType","Ind","Covar")
tab<-meta %>% group_by(CellType,Ind,Cond,Covar) %>% summarise(Count=length(Ind)) %>% as.data.frame()
tab2<-meta %>% group_by(Ind,Cond,Covar) %>% summarise(Tot=length(Ind)) %>% as.data.frame()
tab=left_join(tab,tab2)
tab["logTot"]=log(tab[,"Tot"])
celltypes=unique(tab[,"CellType"])

res<-as.numeric(lapply(celltypes,function(x){

cur=tab[tab[,"CellType"]==x,]
fit<-glm(Count~offset(logTot)+Cond+Covar,data=cur,family="poisson")
mat=summary(fit)$coefficients
print(mat)
return(mat[grep("^Cond",rownames(mat)),4])
}))
names(res)=celltypes
return(res)

}

##logistic mixed model
test_covar_logistic_mixed<-function(meta,condition,celltype,individual,covar)
{
meta=meta[,c(condition,celltype,individual,covar)]
colnames(meta)=c("Cond","CellType","Ind","Covar")
celltypes=unique(meta[,celltype])
res=c()
for(i in celltypes)
{
meta["CellType_cur"]=meta[,celltype]==i

fit=glmer(CellType_cur~Cond+Covar+(1|Ind),data=meta,family="binomial")
print(summary(fit)$coefficients)
#print(summary(fit)$coefficients)
pval=summary(fit)$coefficients[grep("^Cond",rownames(summary(fit)$coefficients)),4]
res=c(res,pval)
}
names(res)=celltypes
return(res)
}






##Wilcoxen test with base
test_covar_wilcox_base<-function(meta,condition,celltype,individual,basetype)
{
celltypes=unique(meta[,celltype])
celltypes=celltypes[celltypes!=basetype]
meta=meta[,c(condition,celltype,individual)]
colnames(meta)=c("cond","celltype","indiv")
tab<-meta %>% group_by(cond,celltype,indiv) %>% summarise(Num=length(cond)) %>% as.data.frame()
tab2<-tab %>% group_by(cond,indiv) %>% summarise(base=Num[celltype==basetype]) %>% as.data.frame()
tab<-left_join(tab,tab2)
tab["Rat"]=tab[,"Num"]/tab[,"base"]
res<-lapply(celltypes,
function(x){
cur=tab[tab[,"celltype"]==x,]
conds=unique(tab[,"cond"])
pval=wilcox.test(cur[cur[,"cond"]==conds[1],"Rat"],cur[cur[,"cond"]==conds[2],"Rat"])$p.value
return(pval)
}
)

res=as.numeric(res)
names(res)=celltypes

return(res)

}


##Wilcoxen test without base
test_covar_wilcox<-function(meta,condition,celltype,individual)
{
celltypes=unique(meta[,celltype])
meta=meta[,c(condition,celltype,individual)]
colnames(meta)=c("cond","celltype","indiv")
tab<-meta %>% group_by(cond,celltype,indiv) %>% summarise(Num=length(cond)) %>% as.data.frame()
tab2<-tab %>% group_by(cond,indiv) %>% summarise(base=sum(Num)) %>% as.data.frame()
tab<-left_join(tab,tab2)
tab["Rat"]=tab[,"Num"]/tab[,"base"]
res<-lapply(celltypes,
function(x){
cur=tab[tab[,"celltype"]==x,]
conds=unique(tab[,"cond"])
pval=wilcox.test(cur[cur[,"cond"]==conds[1],"Rat"],cur[cur[,"cond"]==conds[2],"Rat"])$p.value
return(pval)
}
)
res=as.numeric(res)
names(res)=celltypes

return(res)

}




##multiple regression on cell # with scaling by base


test_covar_propel<-function(meta,condition,celltype,individual,covar,transform="asin")
{
out=propeller(clusters=meta[,celltype],sample=meta[,individual],group=meta[,condition],transform=transform,covar=meta[,covar])
pvals=out[,"P.Value"]
names(pvals)=out[,1]
return(pvals)

}

