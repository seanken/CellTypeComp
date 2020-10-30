library(lme4)
library(dplyr)
library(tidyr)
library (DirichletReg)
library(mclogit)




##Fisher Exact
test_fisher<-function(meta,condition,celltype)
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
test_chi<-function(meta,condition,celltype)
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
test_logistic<-function(meta,condition,celltype)
{
meta["CellType"]=meta[,celltype]
meta["Cond"]=meta[,condition]
celltypes=unique(meta[,celltype])
res=c()
for(i in celltypes)
{
meta["Cur"]=meta[,celltype]==i

fit=glm(Cur~Cond,data=meta,family="binomial")

#print(summary(fit)$coefficients)
pval=summary(fit)$coefficients[grep("^Condwt",rownames(summary(fit)$coefficients)),4]
res=c(res,pval)
}
names(res)=celltypes
return(res)
}


##Direchlet
test_direchlet<-function(meta,condition,celltype,individual,method="common",basetype="")
{
meta=meta[,c(condition,celltype,individual)]
colnames(meta)=c("Cond","CellType","Ind")
tab<-meta %>% group_by(CellType,Cond,Ind) %>% summarise(Count=length(Cond)) %>% as.data.frame()
tab2<-meta %>% group_by(Ind,Cond) %>% summarise(Tot=length(Cond)) %>% as.data.frame()

tab=left_join(tab,tab2)
tab["Percent"]=tab["Count"]/tab["Tot"]
tab["CellType"]=sub("^","CellType_",tab[,"CellType"])

print(levels(tab[,"CellType"]))
if(nchar(basetype)>0)
{
tab["CellType"]=factor(tab[,"CellType"])
tab["CellType"]=relevel(tab[,"CellType"],ref=paste("CellType_",basetype,sep=""))
}

tab=tab[,c("CellType","Ind","Cond","Percent")] %>% spread(CellType,Percent,fill=0)
tab$Y=DR_data(tab[,grep("^CellType",colnames(tab))])
res=DirichReg(Y~Cond,tab,model=method)
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
test_direchlet_base<-function(meta,condition,celltype,individual,basetype)
{
meta[celltype]=factor(meta[,celltype])
meta[celltype]=relevel(meta[,celltype],ref=basetype)
return(test_direchlet(meta,condition,celltype,individual,method="alternative",basetype=basetype))
}


##multinomial model--just calls the multinomial model with no random effects and no overdisplersion
test_multi<-function(meta,condition,celltype,individual,basetype)
{
meta=meta[,c(condition,celltype,individual)]
colnames(meta)=c("Cond","CellType","Ind")
meta["CellType"]=factor(meta[,"CellType"])
meta[,"CellType"]=relevel(meta[,"CellType"],ref=basetype)
res=mblogit(CellType~Cond,data=meta)
mat=summary(res)$coefficients
mat=mat[grep("Cond",rownames(mat)),]
res=data.frame(mat)[,4]
names(res)=as.character(lapply(rownames(mat),function(x){strsplit(x,"~",fixed=T)[[1]][1]}))
return(res)
}


##multinomial mixed model with overdispersion
test_multi_mixed_overdisp<-function(meta,condition,celltype,individual,basetype,overdisp="Afroz")
{
meta=meta[,c(condition,celltype,individual)]
colnames(meta)=c("Cond","CellType","Ind")
meta["CellType"]=factor(meta[,"CellType"])
meta[,"CellType"]=relevel(meta[,"CellType"],ref=basetype)
res=mblogit(CellType~Cond,random=~1|Ind,data=meta,dispersion=overdisp)
mat=summary(res)$coefficients
mat=mat[grep("Cond",rownames(mat)),]
print(summary(res))
res=data.frame(mat)[,4]
names(res)=as.character(lapply(rownames(mat),function(x){strsplit(x,"~",fixed=T)[[1]][1]}))
return(res)
}


test_multi_mixed_table<-function(meta,condition,celltype,individual,basetype,overdisp="Afroz")
{
meta=meta[,c(condition,celltype,individual)]
colnames(meta)=c("Cond","CellType","Ind")
meta["CellType"]=factor(meta[,"CellType"])
meta[,"CellType"]=relevel(meta[,"CellType"],ref=basetype)
for(i in colnames(meta)){meta[i]=factor(meta[,i])}

print(as.data.frame(table(meta)))
meta<-meta %>% group_by(CellType,Cond,Ind) %>% summarise(Freq=length(Cond)) %>% as.data.frame()
print(head(meta))
print(class(meta[,"Freq"]))
res=mblogit(CellType~Cond,random=~1|Ind,data=as.data.frame(table(meta)),dispersion=overdisp,weights = Freq,from.table=T)
mat=summary(res)$coefficients
mat=mat[grep("Cond",rownames(mat)),]
print(summary(res))
res=data.frame(mat)[,4]
names(res)=as.character(lapply(rownames(mat),function(x){strsplit(x,"~",fixed=T)[[1]][1]}))
return(res)
}



##multinomial mixed model with overdispersion v2
test_multi_mixed_overdisp_v2<-function(meta,condition,celltype,individual,basetype)
{
return(test_multi_mixed_overdisp(meta,condition,celltype,individual,basetype,overdisp="Fletcher"))
}

##multinomial mixed model with overdispersion v3
test_multi_mixed_overdisp_v3<-function(meta,condition,celltype,individual,basetype)
{
return(test_multi_mixed_overdisp(meta,condition,celltype,individual,basetype,overdisp="Pearson"))
}

##multinomial mixed model with overdispersion v4
test_multi_mixed_overdisp_v4<-function(meta,condition,celltype,individual,basetype)
{
return(test_multi_mixed_overdisp(meta,condition,celltype,individual,basetype,overdisp="Deviance"))
}


##multinomial mixed model with no overdispersion
test_multi_mixed<-function(meta,condition,celltype,individual,basetype)
{
return(test_multi_mixed_overdisp(meta,condition,celltype,individual,basetype,overdisp=F))
}





##Poisson Regression
test_poisson<-function(meta,condition,celltype,individual)
{
meta=meta[,c(condition,celltype,individual)]
colnames(meta)=c("Cond","CellType","Ind")
tab<-meta %>% group_by(CellType,Ind,Cond) %>% summarise(Count=length(Ind)) %>% as.data.frame()
tab2<-meta %>% group_by(Ind,Cond) %>% summarise(Tot=length(Ind)) %>% as.data.frame()
tab=left_join(tab,tab2)
tab["logTot"]=log(tab[,"Tot"])
celltypes=unique(tab[,"CellType"])

res<-as.numeric(lapply(celltypes,function(x){

cur=tab[tab[,"CellType"]==x,]
fit<-glm(Count~offset(logTot)+Cond,data=cur,family="poisson")
mat=summary(fit)$coefficients
print(mat)
return(mat[grep("^Cond",rownames(mat)),4])
}))
names(res)=celltypes
return(res)

}

##logistic mixed model
test_logistic_mixed<-function(meta,condition,celltype,individual)
{
meta=meta[,c(condition,celltype,individual)]
colnames(meta)=c("Cond","CellType","Ind")
celltypes=unique(meta[,celltype])
res=c()
for(i in celltypes)
{
meta["CellType_cur"]=meta[,celltype]==i

fit=glmer(CellType_cur~Cond+(1|Ind),data=meta,family="binomial")
print(summary(fit)$coefficients)
#print(summary(fit)$coefficients)
pval=summary(fit)$coefficients[grep("^Cond",rownames(summary(fit)$coefficients)),4]
res=c(res,pval)
}
names(res)=celltypes
return(res)
}






##Wilcoxen test with base
test_wilcox_base<-function(meta,condition,celltype,individual,basetype)
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
test_wilcox<-function(meta,condition,celltype,individual)
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




