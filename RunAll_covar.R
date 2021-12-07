source("Ind.Test.with.covar.R")
library(tictoc)

RunAll_covar<-function(meta,condition,celltype,basetype,individual,covar)
{
times=c()
tic();res<-test_covar_logistic(meta,condition,celltype,covar);nams=names(res);pvals=as.numeric(res);dat=data.frame(nams);dat["pval"]=pvals;colnames(dat)[1]="CellType";dat["method"]="logistic";out=toc();dat["time"]=out$toc-out$tic;ret=dat;times<-c(times,out$tic-out$toc);
tic();res<-test_covar_direchlet(meta,condition,celltype,individual,covar,method="common",basetype="");nams=names(res);pvals=as.numeric(res);dat=data.frame(nams);dat["pval"]=pvals;colnames(dat)[1]="CellType";dat["method"]="dirichlet";out=toc();dat["time"]=out$toc-out$tic;ret=rbind(ret,dat);times<-c(times,out$tic-out$toc);
tic();res<-test_covar_direchlet_base(meta,condition,celltype,individual,basetype,covar);nams=names(res);pvals=as.numeric(res);dat=data.frame(nams);dat["pval"]=pvals;colnames(dat)[1]="CellType";dat["method"]="dirich_alt";out=toc();dat["time"]=out$toc-out$tic;ret=rbind(ret,dat);times<-c(times,out$tic-out$toc);
tic();res<-test_covar_multi(meta,condition,celltype,individual,basetype,covar);nams=names(res);pvals=as.numeric(res);dat=data.frame(nams);dat["pval"]=pvals;colnames(dat)[1]="CellType";dat["method"]="multinomial";out=toc();dat["time"]=out$toc-out$tic;ret=rbind(ret,dat);times<-c(times,out$tic-out$toc);
tic();res<-test_covar_multi_mixed_overdisp(meta,condition,celltype,individual,basetype,covar);nams=names(res);pvals=as.numeric(res);dat=data.frame(nams);dat["pval"]=pvals;colnames(dat)[1]="CellType";dat["method"]="multinomial_mixed_overdisp";out=toc();dat["time"]=out$toc-out$tic;ret=rbind(ret,dat);times<-c(times,out$tic-out$toc);
tic();res<-test_covar_multi_mixed(meta,condition,celltype,individual,basetype,covar);nams=names(res);pvals=as.numeric(res);dat=data.frame(nams);dat["pval"]=pvals;colnames(dat)[1]="CellType";dat["method"]="multinomial_mixed";out=toc();dat["time"]=out$toc-out$tic;ret=rbind(ret,dat);times<-c(times,out$tic-out$toc);
tic();res<-test_covar_poisson(meta,condition,celltype,individual,covar);nams=names(res);pvals=as.numeric(res);dat=data.frame(nams);dat["pval"]=pvals;colnames(dat)[1]="CellType";dat["method"]="poisson";out=toc();dat["time"]=out$toc-out$tic;ret=rbind(ret,dat);times<-c(times,out$tic-out$toc);
tic();res<-test_covar_nb(meta,condition,celltype,individual,covar);nams=names(res);pvals=as.numeric(res);dat=data.frame(nams);dat["pval"]=pvals;colnames(dat)[1]="CellType";dat["method"]="nb";out=toc();dat["time"]=out$toc-out$tic;ret=rbind(ret,dat);times<-c(times,out$tic-out$toc);
tic();res<-test_covar_logistic_mixed(meta,condition,celltype,individual,covar);nams=names(res);pvals=as.numeric(res);dat=data.frame(nams);dat["pval"]=pvals;colnames(dat)[1]="CellType";dat["method"]="logistic_mixed";out=toc();dat["time"]=out$toc-out$tic;ret=rbind(ret,dat);times<-c(times,out$tic-out$toc);
tic();res<-test_covar_wilcox_base(meta,condition,celltype,individual,basetype);nams=names(res);pvals=as.numeric(res);dat=data.frame(nams);dat["pval"]=pvals;colnames(dat)[1]="CellType";dat["method"]="wilcox_base";out=toc();dat["time"]=out$toc-out$tic;ret=rbind(ret,dat);times<-c(times,out$tic-out$toc);
tic();res<-test_covar_wilcox(meta,condition,celltype,individual);nams=names(res);pvals=as.numeric(res);dat=data.frame(nams);dat["pval"]=pvals;colnames(dat)[1]="CellType";dat["method"]="wilcox_percent";out=toc();dat["time"]=out$toc-out$tic;ret=rbind(ret,dat);times<-c(times,out$tic-out$toc);
tic();res<-test_covar_chi(meta,condition,celltype);nams=names(res);pvals=as.numeric(res);dat=data.frame(nams);dat["pval"]=pvals;colnames(dat)[1]="CellType";dat["method"]="chi";out=toc();dat["time"]=out$toc-out$tic;ret=rbind(ret,dat);times<-c(times,out$tic-out$toc);
tic();res<-test_covar_fisher(meta,condition,celltype);nams=names(res);pvals=as.numeric(res);dat=data.frame(nams);dat["pval"]=pvals;colnames(dat)[1]="CellType";dat["method"]="fisher";out=toc();dat["time"]=out$toc-out$tic;ret=rbind(ret,dat);times<-c(times,out$tic-out$toc);
tic();res<-test_covar_propel(meta,condition,celltype,individual,covar);nams=names(res);pvals=as.numeric(res);dat=data.frame(nams);dat["pval"]=pvals;colnames(dat)[1]="CellType";dat["method"]="propel_asin";out=toc();dat["time"]=out$toc-out$tic;ret=rbind(ret,dat);times<-c(times,out$tic-out$toc);
tic();res<-test_covar_propel(meta,condition,celltype,individual,covar,"logit");nams=names(res);pvals=as.numeric(res);dat=data.frame(nams);dat["pval"]=pvals;colnames(dat)[1]="CellType";dat["method"]="propel_logit";out=toc();dat["time"]=out$toc-out$tic;ret=rbind(ret,dat);times<-c(times,out$tic-out$toc);


#tic();res<-test_covar_multi_overdisp_v2(meta,condition,celltype,individual,basetype,covar);nams=names(res);pvals=as.numeric(res);dat=data.frame(nams);dat["pval"]=pvals;colnames(dat)[1]="CellType";dat["method"]="multinomial_overdisp_v2";out=toc();dat["time"]=out$toc-out$tic;ret=rbind(ret,dat);times<-c(times,out$tic-out$toc);
#print(15)
#tic();res<-test_covar_multi_overdisp_v3(meta,condition,celltype,individual,basetype,covar);nams=names(res);pvals=as.numeric(res);dat=data.frame(nams);dat["pval"]=pvals;colnames(dat)[1]="CellType";dat["method"]="multinomial_overdisp_v3";out=toc();dat["time"]=out$toc-out$tic;ret=rbind(ret,dat);times<-c(times,out$tic-out$toc);
#print(16)
#tic();res<-test_covar_multi_overdisp_v4(meta,condition,celltype,individual,basetype,covar);nams=names(res);pvals=as.numeric(res);dat=data.frame(nams);dat["pval"]=pvals;colnames(dat)[1]="CellType";dat["method"]="multinomial_overdisp_v4";out=toc();dat["time"]=out$toc-out$tic;ret=rbind(ret,dat);times<-c(times,out$tic-out$toc);
#print(17)
#tic();res<-test_covar_multi_overdisp_v1(meta,condition,celltype,individual,basetype,covar);nams=names(res);pvals=as.numeric(res);dat=data.frame(nams);dat["pval"]=pvals;colnames(dat)[1]="CellType";dat["method"]="multinomial_overdisp_v1";out=toc();dat["time"]=out$toc-out$tic;ret=rbind(ret,dat);times<-c(times,out$tic-out$toc);
return(ret)
}
