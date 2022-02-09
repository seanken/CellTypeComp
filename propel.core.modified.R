##
##This code was taken from https://github.com/Oshlack/speckle and modified to work with our pipeline.
##

library(Matrix)
library(limma)


getTransformedProps <- function(clusters=clusters, sample=sample, 
                                transform=NULL)
{
  if(is.null(transform)) transform <- "asin"
  
  tab <- table(sample, clusters)
  props <- tab/rowSums(tab)
  if(transform=="asin"){
      message("Performing arcsin square root transformation of proportions")
      prop.trans <- asin(sqrt(props))
  }
  else if(transform=="logit"){
      message("Performing logit transformation of proportions")
      props.pseudo <- (tab+0.5)/rowSums(tab+0.5)
      prop.trans <- log(props.pseudo/(1-props.pseudo))
  }
  list(Counts=t(tab), TransformedProps=t(prop.trans), Proportions=t(props))
}
propeller.anova <- function(prop.list=prop.list, design=design, coef = coef,
                            robust=robust, trend=trend, sort=sort)
{
    prop.trans <- prop.list$TransformedProps
    prop <- prop.list$Proportions

    # get cell type mean proportions ignoring other variables
    # this assumes that the design matrix is not in Intercept format
    fit.prop <- lmFit(prop, design[,coef])

    # Change design matrix to intercept format
    design[,1] <- 1
    colnames(design)[1] <- "Int"

    # Fit linear model taking into account all confounding variables
    fit <- lmFit(prop.trans,design)

    # Get F statistics corresponding to group information only
    # You have to remove the intercept term for this to work
    fit <- eBayes(fit[,coef[-1]], robust=robust, trend=trend)

    # Extract F p-value
    p.value <- fit$F.p.value
    # and perform FDR adjustment
    fdr <- p.adjust(fit$F.p.value, method="BH")

    out <- data.frame(PropMean=fit.prop$coefficients, Fstatistic= fit$F,
                    P.Value=p.value, FDR=fdr)
    if(sort){
        o <- order(out$P.Value)
        out[o,]
    }
    else out
}
propeller.ttest <- function(prop.list=prop.list, design=design,
                            contrasts=contrasts, robust=robust, trend=trend,
                            sort=sort)
{
    prop.trans <- prop.list$TransformedProps
    prop <- prop.list$Proportions

    fit <- lmFit(prop.trans, design)
    fit.cont <- contrasts.fit(fit, contrasts=contrasts)
    fit.cont <- eBayes(fit.cont, robust=robust, trend=trend)

    # Get mean cell type proportions and relative risk for output
    # If no confounding variable included in design matrix
    if(length(contrasts)==2){
        fit.prop <- lmFit(prop, design)
        z <- apply(fit.prop$coefficients, 1, function(x) x^contrasts)
        RR <- apply(z, 2, prod)
    }
    # If confounding variables included in design matrix exclude them
    else{
        new.des <- design[,contrasts!=0]
        fit.prop <- lmFit(prop,new.des)
        new.cont <- contrasts[contrasts!=0]
        z <- apply(fit.prop$coefficients, 1, function(x) x^new.cont)
        RR <- apply(z, 2, prod)
    }

    fdr <- p.adjust(fit.cont$p.value, method="BH")

    out <- data.frame(PropMean=fit.prop$coefficients, PropRatio=RR,
                    Tstatistic=fit.cont$t[,1], P.Value=fit.cont$p.value[,1],
                    FDR=fdr)
    if(sort){
        o <- order(out$P.Value)
        out[o,]
    }
    else out
}
propeller <- function(x=NULL, clusters=NULL, sample=NULL, group=NULL,
                      trend=FALSE, robust=TRUE, transform="asin")
{

    
    if(is.null(transform)) transform <- "asin"

    # Get transformed proportions
    prop.list <- getTransformedProps(clusters, sample, transform)

    # Calculate baseline proportions for each cluster
    baseline.props <- table(clusters)/sum(table(clusters))

    # Collapse group information
    group.coll <- table(sample, group)

    design <- matrix(as.integer(group.coll != 0), ncol=ncol(group.coll))
    colnames(design) <- colnames(group.coll)

    if(ncol(design)==2){
        message("group variable has 2 levels, t-tests will be performed")
        contrasts <- c(1,-1)
        out <- propeller.ttest(prop.list, design, contrasts=contrasts,
                               robust=robust, trend=trend, sort=FALSE)
        out <- data.frame(BaselineProp=baseline.props,out)
        o <- order(out$P.Value)
        out[o,]
    }
    else if(ncol(design)>=2){
        message("group variable has > 2 levels, ANOVA will be performed")
        coef <- seq_len(ncol(design))
        out <- propeller.anova(prop.list, design, coef=coef, robust=robust,
                               trend=trend, sort=FALSE)
        out <- data.frame(BaselineProp=as.vector(baseline.props),out)
        o <- order(out$P.Value)
        out[o,]
  }

}




