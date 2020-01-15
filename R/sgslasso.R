sgslasso <- function(x,y,type=c("lasso")[1],
                     index=NULL,index.subgroup=NULL,p.group=NULL,tau=NULL,
                     file.group="file",plots.group=TRUE,plots.subgroup=TRUE,
                     plots.ind=TRUE,delta.group=2,delta.subgroup=2,delta.ind=2,
                     use.Gram=TRUE,format.data=TRUE,standardize=TRUE){

    if (type == "lasso"){
        out <- lasso.computations(microbes=x,phenotypes=y,plots=FALSE,file="lasso_",
                                  format.data=FALSE,
                                  delta=delta,use.Gram=use.Gram)$interest
    }

    return(out)

}
