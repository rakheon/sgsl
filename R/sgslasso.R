sgsl <- function(x,y,type=c("lasso", "group", "ggroup", "ggroupind", "sgsl")[1],
                 index=NULL,index.subgroup=NULL,p.group=NULL,tau=NULL,
                 file.group="file",plots.group=FALSE,plots.subgroup=FALSE,
                 plots.ind=FALSE,delta.group=2,delta.subgroup=2,delta.ind=2,
                 use.Gram=TRUE,format.data=FALSE,standardize=TRUE){

    ## Methods ##
    if (type == "lasso"){
        ## Apply Lasso ##
        out <- lasso.computations(microbes=x,phenotypes=y,plots=FALSE,file="lasso_",
                                  format.data=FALSE,
                                  delta=delta,use.Gram=use.Gram)$interest
    } else if (type == "group"){
        ## Apply pure group Lasso ##
        out <- pure.grplasso.computations(microbes=x,phenotypes=y,
                                          index=index,p.group=p.group,tau=tau,file.group="pure_group_lasso_",
                                          plots.group=plots.group,delta.group=delta,
                                          format.data=FALSE,standardize=standardize)$interest
    } else if (type == "ggroup"){
        ## group group Lasso ##
        out <- group.group.lasso.computations(microbes=x,phenotypes=y,index,index.subgroup,p.group,tau=tau,
                                              file.group="group_group_lasso_",plots.group=plots.group,plots.subgroup=plots.subgroup,
                                              delta.group=2,delta.subgroup=2,format.data=FALSE,standardize=standardize)
    } else if (type == "ggroupind"){
        out <- group.group.indlasso.lasso.computations(microbes=x,phenotypes=y,index,index.subgroup,p.group,tau=tau,
                                                       file.group="group_group_indlasso_",plots.group=plots.group,plots.subgroup=plots.subgroup,
                                                       plots.ind=plots.ind,delta.group=2,delta.subgroup=2,delta.ind=2,
                                                       use.Gram=TRUE,format.data=FALSE,standardize=standardize)
    } else if (type == "sgsl"){
        out <- sparse.group.subgroup.computations(microbes=x,phenotypes=y,group.index=index,subgroup.index=index.subgroup,
                                                  tau=tau,alpha1,alpha2,alpha3,nlam,lambdas,lambda.accuracy,
                                                  file.group="sgsl_",plots.group=plots.group,delta.group=2,format.data=FALSE,
                                                  cv.criterion=FALSE,nfold=10,alphas.cv.range=seq(0.1,0.95,by=0.05))
    }




    return(out)

}
