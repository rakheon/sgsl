#' Fit the sparse group-subgroup lasso (SGSL)
#'
#' @param x p by N matrix of predictors (N: sample size, p: number of predictors)
#' @param y 1 by N matrix of response variable
#' @param type One of "lasso", "group", "ggroup", "ggroupind" or "sgsl".
#' @param index.subgroup index for subgroups
#' @param tau tau
#' @param delta delta for cross-validation
#' @param delta.group delta for cross-validation for groups
#' @param delta.subgroup delta for cross-validation for subgroups
#' @param delta.ind delta for cross-validation for ind
#' @param standardize logical. TRUE for standardizing the data.
#'
#' @return out: indicators of the selected predictors. 1 for selected predictors and 0 for not selected predictors
#' @export
#'
#' @examples
#' tau=0.94
#' set.seed(1)
#' x <- matrix(rnorm(2400), nrow=80)
#' y <- 0.25*x[1,] + 0.5*x[2,] + 0.75*x[3,] + 1.0*x[4,] + 0.25*x[5,] + 0.5*x[6,] + 0.75*x[7,] + 1.0*x[8,]  + matrix(rnorm(30), nrow=1)
#' index.subgroup <- matrix(NA,nrow=10,ncol=80)
#' index.subgroup[1,1:4]=1; index.subgroup[1,5:8]=2
#' index.subgroup[2,9:12]=3; index.subgroup[2,13:16]=4
#' index.subgroup[3,17:20]=5; index.subgroup[3,21:24]=6
#' index.subgroup[4,25:28]=7; index.subgroup[4,29:32]=8
#' index.subgroup[5,33:36]=9; index.subgroup[5,37:40]=10
#' index.subgroup[6,41:44]=11; index.subgroup[6,45:48]=12
#' index.subgroup[7,49:52]=13; index.subgroup[7,53:56]=14
#' index.subgroup[8,57:60]=15; index.subgroup[8,61:64]=16
#' index.subgroup[9,65:68]=17; index.subgroup[9,69:72]=18
#' index.subgroup[10,73:76]=19; index.subgroup[10,77:80]=20
#' sgslout <- sgsl(x,y,type="lasso",index.subgroup = index.subgroup)
sgsl <- function(x,y,type=c("lasso", "group", "ggroup", "ggroupind", "sgsl")[1],
                 index.subgroup,#p.group=NULL,index=NULL,
                 #file.group="file",plots.ind=FALSE,plots.group=FALSE,plots.subgroup=FALSE,
                 tau=NULL,delta=2,delta.group=2,delta.subgroup=2,delta.ind=2,
                 #use.Gram=TRUE,format.data=FALSE,
                 standardize=TRUE){

    tools <- form.tools(index.subgroup)
    index <- tools$index
    p.group <- tools$p.group
    p.subgroup <- tools$p.subgroup
    x <- as.data.frame(x)
    y <- as.data.frame(y)

    ## Methods ##
    if (type == "lasso"){
        ## Apply Lasso ##
        out <- lasso.computations(microbes=x,phenotypes=y,plots=FALSE,file="lasso_",
                                  format.data=FALSE,
                                  delta=delta,use.Gram=TRUE)
    } else if (type == "group"){
        ## Apply pure group Lasso ##
        out <- pure.grplasso.computations(microbes=x,phenotypes=y,
                                          index=index,p.group=p.group,tau=tau,file.group="pure_group_lasso_",
                                          plots.group=FALSE,delta.group=delta.group,
                                          format.data=FALSE,standardize=standardize)
    } else if (type == "ggroup"){
        ## group group Lasso ##
        out <- group.group.lasso.computations(microbes=x,phenotypes=y,index,index.subgroup,p.group,tau=tau,
                                              file.group="group_group_lasso_",plots.group=FALSE,plots.subgroup=FALSE,
                                              delta.group=2,delta.subgroup=2,format.data=FALSE,standardize=standardize)
    } else if (type == "ggroupind"){
        out <- group.group.indlasso.lasso.computations(microbes=x,phenotypes=y,index,index.subgroup,p.group,tau=tau,
                                                       file.group="group_group_indlasso_",plots.group=FALSE,plots.subgroup=FALSE,
                                                       plots.ind=FALSE,delta.group=2,delta.subgroup=2,delta.ind=2,
                                                       use.Gram=TRUE,format.data=FALSE,standardize=standardize)
    } else {
        out <- NULL
        #  out <- sparse.group.subgroup.computations(microbes=x,phenotypes=y,group.index=index,subgroup.index=index.subgroup,
        #                                            tau=tau,alpha1,alpha2,alpha3,nlam,lambdas,lambda.accuracy,
        #                                            file.group="sgsl_",plots.group=FALSE,delta.group=2,format.data=FALSE,
        #                                            cv.criterion=FALSE,nfold=10,alphas.cv.range=seq(0.1,0.95,by=0.05))
    }

    return(out)

}
