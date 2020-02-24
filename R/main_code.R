###########
# Library #
###########

library(lars)		## for Lasso
library(grplasso)	## for Group Lasso
library(SGL)		## for Sparse Group Lasso
library(xtable)         ## for printing tables
library(abind)          ## for combining results in realdata_analysis_main.R
library(Matrix)         ## creates block diagonal matrix
library(MASS)           ## generate from multivariate normals

#################
## Source code ##
#################

#source("../main/subgroup_SGL.R")  ## for sparse group-subgroup lasso
#source("../main/mySGLcode.R")     ## SGL code


####################################
## Function for organizing output ##
####################################

## function to compute mean and alpha and (1-alpha) quantiles
summarize.function <- function(x){
  nsimu <- ncol(x)

  lower.bound <- function(y){
    alpha <- 0.05
    sort.y <- sort(y)
    lo.bound <- sort.y[ceiling(length(y) * alpha/2)]
    return(lo.bound)
  }

  upper.bound <- function(y){
    alpha <- 0.05
    sort.y <- sort(y)
    hi.bound <- sort.y[ceiling(length(y) * (1-alpha/2))]
    return(hi.bound)
  }

  mean.x <- apply(x,1,mean)
  low.x <- apply(x,1,lower.bound)
  hi.x <- apply(x,1,upper.bound)

  sd.x <- apply(x,1,sd)
  qnorm.lo <- mean.x - qnorm(1-alpha/2) * sd.x / sqrt(nsimu)
  qnorm.hi <- mean.x + qnorm(1-alpha/2) * sd.x / sqrt(nsimu)

  out.summary <- data.frame(cbind(mean.x,low.x,hi.x,sd.x,qnorm.lo,qnorm.hi))
  colnames(out.summary) <- c("mean","lo_CI","hi_CI","sd","qnorm_lo_CI","qnorm_hi_CI")
  rownames(out.summary) <- rownames(x)

  return(out.summary)
}



###################################
## Functions for organizing data ##
###################################

## function to count number of non-NA elements in a vector
number.non.na <- function(x){
  return(sum(!is.na(x)))
}



## function to form p.group,p.subgroup, and index variables
## p.group : number of covariates in each group
## p.subgroup : matrix of number of covariates in each subgroup (each row corresponds to a group)
## index : group assignments
form.tools <- function(index.subgroup){
  p.group <- apply(index.subgroup,1,number.non.na)
  p.subgroup <- as.numeric(table(index.subgroup))
  index <- rep(1:length(p.group),p.group)
  list(p.group=p.group,p.subgroup=p.subgroup,index=index)
}





####################
## LASSO Approach ##
####################



## function to standardize a vector x
make.std <- function(x){
  N <- length(x)
  ( x-mean(x) ) / ( sd(as.vector(x)) * sqrt( N / (N-1) ) )
}

## function to center a vector x
make.center <- function(x){
  return(x-mean(x))
}

lasso <- function(phenotypes,microbes,file="file",plots=TRUE,delta=2,std.y=TRUE,use.Gram=TRUE){
  y <- t(phenotypes)
  X <- t(microbes)

  N <- length(y)

  if(std.y==TRUE){
    ## Standardize y
    y1 <- make.std(y)
  } else {
    ## Centered response
    y1 <-  make.center(y)
  }


  ## Standardized design matrix X
  X1 <- apply(X,2,make.std)

  ## Run Lasso
  Lasso.out <-  lars(X1, y1, type = c("lasso"),
                     trace = FALSE, normalize = FALSE, intercept = FALSE,
                     use.Gram=use.Gram)

  order.variables <- 0

  ## use Cp-like criterion to find best descriptive model
  p = dim(X1)[2]
  s = length(Lasso.out$df)
  p.pos = NULL

  RSS = NULL
  for (i in 1:s){
    RSS[i] = sum((y1-predict(Lasso.out, X1, s=i, type = c("fit"))$fit)**2)
    p.pre = predict(Lasso.out, X1, s=i, type = c("coefficients"))$coefficients
    p.pos = c(p.pos,length(p.pre[abs(p.pre)>1e-10]))
  }

  # Estimated MSE
  MSE <- sd(as.vector(y1)) * sqrt( N / (N-1) )
  MSE <- MSE^2

  p.min = which.min(RSS/MSE+delta*p.pos)

  ## final best descriptive model
  predict.out <- predict(Lasso.out, X1, s=p.min, type = c("coefficients"))

  ind <- which(abs(predict.out$coefficients)>1e-10)
  sig.variables <- rep(0,nrow(microbes))
  sig.variables[ind] <- 1

  if(plots==TRUE){
    postscript(paste(file,"_lasso1.eps",sep=""))
    plot(Lasso.out,cex.axis=1.5,cex.lab=1.5)
    dev.off()

    postscript(paste(file,"_lasso2.eps",sep=""))
    par(mar=c(5, 4, 4, 2)+1)
    plot(1:s,RSS+2*(p.pos),type="l",cex.axis=1.5,cex.lab=1.5,ylab=substitute(M[n](that,p),list(that=delta)), xlab="Steps")

    abline(v=p.min,lty=2)
    dev.off()
  }
  list(order.variables=order.variables,sig.variables=sig.variables)
}

lasso.computations <- function(microbes,phenotypes,plots=TRUE,file="name",format.data = TRUE,delta=2,
			use.Gram=TRUE){
  if(format.data==TRUE){
    phenotypes <- phenotypes[-c(1,2),]
  } else {
    phenotypes <- phenotypes
  }
  interest <- matrix(0,nrow=nrow(microbes),ncol=nrow(phenotypes))
  interest <- as.data.frame(interest)
  rownames(interest) <- rownames(microbes)
  colnames(interest) <- rownames(phenotypes)

  order.var <- array(0,dim=c(nrow(phenotypes),500))

  for(i in 1:ncol(interest)){
    ##print(i)
    lasso.out <- lasso(phenotypes[i,],microbes,file=paste(file,rownames(phenotypes)[i],sep=""),
                       plots=plots,delta=delta,use.Gram=use.Gram)
    interest[,i] <- lasso.out$sig.variables
    order.var[i,1:length(lasso.out$order.variables)] <- lasso.out$order.variables
  }
  list(interest=interest,order.var=order.var)
}




#########################
## GROUP LASSO Approach #
#########################

group.lasso <- function(phenotypes,microbes,index,tau,file="file",plots=TRUE,delta=2,std.y=TRUE,
                        standardize=TRUE){
  y <- t(phenotypes)
  X <- t(microbes)

  N <- length(y)

  if(std.y==TRUE){
    ## Standardize y
    y1 <- make.std(y)
  } else {
    ## Centered response
    y1 <-  make.center(y)
  }


  ## Standardized design matrix X
  X1 <- apply(X,2,make.std)

  ## Use a multiplicative grid for penalty parameter lambda, starting at maximal lambda value
  lambda <- lambdamax(X1,y=y1, index=index,penscale=sqrt,model=LinReg(),
                      center=FALSE,
                      standardize=standardize)  * c(tau^(0:100),0)

  ## Run Group Lasso
  Lasso.out <-  grplasso(X1, y = y1, index = index, lambda = lambda,
                         model = LinReg(),
                         penscale = sqrt,
                         control = grpl.control(update.hess = "lambda",
                           trace = 0),center=FALSE,
                         standardize=standardize)

  ## use Cp-like criterion to find best descriptive model
  p = dim(X1)[2]
  s = length(Lasso.out$lambda)
  p.pos = NULL

  RSS = NULL
  for (i in 1:s){
    RSS[i] = sum((y1-Lasso.out$fitted[,i])**2)
    p.pre = Lasso.out$coefficients[,i]
    p.pos = c(p.pos,length(p.pre[abs(p.pre)>1e-10]))
  }

  ## Estimated MSE
  MSE <- sd(as.vector(y1)) * sqrt( N / (N-1) )
  MSE <- MSE^2


  p.min = which.min(RSS/MSE+delta*p.pos)

  ## final best descriptive model
  predict.out <- Lasso.out$coefficients[,p.min]

  ind <- which(abs(predict.out)>1e-10)
  sig.variables <- rep(0,nrow(microbes))
  sig.variables[ind] <- 1

  if(plots==TRUE){
    postscript(paste(file,"_group_lasso1.eps",sep=""))
    plot(Lasso.out,cex.axis=1.5,cex.lab=1.5)
    dev.off()

    postscript(paste(file,"_group_lasso2.eps",sep=""))
    par(mar=c(5, 4, 4, 2)+1)
    plot(1:s,RSS+2*(p.pos),type="l",cex.axis=1.5,cex.lab=1.5,ylab=substitute(M[n](that,p),list(that=delta)), xlab="lambda")

    abline(v=p.min,lty=2)
    dev.off()
  }
  list(sig.variables=sig.variables)
}

##############################
## PURE GROUP LASSO Approach #
##############################

pure.grplasso.computations <- function(microbes,phenotypes,index,p.group,tau,
                                       file.group="file",plots.group=TRUE,delta.group=2,format.data=TRUE,
                                       standardize=TRUE){
  if(format.data==TRUE){
    phenotypes <- phenotypes[-c(1,2),]
  } else {
    phenotypes <- phenotypes
  }
  interest <- matrix(0,nrow=nrow(microbes),ncol=nrow(phenotypes))
  interest <- as.data.frame(interest)
  rownames(interest) <- rownames(microbes)
  colnames(interest) <- rownames(phenotypes)

  for(i in 1:ncol(interest)){
    ##print(i)
    lasso.out <- group.lasso(phenotypes[i,],microbes,index,tau,
                             file=file.group,plots=plots.group,delta=delta.group,standardize=standardize)

    interest[,i] <- lasso.out$sig.variables
  }
  list(interest=interest)
}



############################################
## PURE GROUP LASSO & GROUP LASSO Approach #
############################################

## delta.group : delta applied to C_p criterion for group lasso
## delta.subgroup : delta applied to C_p critierian for group lasso among subgroups

group.group.lasso <- function(phenotypes,microbes,index,index.subgroup,p.group,tau,
                              file.group="file",plots.group=TRUE,plots.subgroup=TRUE,delta.group=2,delta.subgroup=2,
                              standardize=TRUE){

  ## group lasso to groups
  group.lasso.out <- group.lasso(phenotypes,microbes,index,tau,
                                 file=file.group,plots=plots.group,
                                 delta=delta.group,standardize=standardize)

  sig.variables <- group.lasso.out$sig.variables

  main.ind <- which(sig.variables==1)

  ## setting up to apply group Lasso to subgroups
  tmp.subgroup <- as.vector(t(index.subgroup))
  tmp.subgroup <- tmp.subgroup[!is.na(tmp.subgroup)]
  new.index.subgroup <- tmp.subgroup[main.ind]

  if(sum(main.ind)!=0){
    X.tmp <- microbes[main.ind,]

    ##if(length(main.ind)>1){
    ##   X.tmp <- microbes[main.ind,]
    ##} else {
    ##  X.tmp <- t(data.frame(microbes[main.ind,]))
    ##}

    subgroup.lasso.out <- group.lasso(phenotypes,X.tmp,new.index.subgroup,tau,
                                      file=paste(file.group,"_group_group_lasso_no_",k,sep=""),
                                      plots=plots.subgroup,
                                      delta=delta.subgroup,standardize=standardize)
    sig.variables[main.ind] <- subgroup.lasso.out$sig.variables
  }

  list(sig.variables=sig.variables)

}

group.group.lasso.computations <- function(microbes,phenotypes,index,index.subgroup,p.group,tau,
                                           file.group="file",plots.group=TRUE,plots.subgroup=TRUE,delta.group=2,
                                           delta.subgroup=2,format.data=TRUE,standardize=TRUE){

  if(format.data==TRUE){
    phenotypes <- phenotypes[-c(1,2),]
  } else {
    phenotypes <- phenotypes
  }
  interest <- matrix(0,nrow=nrow(microbes),ncol=nrow(phenotypes))
  interest <- as.data.frame(interest)
  rownames(interest) <- rownames(microbes)
  colnames(interest) <- rownames(phenotypes)

  for(i in 1:ncol(interest)){
    ##print(i)
    lasso.out <- group.group.lasso(phenotypes[i,],microbes,index,index.subgroup,p.group,tau,
                                   file.group=paste(file,rownames(phenotypes)[i],sep=""),
                                   plots.group,plots.subgroup,delta.group,delta.subgroup,
                                   standardize=standardize)

    interest[,i] <- lasso.out$sig.variables
  }
  list(interest=interest)
}






######################################################
## PURE GROUP LASSO, GROUP LASSO, and LASSO Approach #
######################################################
## delta.group : delta applied to C_p criterion for group lasso
## delta.subgroup : delta applied to C_p critierian for group lasso among subgroups
## delta.ind    : delta applied to C_p criterion for lasso with individual features

group.group.indlasso.lasso <- function(phenotypes,microbes,index,index.subgroup,p.group,tau,
                                       file.group="file",plots.group=TRUE,plots.subgroup=TRUE,plots.ind=TRUE,
                                       delta.group=2,delta.subgroup=2,delta.ind=2,use.Gram=TRUE,
                                       standardize=TRUE){

  ## group lasso to groups
  group.lasso.out <- group.lasso(phenotypes,microbes,index,tau,
                                 file=file.group,plots=plots.group,
                                 delta=delta.group,standardize=standardize)

  sig.variables <- group.lasso.out$sig.variables

  main.ind <- which(sig.variables==1)

  ################################################
  ## setting up to apply group Lasso to subgroups #
  #################################################

  tmp.subgroup <- as.vector(t(index.subgroup))
  tmp.subgroup <- tmp.subgroup[!is.na(tmp.subgroup)]
  new.index.subgroup <- tmp.subgroup[main.ind]

  if(sum(main.ind)!=0){
    X.tmp <- microbes[main.ind,]
    ##if(length(main.ind)>1){
    ##  X.tmp <- microbes[main.ind,]
    ##} else {
    ##   X.tmp <- t(data.frame(microbes[main.ind,]))
    ##}
    subgroup.lasso.out <- group.lasso(phenotypes,X.tmp,new.index.subgroup,tau,
                                      file=paste(file.group,"_group_group_lasso_no_",k,sep=""),
                                      plots=plots.subgroup,
                                      delta=delta.subgroup,standardize=standardize)
    sig.variables[main.ind] <- subgroup.lasso.out$sig.variables
  }

######################################################
  ## setting up to apply lasso among variables selected #
#######################################################

  new.main.ind <- which(sig.variables==1)

  if(sum(new.main.ind)!=0){
    X.tmp <- microbes[new.main.ind,]

    ##if(length(new.main.ind)>1){
    ##  X.tmp <- microbes[new.main.ind,]
    ##} else {
    ##   X.tmp <- t(data.frame(microbes[new.main.ind,]))
    ##}
    lasso.out <-  lasso(phenotypes,X.tmp,
                        file=paste(file.group,"_group_group_indlasso_lasso_no_",k,sep=""),plots=plots.ind,
                        delta=delta.ind,use.Gram=use.Gram)
    sig.variables[new.main.ind] <- lasso.out$sig.variables
  }

  list(sig.variables=sig.variables)

}

group.group.indlasso.lasso.computations <- function(microbes,phenotypes,index,index.subgroup,p.group,tau,
                                                    file.group="file",plots.group=TRUE,plots.subgroup=TRUE,
                                                    plots.ind=TRUE,delta.group=2,delta.subgroup=2,delta.ind=2,
                                                    use.Gram=TRUE,format.data=TRUE,standardize=TRUE){
  if(format.data==TRUE){
    phenotypes <- phenotypes[-c(1,2),]
  } else {
    phenotypes <- phenotypes
  }
  interest <- matrix(0,nrow=nrow(microbes),ncol=nrow(phenotypes))
  interest <- as.data.frame(interest)
  rownames(interest) <- rownames(microbes)
  colnames(interest) <- rownames(phenotypes)

  for(i in 1:ncol(interest)){
    ##print(i)
    lasso.out <- group.group.indlasso.lasso(phenotypes[i,],microbes,index,index.subgroup,p.group,tau,
                                            file.group=paste(file,rownames(phenotypes)[i],sep=""),
                                            plots.group,plots.subgroup,
                                            plots.ind,delta.group,delta.subgroup,delta.ind,
                                            use.Gram=use.Gram,standardize=standardize)

    interest[,i] <- lasso.out$sig.variables
  }
  list(interest=interest)
}


################################
## SPARSE GROUP SUBGROUP LASSO #
################################

sparse.group.subgroup <- function(phenotypes,microbes,group.index,subgroup.index,tau,alpha1=0.45,alpha2=0.45,
                                  alpha3=1-alpha1-alpha2,
                                  nlam=100,lambdas=NULL,lambda.accuracy=1e-4,
                                  file="file",plots=TRUE,delta=2,std.y=TRUE){
  y <- t(phenotypes)
  X <- t(microbes)

  N <- length(y)

  if(std.y==TRUE){
    ## Standardize y
    y1 <- make.std(y)
  } else {
    ## Centered response
    y1 <-  make.center(y)
  }


  ## Standardized design matrix X
  X1 <- apply(X,2,make.std)

  ## Put data in a list
  data.list <- list(x=X1,y=y1)


  ## Run sparse group subgroup lasso
  Lasso.out <-  subgroup.SGL(data.list, group.index=group.index,
                             subgroup.index=subgroup.index,
                             type = "linear",nlam=nlam,standardize=FALSE,
                             alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,
                             lambdas=lambdas,
                             lambda.accuracy=lambda.accuracy)
  if(ncol(X1)==1){
    Lasso.out$beta <- t(as.matrix(Lasso.out$beta))
  }

  order.variables <- 0

  ## Samuel's corrections 11/9/2011
  ## use Cp-like criterion to find best descriptive model
  p = dim(X1)[2]
  s = length(Lasso.out$lambdas)
  p.pos = NULL

  RSS = NULL
  for (i in 1:s){
    RSS[i] = sum((y1-predictSGL(Lasso.out, X1, lam=i))**2)
    p.pre = Lasso.out$beta[,i]
    p.pos = c(p.pos,length(p.pre[abs(p.pre)>1e-10]))
  }

  ## Estimated MSE
  MSE <- sd(as.vector(y1)) * sqrt( N / (N-1) )
  MSE <- MSE^2

  p.min = which.min(RSS/MSE+delta*p.pos)
  ##print("subgroup SGL")
  ##print(Lasso.out$lambdas[p.min])

  ## final best descriptive model
  predict.out <- Lasso.out$beta[,p.min]

  ind <- which(abs(predict.out)>1e-10)
  sig.variables <- rep(0,nrow(microbes))
  sig.variables[ind] <- 1

  if(plots==TRUE){
    postscript(paste(file,"_sgl2.eps",sep=""))
    par(mar=c(5, 4, 4, 2)+1)
    plot(1:s,RSS+2*(p.pos),type="l",cex.axis=1.5,cex.lab=1.5,ylab=substitute(M[n](that,p),list(that=delta)), xlab="Steps")

    abline(v=p.min,lty=2)
    dev.off()
  }
  list(sig.variables=sig.variables,predict.out=predict.out)
}

## cross-validation code to get one set of alphas for sparse.group.subgroup
cv.sparse.group.subgroup.alphas <- function(phenotypes,microbes,group.index,subgroup.index,tau,alphas=NULL,
                                     nfold=10,ratio=1,
                                     nlam=100,lambdas=NULL,lambda.accuracy=1e-4,
                                     file="file",plots=FALSE,delta=2,std.y=TRUE){

  ## Matrix to store results
  residmat <- matrix(NA, nrow(alphas), nfold)

  ## do transpose to get correct dimension
  y <- t(phenotypes)
  X <- t(microbes)

  ## Randomly partition the data
  all.folds <- cv.folds(length(y),nfold)

  for(a in 1:nrow(alphas)){
    for(j in seq(nfold)){
      ## data to omit
      omit <- all.folds[[j]]

      phenotypes.fold <- t(y[-omit])
      microbes.fold <- t(X[-omit,,drop=FALSE])
      beta.omit <- sparse.group.subgroup(phenotypes.fold,microbes.fold,group.index,subgroup.index,tau,
                                         alpha1=as.numeric(alphas[a,1]),
                                         alpha2=as.numeric(alphas[a,2]),
                                         alpha3=as.numeric(alphas[a,3]),
                                         nlam=nlam,lambdas=lambdas,lambda.accuracy=lambda.accuracy,
                                         file="file",plots=plots,delta=delta,std.y=std.y)

      ## Find final fit with data omitted
      fit <- X[omit,,drop=FALSE] %*% beta.omit$predict.out

      ## Store residual
      residmat[a,j] <- apply((y[omit]-fit)^2,2,sum)
    }
  }

  cv <- apply(residmat,1,mean)

  ## Check which alpha's lead to min(cv)
  alpha.ind <- which(cv==min(cv,na.rm=TRUE))  ## Ignore NA values
  alpha.opt <- alphas[alpha.ind,]

  ## Is alpha.opt unique? If not, take average.
  if(length(alpha.ind)>1){
    alpha.opt <- apply(alpha.opt,2,mean)
  }

  list(alpha.opt=alpha.opt)
}

make.alphas <- function(alphas.cv=NULL,ratio=1){
  ## Index of fixed alpha
  fixed <- which(lapply(alphas.cv,length)==1)

  ## Index of NON-varying alpha
  other <- which(lapply(alphas.cv,is.null)==TRUE)

  ## Index of varying alpha
  vary <- setdiff(1:length(alphas.cv),c(fixed,other))

  alphas <- matrix(0,nrow=length(alphas.cv[[vary]]),ncol=3)
  colnames(alphas) <- c("alpha1","alpha2","alpha3")

  alphas[,vary] <- alphas.cv[[vary]]
  alphas[,other] <- ratio * alphas.cv[[vary]]
  alphas[,fixed] <- alphas.cv[[fixed]]

  ## Only keep rows such that sum(alphas)<=1
  alphas <- good.alphas(alphas)

  return(alphas)
}

## function to keep those alphas such that sum(alphas) <=1
good.alphas <- function(alphas){
  ## Only keep rows such that sum(alphas)<=1
  ind <- which(apply(alphas,1,sum)==1)
  alphas <- alphas[ind,]

  return(alphas)
}

cv.sparse.group.subgroup <- function(phenotypes,microbes,group.index,subgroup.index,tau,
                                     alphas.cv.range=c(seq(0.01,0.1,by=0.01),seq(0.15,0.95,by=0.05)),
                                     nfold=10,nlam=100,lambdas=NULL,lambda.accuracy=1e-4,
                                     file="file",plots=FALSE,delta=2,std.y=TRUE,sam.implement=FALSE){

  if(sam.implement==FALSE){
    alphas.cv <- list(alpha1.cv = alphas.cv.range, alpha2.cv = alphas.cv.range, alpha3.cv = alphas.cv.range)
    alphas <- do.call(expand.grid,alphas.cv)
    alphas <- good.alphas(alphas)

    cv.out <- cv.sparse.group.subgroup.alphas(phenotypes,microbes,group.index,subgroup.index,
                                              tau,alphas=alphas,
                                              nfold=nfold,ratio=ratio,
                                              nlam=nlam,lambdas=lambdas,lambda.accuracy=lambda.accuracy,
                                              plots=FALSE,delta=delta,std.y=std.y)
    alpha.opt <- as.numeric(cv.out$alpha.opt)

  } else {


###########################################################
    ## First application: alpha3=1/3, alpha1=ratio * alpha2, ratio=1 ##
###########################################################

    alphas.cv <- list(alpha1.cv = NULL, alpha2.cv = alphas.cv.range, alpha3.cv = 1/3)
    ratio <- 1
    alphas <- make.alphas(alphas.cv,ratio)

    cv.out <- cv.sparse.group.subgroup.alphas(phenotypes,microbes,group.index,subgroup.index,
                                              tau,alphas=alphas,
                                              nfold=nfold,ratio=ratio,
                                              nlam=nlam,lambdas=lambdas,lambda.accuracy=lambda.accuracy,
                                              plots=FALSE,delta=delta,std.y=std.y)
    alpha.opt <- cv.out$alpha.opt

#############################################################################################
    ## Second application: alpha2= alphas[2], alpha1=ratio * alpha3, ratio=alphas[1]/alphas[3] ##
#############################################################################################

    alphas.cv <- list(alpha1.cv = NULL, alpha2.cv = alpha.opt[2], alpha3.cv = alphas.cv.range)
    ratio <- alpha.opt[1] / alpha.opt[3]
    alphas <- make.alphas(alphas.cv,ratio)

    cv.out <- cv.sparse.group.subgroup.alphas(phenotypes,microbes,group.index,subgroup.index,
                                              tau,alphas=alphas,
                                              nfold=nfold,ratio=ratio,
                                              nlam=nlam,lambdas=lambdas,lambda.accuracy=lambda.accuracy,
                                              plots=FALSE,delta=delta,std.y=std.y)
    alpha.opt <- cv.out$alpha.opt


############################################################################################
    ## Third application: alpha1= alphas[1], alpha2=ratio * alpha3, ratio=alphas[2]/alphas[3] ##
############################################################################################

    alphas.cv <- list(alpha1.cv = alpha.opt[1], alpha2.cv = NULL, alpha3.cv = alphas.cv.range)
    ratio <- alpha.opt[2] / alpha.opt[3]
    alphas <- make.alphas(alphas.cv,ratio)

    cv.out <- cv.sparse.group.subgroup.alphas(phenotypes,microbes,group.index,subgroup.index,
                                              tau,alphas=alphas,
                                              nfold=nfold,ratio=ratio,
                                              nlam=nlam,lambdas=lambdas,lambda.accuracy=lambda.accuracy,
                                              plots=FALSE,delta=delta,std.y=std.y)
    alpha.opt <- cv.out$alpha.opt
  }
############################################################
  ## Re-run sparse group-subgroup with optimal alpha values ##
############################################################
  Lasso.out <-  sparse.group.subgroup(phenotypes,microbes,group.index,subgroup.index,tau,
                                      alpha1=as.numeric(alpha.opt[1]),
                                      alpha2=as.numeric(alpha.opt[2]),
                                      alpha3=as.numeric(alpha.opt[3]),
                                      nlam=nlam,lambdas=lambdas,lambda.accuracy=lambda.accuracy,
                                      plots=FALSE,delta=delta,std.y=std.y)

  list(sig.variables = Lasso.out$sig.variables, alphas=alpha.opt)
}

sparse.group.subgroup.computations <- function(microbes,phenotypes,group.index,subgroup.index,tau,alpha1,alpha2,
                                               alpha3,nlam,lambdas,
                                               lambda.accuracy,
                                               file.group="file",plots.group=TRUE,delta.group=2,format.data=TRUE,
                                               cv.criterion=FALSE,nfold=10,alphas.cv.range=seq(0.1,0.95,by=0.05)){
  if(format.data==TRUE){
    phenotypes <- phenotypes[-c(1,2),]
  } else {
    phenotypes <- phenotypes
  }
  interest <- matrix(0,nrow=nrow(microbes),ncol=nrow(phenotypes))
  interest <- as.data.frame(interest)
  rownames(interest) <- rownames(microbes)
  colnames(interest) <- rownames(phenotypes)

  ## Store alpha values
  alpha.out <- matrix(0,nrow=3,ncol=nrow(phenotypes))

  for(i in 1:ncol(interest)){
    ##print(i)
    if(cv.criterion==TRUE){

      lasso.out <- cv.sparse.group.subgroup(phenotypes[i,],microbes,group.index,subgroup.index,tau,
                                            alphas.cv.range=alphas.cv.range,
                                            nfold=nfold,nlam=nlam,lambdas=lambdas,
                                            lambda.accuracy=lambda.accuracy,
                                            file=file.group,plots=plots.group,delta=delta.group)
      alpha.out[,i] <- lasso.out$alphas

    } else {

      lasso.out <- sparse.group.subgroup(phenotypes[i,],microbes,group.index,subgroup.index,
                                         tau,alpha1,alpha2,alpha3,nlam=nlam,lambdas=lambdas,
                                         lambda.accuracy=lambda.accuracy,
                                         file=file.group,plots=plots.group,delta=delta.group)
    }
    interest[,i] <- lasso.out$sig.variables
  }
  list(interest=interest,alpha.out=alpha.out)
}




#######################
## SPARSE GROUP LASSO #
#######################

sparse.group <- function(phenotypes,microbes,index,tau,alpha=0.95,file="file",plots=TRUE,delta=2,std.y=TRUE){
  y <- t(phenotypes)
  X <- t(microbes)

  N <- length(y)

  if(std.y==TRUE){
    ## Standardize y
    y1 <- make.std(y)
  } else {
    ## Centered response
    y1 <-  make.center(y)
  }


  ## Standardized design matrix X
  X1 <- apply(X,2,make.std)

  ## Put data in a list
  data.list <- list(x=X1,y=y1)

  ## Run sparse group lasso
  Lasso.out <-  mySGL(data.list, index=index, type = "linear",nlam=100,
                      standardize=FALSE,alpha=alpha,lambdas=NULL)

  if(ncol(X1)==1){
    Lasso.out$beta <- t(as.matrix(Lasso.out$beta))
  }

  order.variables <- 0

  ## Samuel's corrections 11/9/2011
  ## use Cp-like criterion to find best descriptive model
  p = dim(X1)[2]
  s = length(Lasso.out$lambdas)
  p.pos = NULL

  RSS = NULL
  for (i in 1:s){
    ##print(i)
    RSS[i] = sum((y1-predictSGL(Lasso.out, X1, lam=i))**2)
    p.pre = Lasso.out$beta[,i]
    p.pos = c(p.pos,length(p.pre[abs(p.pre)>1e-10]))
  }

  ## Estimated MSE
  MSE <- sd(as.vector(y1)) * sqrt( N / (N-1) )
  MSE <- MSE^2

  p.min = which.min(RSS/MSE+delta*p.pos)


  ## final best descriptive model
  predict.out <- Lasso.out$beta[,p.min]


  ind <- which(abs(predict.out)>1e-10)
  sig.variables <- rep(0,nrow(microbes))
  sig.variables[ind] <- 1

  if(plots==TRUE){
    ##postscript(paste(file,"_sgl1.eps",sep=""))
    ##plot(Lasso.out,cex.axis=1.5,cex.lab=1.5)
    ##dev.off()

    postscript(paste(file,"_sgl2.eps",sep=""))
    par(mar=c(5, 4, 4, 2)+1)
    plot(1:s,RSS+2*(p.pos),type="l",cex.axis=1.5,cex.lab=1.5,ylab=substitute(M[n](that,p),list(that=delta)), xlab="Steps")

    abline(v=p.min,lty=2)
    dev.off()
  }
  list(sig.variables=sig.variables,predict.out=predict.out)
}

## Code for cross-validation to select \alpha
## nfold : number of folds in cross-validation

cv.sparse.group <- function(phenotypes,microbes,index,tau,delta,nfold=10,alpha.cv=seq(0.05,0.95,by=0.1),
                            plots=FALSE, std.y=TRUE){
  ## do transpose to get correct dimension
  y <- t(phenotypes)
  X <- t(microbes)

  ## Randomly partition the data
  all.folds <- cv.folds(length(y),nfold)

  ## Matrix to store results
  residmat <- matrix(0, length(alpha.cv), nfold)

  for(a in 1:length(alpha.cv)){
    ##print(a)
    alpha <- alpha.cv[a]
    for(j in seq(nfold)){
      ##print(j)
      ## data to omit
      omit <- all.folds[[j]]

      ## put data into a list (remove omitted data)
      phenotypes.fold <- t(y[-omit])
      microbes.fold <- t(X[-omit,,drop=FALSE])
      beta.omit <- sparse.group(phenotypes.fold,microbes.fold,index,tau,alpha=alpha,plots=plots,delta=delta,
                                std.y=std.y)

      ## Find final fit with data omitted
      fit <- X[omit,,drop=FALSE] %*% beta.omit$predict.out

      ## Store residual
      residmat[a,j] <- apply((y[omit]-fit)^2,2,sum)

    }
  }

  cv <- apply(residmat,1,mean)

  ## Check which alpha's lead to min(cv)
  alpha.ind <- which(cv==min(cv))
  alpha.opt <- alpha.cv[alpha.ind]

  ## Is alpha.opt unique? If not, take average.
  if(length(alpha.opt)>1){
    alpha.opt <- mean(alpha.opt)
  }

  ## Re-run mySGL with optimal alpha value
  Lasso.out <- sparse.group(phenotypes,microbes,index,tau,alpha=alpha.opt,plots=plots,delta=delta,std.y=std.y)
  list(sig.variables = Lasso.out$sig.variables,alpha=alpha.opt)
}



########################################################################
## GROUP LASSO at GROUP LEVEL and SPARSE GROUP LASSO at subgroup level #
########################################################################
## delta.group : delta applied to C_p criterion for group lasso
## delta.subgroup   : delta applied to C_p critierian for group lasso among subgroups

group.sgl.lasso <- function(phenotypes,microbes,index,index.subgroup,p.group,tau,alpha,
                            file.group="file",plots.group=TRUE,plots.subgroup=TRUE,delta.group=2,
                            delta.subgroup=2,standardize=TRUE,
                            cv.criterion=FALSE, nfold=10, alpha.cv=seq(0.05,0.95,by=0.1)){

  ## group lasso to groups
  group.lasso.out <- group.lasso(phenotypes,microbes,index,tau,file=file.group,
                                 plots=plots.group,
                                 delta=delta.group,standardize=standardize)

  sig.variables <- group.lasso.out$sig.variables

  main.ind <- which(sig.variables==1)

  ## setting up to apply group Lasso to subgroups
  tmp.subgroup <- as.vector(t(index.subgroup))
  tmp.subgroup <- tmp.subgroup[!is.na(tmp.subgroup)]
  new.index.subgroup <- tmp.subgroup[main.ind]

  ## Initial set of alpha value. Will change if cv.criterion==TRUE
  alpha <- alpha

  if(sum(main.ind)!=0){
    X.tmp <- microbes[main.ind,]
    if(cv.criterion==TRUE){

      subgroup.lasso.out <- cv.sparse.group(phenotypes,X.tmp,new.index.subgroup,
                                            tau,delta=delta.subgroup,nfold=nfold,alpha.cv=alpha.cv,
                                            plots=plots.subgroup)
      alpha <- subgroup.lasso.out$alpha

    } else {
      subgroup.lasso.out <- sparse.group(phenotypes,X.tmp,new.index.subgroup,
                                         tau,alpha,
                                         file=paste(file.group,"_sgl_no_",k,
                                           sep=""),plots=plots.subgroup,
                                         delta=delta.subgroup)
    }
    sig.variables[main.ind] <- subgroup.lasso.out$sig.variables
  }

  list(sig.variables=sig.variables,alpha=alpha)

}

group.sgl.lasso.computations <- function(microbes,phenotypes,index,index.subgroup,p.group,tau,alpha,
                                         file.group="file",plots.group=TRUE,plots.subgroup=TRUE,delta.group=2,
                                         delta.subgroup=2,format.data=TRUE,standardize=TRUE,
                                         cv.criterion=FALSE,nfold=10, alpha.cv=seq(0.05,0.95,by=0.1)){

  if(format.data==TRUE){
    phenotypes <- phenotypes[-c(1,2),]
  } else {
    phenotypes <- phenotypes
  }
  interest <- matrix(0,nrow=nrow(microbes),ncol=nrow(phenotypes))
  interest <- as.data.frame(interest)
  rownames(interest) <- rownames(microbes)
  colnames(interest) <- rownames(phenotypes)

  ## store alpha values
  alpha.out <- matrix(0,nrow=1,ncol=nrow(phenotypes))

  for(i in 1:ncol(interest)){
    ##print(i)

    lasso.out <- group.sgl.lasso(phenotypes[i,],microbes,index,index.subgroup,p.group,tau,alpha,
                                 file=file.group,plots.group=plots.group,plots.subgroup=plots.subgroup,
                                 delta.group=delta.group,delta.subgroup=delta.subgroup,standardize=standardize,
                                 cv.criterion=cv.criterion,nfold=nfold,alpha.cv=alpha.cv)

    interest[,i] <- lasso.out$sig.variables
    alpha.out[,i] <- lasso.out$alpha
  }
  list(interest=interest,alpha.out=alpha.out)
}





###################################################
## SPARSE GROUP LASSO for group Lasso simulations #
###################################################

sparse.group.computations <- function(microbes,phenotypes,index,p.group,tau,alpha,
                                      file.group="file",plots.group=TRUE,delta.group=2,format.data=TRUE,
                                      cv.criterion=FALSE, nfold=10, alpha.cv=seq(0.05,0.95,by=0.1)){
  if(format.data==TRUE){
    phenotypes <- phenotypes[-c(1,2),]
  } else {
    phenotypes <- phenotypes
  }
  interest <- matrix(0,nrow=nrow(microbes),ncol=nrow(phenotypes))
  interest <- as.data.frame(interest)
  rownames(interest) <- rownames(microbes)
  colnames(interest) <- rownames(phenotypes)

  ## Store alpha values
  alpha.out <- matrix(0,nrow=1,ncol=nrow(phenotypes))

  for(i in 1:ncol(interest)){
    ##print(i)
    if(cv.criterion==TRUE){
      lasso.out <- cv.sparse.group(phenotypes[i,],microbes,index,tau,delta=delta.group,nfold=nfold,alpha.cv=alpha.cv,
                      plots=plots.group)
      alpha.out[,i] <- lasso.out$alpha

    } else {
      lasso.out <- sparse.group(phenotypes[i,],microbes,index,tau,alpha,file=file.group,
                                plots=plots.group,delta=delta.group)
    }
    interest[,i] <- lasso.out$sig.variables
  }
  list(interest=interest, alpha.out=alpha.out)
}







#####################################
## Data generation for Group lasso ##
#####################################
## N 	    : sample size
## p.group : number of covariates in each group
## p.subgroup : matrix of number of covariates in each subgroup (each row corresponds to a group)
## beta.coef : matrix of coefficient vectors (rows: groups, columns: coefficient values)
## sigma   : model error covariance
## index   : vector indicating which covariates are in which groups

data.group <- function(N,p.subgroup,beta.coef,sigma){

    ########################################
    ## Form beta vector ##
    ########################################

    beta <- as.vector(t(beta.coef))	## Make into a vector
    beta <- beta[!is.na(beta)]		## Remove NA values

    ##################################################################
    ## Make p.subgroup into a vector ##
    ##################################################################

    p.subgroup.vector <- p.subgroup

    ######################################
    ## Form covariates ##
    ######################################
    ## X \sim Normal(mu,Sigma)

    no.subgroup <- length(p.subgroup.vector)
    Sigma.list <- vector("list",no.subgroup)
    for(j in 1:no.subgroup){
        Sigma.list[[j]] <- 0.7 * matrix(1,nrow=p.subgroup.vector[j],
                                        ncol=p.subgroup.vector[j]) +
            0.3 * diag(1,nrow=p.subgroup.vector[j])
    }

    ## Form Sigma covariance matrix
    X.Sigma <- as.matrix(bdiag(Sigma.list))

    ## Form mean of normal distribution
    X.mu <- rep(0,length(beta))

    ## Form covariates
    X <- mvrnorm(N,X.mu,X.Sigma)
    colnames(X) <- paste("X_",seq(1,sum(p.group)),sep="")


    ################################################
    ## Form response vector ##
    ################################################

    y <- X %*% beta + rnorm(N,mean=0,sd=sigma)


    list(y=y,X=X)
}
