source("/Users/junkyungkim/Desktop/ResearchWork/R_documentation/sparse_gplasso/sgsl/R/main_code.R")

###############
# Source code #
###############

dyn.load("/Users/junkyungkim/Desktop/ResearchWork/R_documentation/sparse_gplasso/Cpp_code/mylin.dll")
dyn.load("/Users/junkyungkim/Desktop/ResearchWork/R_documentation/sparse_gplasso/Code-on-Website/Cpp_code/mylin.dll")
dyn.load("../Cpp_code/mylin.dll")
source("../main/subgroup_SGL.R")
source("../main/main_code.R")

## Form beta coefficients for each group
delta=2; nsimu=1; ncv=10; percents.range=c(50,60,70,80,90,100)
tau=0.94;alpha.SGL=0.95;alpha=0.95;alpha1=0.45;alpha2=0.45;
alpha3=1-alpha1-alpha2;lambdas=NULL;nlam=100
delta.sub =2; delta.ind=2
two.alphas=TRUE
use.Gram <-  TRUE
group.standardize=TRUE
nfold=10
alphas.cv.range=seq(0.1,0.95,by=0.05)
run.alphas.range=FALSE
lambda.accuracy=1e-4

N=30;
L=10;
p.in.group =8;
p=L * p.in.group;
sigma <- sqrt(1);
beta.coef <- matrix(0,nrow=2*L,ncol=(p/L)/2)
beta.coef[1,] <- c(0.1,0.2,0.3,0.4)
beta.coef[2,] <- c(0.1,0.2,0.3,0.4)
## index for sub-groups, subgroup labels are unique!
p.group <- rep(p/L,L)
index.subgroup <- matrix(NA,nrow=L,ncol=p)
tmp <- 0
for(k in 1:L){
    if(k==1){
        index.subgroup[k,1:p.group[k]] <- c(rep(1,(p/L)/2),rep(2,(p/L)/2))
    } else {
        ind <- 1:p.group[k] + sum(p.group[(k-1):1])
        index.subgroup[k,ind] <- c(rep(k+tmp,(p/L)/2),rep(k+tmp+1,(p/L)/2))
    }
    tmp <- tmp + 1
}
tools <- form.tools(index.subgroup)
index <- tools$index
p.group <- tools$p.group
p.subgroup <- tools$p.subgroup
p <- sum(p.group) 		## Total number of covariates







out.rownames <- paste("X_",seq(1,p),sep="")  ## Names of microbes


## results from pure Lasso
out.lasso <- as.data.frame(matrix(0,nrow = p, ncol=length(delta)))
colnames(out.lasso) <- paste("lasso.delta.",delta,sep="")
rownames(out.lasso) <- paste("X_",seq(1,p),sep="")


## results from pure group Lasso
out.pure.gp.lasso <-  as.data.frame(matrix(0,nrow = p, ncol=length(delta)))
colnames(out.pure.gp.lasso) <- paste("pure.gp.lasso.delta.",delta,sep="")
rownames(out.pure.gp.lasso) <- paste("X_",seq(1,p),sep="")


## results from SGL
out.sgl <-  as.data.frame(matrix(0,nrow = p, ncol=length(delta)))
colnames(out.sgl) <- paste("sgl.delta.",delta,sep="")
rownames(out.sgl) <- paste("X_",seq(1,p),sep="")

## results from SGL.cv (we fix delta=2)
mult.cv.sgl <- as.data.frame(matrix(0,nrow=p,ncol=nsimu,
                                    dimnames = list(paste("X_",seq(1,p),sep=""),
                                                    paste("mult.cv.sgl.all.",seq(1,nsimu),sep=""))))

mult.cv.alpha.sgl <- as.data.frame(matrix(0, nrow=1, ncol=ncv,
                                          dimnames=list("alpha",seq(1,ncv))))

mult.cv.sgl.summary <- as.data.frame(matrix(0,nrow=p,ncol=length(percents.range),
                                            dimnames=list(paste("X_",seq(1,p),sep=""),
                                                          paste("mult.cv.sgl.summary.",percents.range,sep=""))))


## results from SGL over alpha.SGL range
out.sgl.alpharange <-  as.data.frame(matrix(0,nrow = p, ncol=length(delta) * length(alpha.SGL)))

tmp <- paste(".alpha.",alpha.SGL,sep="")
tmp.names <- NULL

for(k in 1:length(delta)){
    tmp.names <- c(tmp.names,paste("sgl.delta.",delta[k],tmp,
                                   sep=""))
}

colnames(out.sgl.alpharange) <- tmp.names
rownames(out.sgl.alpharange) <- paste("X_",seq(1,p),sep="")




## results for group lasso at group level and group lasso at subgroup level
out.gp.gp.lasso <- as.data.frame(matrix(0,nrow = p,ncol=length(delta) * length(delta.sub)))

## results for group lasso at group level and SGL at subgroup level
out.group.sgl.lasso <- as.data.frame(matrix(0,nrow = p, ncol=length(delta) * length(delta.sub)))

## results from group lasso at group level and SGL at subgroup level with CV (we fix delta=2)
mult.cv.group.sgl.lasso <- as.data.frame(matrix(0,nrow=p,ncol=nsimu,
                                                dimnames = list(paste("X_",seq(1,p),sep=""),
                                                                paste("mult.cv.group.sgl.lasso.all.",seq(1,nsimu),sep=""))))

mult.cv.alpha.group.sgl.lasso <- as.data.frame(matrix(0, nrow=1, ncol=ncv,
                                                      dimnames=list("alpha",seq(1,ncv))))

mult.cv.group.sgl.lasso.summary <- as.data.frame(matrix(0,nrow=p,ncol=length(percents.range),
                                                        dimnames=list(paste("X_",seq(1,p),sep=""),
                                                                      paste("mult.cv.group.sgl.lasso.summary.",
                                                                            percents.range,sep=""))))
## column names for group-group lasso
tmp <- paste(".sub.delta.",delta.sub,sep="")

tmp.names <- NULL
tmp.sgl.names <- NULL

for(k in 1:length(delta)){
    tmp.names <- c(tmp.names,paste("gp.lasso.delta.",delta[k],tmp,sep=""))
    tmp.sgl.names <- c(tmp.sgl.names,paste("gp.sgl.lasso.delta.",delta[k],tmp,
                                           sep=""))
}
colnames(out.gp.gp.lasso) <- 	tmp.names
rownames(out.gp.gp.lasso) <- paste("X_",seq(1,p),sep="")

colnames(out.group.sgl.lasso) <- tmp.sgl.names
rownames(out.group.sgl.lasso) <- paste("X_",seq(1,p),sep="")


## results for group lasso at group level, group lasso at subgroup level, and lasso among all subgroups remaining

out.gp.gp.indlasso.lasso <- as.data.frame(matrix(0,nrow = p, ncol=length(delta) * length(delta.sub) *
                                                     length(delta.ind) ))

## column names for group-group lasso
tmp <- paste(".ind.delta.",delta.ind,sep="")
tmp.names.sub <- NULL

for(k in 1:length(delta.sub)){
    tmp.names.sub <- c(tmp.names.sub,paste(".sub.delta.",delta.sub[k],
                                           tmp,sep=""))
}

tmp.names <- NULL
for(k in 1:length(delta)){
    tmp.names <- c(tmp.names,paste("gp.lasso.delta.",delta[k],tmp.names.sub,
                                   sep=""))
}
colnames(out.gp.gp.indlasso.lasso) <- tmp.names
rownames(out.gp.gp.indlasso.lasso) <- paste("X_",seq(1,p),sep="")

## results from group lasso at group level and SGL at subgroup level over a range of alpha values
out.group.sgl.lasso.range <- as.data.frame(matrix(0,nrow = p, ncol=length(delta) * length(delta.sub) *
                                                      length(alpha) ))

## column names for group-group lasso
tmp <- paste(".alpha.",alpha,sep="")
tmp.names.sub <- NULL

for(k in 1:length(delta.sub)){
    tmp.names.sub <- c(tmp.names.sub,paste(".sgl.delta.",delta.sub[k],tmp,
                                           sep=""))
}

tmp.names <- NULL
for(k in 1:length(delta)){
    tmp.names <- c(tmp.names,paste("gp.lasso.delta.",delta[k],tmp.names.sub,
                                   sep=""))
}
colnames(out.group.sgl.lasso.range) <-         tmp.names
rownames(out.group.sgl.lasso.range) <- paste("X_",seq(1,p),sep="")




## results from sparse group subgroup Lasso
out.sparse.gp.subgp.lasso <- as.data.frame(matrix(0,nrow = p, ncol=length(delta)))
colnames(out.sparse.gp.subgp.lasso) <- paste("sparse.gp.subgp.lasso.delta.",delta,sep="")
rownames(out.sparse.gp.subgp.lasso) <- paste("X_",seq(1,p),sep="")

## results from sparse group subgroup Lasso with CV (we fix delta=2)
mult.cv.sparse.gp.subgp.lasso <- as.data.frame(matrix(0,nrow=p,ncol=nsimu,
                                                      dimnames = list(paste("X_",seq(1,p),sep=""),
                                                                      paste("mult.cv.sparse.gp.subgp.lasso.all.",
                                                                            seq(1,nsimu),sep=""))))

mult.cv.alpha.sparse.gp.subgp.lasso <- as.data.frame(matrix(0, nrow=3, ncol=ncv,
                                                            dimnames=list(paste("alpha",seq(1,3),sep="")
                                                                          ,seq(1,ncv))))

mult.cv.sparse.gp.subgp.lasso.summary <- as.data.frame(matrix(0,nrow=p,ncol=length(percents.range),
                                                              dimnames=list(paste("X_",seq(1,p),sep=""),
                                                                            paste("mult.cv.sparse.gp.subgp.lasso.summary.",
                                                                                  percents.range,sep=""))))


## results from sparse group subgroup Lasso at different alpha1, alpha2 combinations
##tmp <- paste(".delta.",delta,sep="")
tmp.names <- NULL

for(k in 1:length(alpha1)){
    for(r in 1:length(alpha2)){
        if(two.alphas==TRUE){
            if(alpha1[k]+alpha2[r]<1){
                tmp.names <- c(tmp.names,paste(".a1.",alpha1[k],".a2.",alpha2[r],sep=""))
            }
        } else {
            for(s in 1:length(alpha3)){
                tmp.names <- c(tmp.names,paste(".a1.",alpha1[k],".a2.",alpha2[r],".a3.",alpha3[s],sep=""))
            }
        }
    }
}

tmp.names2 <- NULL
for(k in 1:length(delta)){
    tmp.names2 <- c(tmp.names2,paste("sub.SGL.del.",delta[k],tmp.names,sep=""))
}

out.subgroup.SGL.all.alphas  <- as.data.frame(matrix(0,nrow = p, ncol=length(tmp.names2)))

colnames(out.subgroup.SGL.all.alphas) <-         tmp.names2
rownames(out.subgroup.SGL.all.alphas) <- paste("X_",seq(1,p),sep="")





data <- data.group(N,p.group,beta.coef,sigma)
phenotypes <- as.data.frame(t(data$y))
microbes   <- as.data.frame(t(data$X))
rownames(phenotypes) <- "phenotype"

tmp.alphas <- 1
tmp.sgl.alphas <- 1


out.lasso <- out.lasso + lasso.computations(microbes,phenotypes,plots=FALSE,file="lasso_",
                                            format.data=FALSE,
                                            delta=delta,use.Gram=use.Gram)$interest
out.pure.gp.lasso <- out.pure.gp.lasso +
    pure.grplasso.computations(microbes,phenotypes,
                               index,p.group,tau,file.group="pure_group_lasso_",
                               plots.group="FALSE",delta.group=delta,
                               format.data=FALSE,standardize=group.standardize)$interest
if(length(alpha.SGL)>1){
    for(qq in 1:length(alpha.SGL)){
        out.sgl.alpharange[,qq] <- out.sgl.alpharange[,qq] +
            sparse.group.computations(microbes,phenotypes,index,p.group,tau,
                                      alpha=alpha.SGL[qq],
                                      file.group="SGL_",plots.group=FALSE,
                                      delta.group=delta,
                                      format.data=FALSE)$interest
        ##print(alpha.SGL[qq])
    }
} else {
    out.sgl <- out.sgl +
        sparse.group.computations(microbes,phenotypes,index,p.group,tau,alpha.SGL,
                                  file.group="SGL_",plots.group=FALSE,
                                  delta.group=delta,format.data=FALSE)$interest
}
for(v in 1:ncv){
    mult.cv.sgl.tmp <- sparse.group.computations(microbes,phenotypes,index,p.group,tau,alpha.SGL,
                                                 file.group="SGL_",plots.group=FALSE,
                                                 delta.group=delta,format.data=FALSE,
                                                 cv.criterion=TRUE,nfold=nfold,
                                                 alpha.cv=setdiff(alphas.cv.range,0.5))

    mult.cv.sgl[,j] <- mult.cv.sgl[,j] + mult.cv.sgl.tmp$interest
    mult.cv.alpha.sgl[,v] <- mult.cv.alpha.sgl[,v] + mult.cv.sgl.tmp$alpha.out
}
if(run.alphas.range==TRUE){
    if(two.alphas==TRUE){
        for(kk in 1:length(alpha1)){
            for(rr in 1:length(alpha2)){
                if(alpha1[kk] + alpha2[rr] < 1){
                    out.subgroup.SGL.all.alphas[,tmp.alphas] <- out.subgroup.SGL.all.alphas[,tmp.alphas] +
                        sparse.group.subgroup.computations(microbes,phenotypes,
                                                           index,index.subgroup,tau,
                                                           alpha1=alpha1[kk],alpha2=alpha2[rr],
                                                           alpha3=1-alpha1[kk]-alpha2[rr],nlam=nlam,
                                                           lambdas=lambdas,lambda.accuracy=lambda.accuracy,
                                                           file.group="sparse_group_subgroup_lasso_",plots.group=FALSE,
                                                           delta.group=delta,format.data=FALSE)$interest
                    tmp.alphas <- tmp.alphas + 1
                    ##print(c("alpha1",alpha1[kk],"alpha2",alpha2[rr]))
                    ##print(c("tmp.alphas",tmp.alphas))
                }
            }
        }
    } else {
        for(kk in 1:length(alpha1)){
            for(rr in 1:length(alpha2)){
                for(ss in 1:length(alpha3)){
                    out.subgroup.SGL.all.alphas[,tmp.alphas] <- out.subgroup.SGL.all.alphas[,tmp.alphas] +
                        sparse.group.subgroup.computations(microbes,phenotypes,
                                                           index,index.subgroup,tau,
                                                           alpha1=alpha1[kk],alpha2=alpha2[rr],alpha3=alpha3[ss],
                                                           nlam=nlam,
                                                           lambdas=lambdas,lambda.accuracy=lambda.accuracy,
                                                           file.group="sparse_group_subgroup_lasso_",plots.group=FALSE,
                                                           delta.group=delta,format.data=FALSE)$interest
                    tmp.alphas <- tmp.alphas + 1
                    ##print(c("alpha1",alpha1[kk],"alpha2",alpha2[rr]))
                    ##print(c("tmp.alphas",tmp.alphas))
                }
            }
        }
    }
} else {

    out.sparse.gp.subgp.lasso <- out.sparse.gp.subgp.lasso +
        sparse.group.subgroup.computations(microbes,phenotypes,
                                           index,index.subgroup,tau,alpha1,alpha2,alpha3,nlam,lambdas,
                                           lambda.accuracy,
                                           file.group="sparse_group_subgroup_lasso_",plots.group=FALSE,
                                           delta.group=delta,format.data=FALSE)$interest
}
