## required packages:
## foreach, doMC, pspline, lokern, reshape2, monomvn

suppressMessages(library(foreach))
suppressMessages(library(doMC))

#' Function to split a vector into a list at specified positions 
#' adapted from http://stackoverflow.com/questions/16357962/r-split-numeric-vector-at-position
#' @param x vetor of data
#' @param pos positions to split the data vector
#' @return a list of split data
splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))

#' Function to smooth data and calculate derivative with glkerns
#' @param x vector of data
#' @param t time stamp of each point
#' @param deriv order of derivative to calculate
#' @param t.out time stamp of output (default: to be the same as input)
#' @return a vector of fitted value after smoothing
smooth.deri.glkerns <- function(x, t, deriv=0, t.out=t, smethod='glkerns', ...){
    suppressMessages(library(lokern))
    fit.tmp <- glkerns(t, x, deriv=deriv, x.out = t.out, bandwidth=10)$est
    res <- abs(x - fit.tmp)
    res.med <- median(res)
    res.mad <- mad(res)
    outliers <- (res - res.med) / (res.mad+1e-16) > 5
    outliers[c(1:3,length(outliers))] <- FALSE
    if(smethod=='glkerns'){
        glkerns(t[!outliers], x[!outliers], deriv=deriv, x.out = t.out,...)$est
    }else{
        lokerns(t[!outliers], x[!outliers], deriv=deriv, x.out = t.out,...)$est
    }
}

#' Function to smooth data and calculate derivative with pspline
#' @param x vector of data
#' @param t time stamp of each point
#' @param deriv order of derivative to calculate
#' @param t.out time stamp of output (default: to be the same as input)
#' @param method spline fitting method used for "pspline" (default: 3)
#' @param norder oder of spline basis (default: 3)
#' @return a vector of fitted value after smoothing
smooth.deri.pspline <- function(x, t, deriv=0, t.out=t, smethod=3, norder=3, ...){
    suppressMessages(library(pspline))
    ## initial fit to find outliers
    fit.tmp <- c(smooth.Pspline(t, x,method=smethod, norder=norder, spar=1e10)$ysmth)
    res <- abs(x - fit.tmp)
    res.med <- median(res)
    res.mad <- mad(res)
    outliers <- (res - res.med) / (res.mad+1e-16) > 5
    outliers[c(1:3,length(outliers))] <- FALSE
    if (2*norder + 1 >= length(t[!outliers])){outliers[1:length(outliers)] <- FALSE}
    c(predict(smooth.Pspline(t[!outliers], x[!outliers],method=smethod, norder=norder, ...), t.out, deriv ) )
}

#' Function to smooth data and calculate derivative with pspline
#' @param m a matrix of data (variable in columns and time in rows)
#' @param t matrix or vector of time points for each data
#' @param breaks break points to re-smooth the data
#' @param t.out time stamps of output (default: to be the same as input)
#' @param deriv order of derivative to calculate
#' @param ncpu number of cpus (default: 10)
#' @param method smoothing method used to calculate gradients (default: pspline)
#' @return a matrix of gradients calculated (deriv=1), a matrix of smoothed data (deriv=0)
SmoothGradient <- function(m, t, breaks, t.out = t, deriv=1, ncpu=10, method='pspline', ...){
    m[is.na(m)] <- 0
    registerDoMC(ncpu) 
    p <- ncol(m)
    if(is.null(ncol(t))){
        t.mat <- matrix(rep(t, p), ncol=p)
    }
    idList <- splitAt(seq_along(t) , breaks)
    ## assertion on method used
    try(if(!method %in% c('pspline', 'glkerns'))  stop('Invalid method!'))
    
    smoothed <- foreach(c = 1:p, .combine='cbind') %dopar% {
        ts <- t.mat[,c]
        x <- m[,c]        
        res <- foreach(id=idList, .combine = 'c') %do%{
            if(method == 'pspline'){
                smooth.deri.pspline(x[id],ts[id],deriv=deriv, t.out=ts[id],...)
            }else if(method == 'glkerns'){
                smooth.deri.glkerns(x[id],ts[id],deriv=deriv, t.out=ts[id],...)
            }
        }        
        res
    }    
    smoothed    
}

#' Function to transform a matrix of compositions with addative log ratio (alr)
#' @param m matrix of data (variables in columns, samples in rows) 
#' @param refSp reference variable
#' @return transformed matrix
alr.transformation <- function(m, refSp){
    if(min(m, na.rm = TRUE)==0) {
        pseudo <- 10^(floor(log10(min(m[m!=0], na.rm = TRUE))))
        m[m==0] <- pseudo
    }
    log(m/m[,refSp])
}

#' Function to inverse transform a matrix from addative log ratio (alr) to compostions
#' @param alr matrix of data (variables in columns, samples in rows) 
#' @param refSp reference variable
#' @return inverse transformed matrix
alr.transformation.inv <- function(alr,refSp){
    tmp <- exp(alr)
    tmp[, refSp] <- 1
    tmp/rowSums(tmp,na.rm=TRUE)
}

#' Function to transform a matrix of compositions with centered log ratio (clr)
#' @param m matrix of data (variables in columns, samples in rows)
#' @return transformed matrix
clr.transformation <- function(m){
    ##m: matrix of data (variables in columns, samples in rows) 
    log(m) - rowMeans(log(m))
}

#' Function to detect outliers in time serires
#' @param x vector of time series
#' @param t vector of time stamps
#' @param dev deviation (dev * mad) from the median to be considered as outliers
tsoutliers <- function(x,t,dev=2){
    NAs <- is.na(x)
    score <- rep(Inf, length(x))
    x <- x[!NAs]
    t <- t[!NAs]
    resid <- abs(residuals(loess(x ~ t)))
    score[!NAs] <- (resid - median(resid))>mad(resid) * dev
    
    ## resid.med <- median(resid)
    ## resid.mad <- mad(resid)
    ## limits <- resid.med+resid.mad * dev* c(-1,1)
    ## score[!NAs] <- abs(pmin((resid-limits[1]),0) + pmax((resid - limits[2]),0))
    return(score)
}


#' Function to detect outliers in time serires
#' @param x vector of time series
#' @param t vector of time stamps
#' @param breaks break points to re-smooth the data
#' @param dev deviation (dev * mad) from the median to be considered as outliers
tsoutliers.w.brks <- function(x,t,breaks,dev=2){
    idList <- splitAt(1:length(x) , breaks)
    score <- foreach(id=idList, .combine = 'c') %do%{
        tsoutliers(x[id],t[id],dev)
    }    
    return(score)
}


#' Function to convert parameter vector alpha and matrix beta to MDSINE's output format
#' @param alpha
#' @param beta
#' @param gamma
#' @param vnames variable names
#' @return a data frame in MDSINE's output format
formatOutput <- function(alpha, beta, gamma=NULL, vnames){
    suppressMessages(library(reshape2))
    p <- length(alpha)
    nPerturbs <- ncol(gamma)
    if(is.null(gamma)) nPerturbs <- 0       
    parm <- data.frame(parameter_type=c(rep("growth_rate",p), rep("interaction",p*p), rep("perturbation", nPerturbs*p)))
    if(nPerturbs>0){
        tmp <- rep(paste0("perturbation",1:nPerturbs), each=p)
    }else{
        tmp <- NULL
    }
    parm$source_taxon <- c(rep(NA,p),rep(vnames, each = p),tmp)
    parm$target_taxon <- vnames
    parm$value <- 0
    parm$significance <- NA
    parm$MCMC_std <- NA
    tmp.beta <- melt(t(beta))
    tmp.beta <- tmp.beta[order(tmp.beta$Var1),]
    parm$value <- c(alpha, tmp.beta$value, c(gamma))
    parm
}

#' Function to estimate parameters using gradient matching bayesian lasso 
#' @param X matrix of data (variables in columns, measurements in rows)
#' @param P matrix of perturbation indicators
#' @param Ys matrix of response (variables in columns, measurements in rows)
#' @param Fs matrix of filter indicating whether the Y-X pair is used
#' @param ncpu number of CPU used
#' @param rmSp the index of variable removed for alr transformed data
#' @param vnames the name of species to format output
#' @param seed the seed
#' @return a list with alpha, beta, and output in MDSINE's format
BLASSO <- function(X, P, Ys, Fs, ncpu, rmSp, vnames, seed=NULL){
    suppressMessages(library(monomvn))
    registerDoMC(ncpu)    
    p <- ncol(X)
    nPerturbs <- ncol(P)
    X[is.na(X)] <- 0
    tmp <- ifelse(rmSp==0, p, p-1)
    Fs[is.na(Fs)] <- TRUE
    theta <- foreach(o = 1:tmp, .combine=cbind) %dopar% {
        Y <- Ys[, o]
        f <- Fs[, o]
        if(!is.null(seed)) set.seed(seed)	
        if(is.null(P)){
            bl.fit <- blasso(X[!f & !is.na(Y),], Y[!f & !is.na(Y)], verb = 0, T=1500)
        }else{
            bl.fit <- blasso(cbind(P[!f & !is.na(Y),],X[!f & !is.na(Y),]), Y[!f & !is.na(Y)], verb = 0, T=1500)
        }
        a <- mean(bl.fit$mu[-(1:500)])
        b <- apply(bl.fit$beta[-c(1:500),],2,mean)
        c(a,b)
    }
    alpha.v <- rep(0, p)
    beta.m <- matrix(0, p,p)
    if(rmSp==0){
        alpha.v <- theta[1,]
    }else{
        alpha.v[-rmSp] <- theta[1,]
    }
    if(is.null(nPerturbs)){
        gamma.m <- NULL
        if(rmSp==0){
            beta.m <- t(theta[-1,])
        }else{
            beta.m[-rmSp,] <- t(theta[-1,])
        }
    }else{
        if(rmSp==0){
            gamma.m <- matrix(0, nrow=p, ncol=nPerturbs)
            gamma.m <- t(theta[-1,][1:nPerturbs,])
            beta.m <- t(theta[-1,][-(1:nPerturbs),])
            gamma.m <- postProcessGamma(alpha.v, gamma.m)            
        }else{
            gamma.m <- matrix(0, nrow=p, ncol=nPerturbs)
            gamma.m[-rmSp,] <- t(theta[-1,][1:nPerturbs,])
            beta.m[-rmSp,] <- t(theta[-1,][-(1:nPerturbs),])
            gamma.m <- postProcessGamma(alpha.v, gamma.m)
        }
    }    
    beta.m <- postProcessBeta(beta.m)

    return(list(alpha.v=alpha.v, beta.m=beta.m, gamma.m=gamma.m,
                mdsine=formatOutput(alpha.v, beta.m, gamma.m, vnames)))
}

#' Function to estimate parameters using gradient matching bayesian variable select
#' @param X matrix of data (variables in columns, measurements in rows)
#' @param P matrix of perturbation indicators
#' @param Ys matrix of response (variables in columns, measurements in rows)
#' @param Fs matrix of filter indicating whether the Y-X pair is used
#' @param vnames the name of species to format output
#' @param selBest select the best model (default: FALSE, i.e. use the best model where the intercept is positive)
#' @param ncpu number of CPU used
#' @return a list with alpha, beta, and output in MDSINE's format
BVS <- function(X, P, Ys, Fs, vnames=colnames(X), selBest=FALSE, ncpu=1){
    suppressMessages(library(BayesVarSel))
    p <- ncol(X)
    tmpnames <- paste0("var",1:p)
    colnames(X) <- tmpnames
    nPerturbs <- ncol(P)
    X[is.na(X)] <- 0
    Fs[is.na(Fs)] <- TRUE    
    theta <- foreach(o = 1:p, .combine=cbind) %do% {
        Y <- Ys[, o]
        f <- Fs[, o]
        if(is.null(P)){
            sink("/dev/null")
            bvs.fit <- Bvs("Y~.", fixed.cov = c("Intercept",tmpnames[o]),
                           data=data.frame(X[!f,],Y=Y[!f]),n.keep=100)
            sink()
        }else{
            sink("/dev/null")
            bvs.fit <- Bvs("Y~.", fixed.cov = c("Intercept",tmpnames[o]),
                            data=data.frame(P[!f,],X[!f,],Y=Y[!f]),n.keep=100)
            sink()
        }
        apply(BMAcoeff.new(bvs.fit,ncpu=ncpu,best=selBest),2,mean)        
    }    
    alpha.v <- rep(0, p)
    beta.m <- matrix(0, p,p)
    alpha.v <- theta[1,]
    if(is.null(nPerturbs)){
        gamma.m <- NULL
        beta.m <- t(theta[-1,])
    }else{
        gamma.m <- matrix(0, nrow=p, ncol=nPerturbs)
        gamma.m <- t(theta[-1,][1:nPerturbs,])
        beta.m <- t(theta[-1,][-(1:nPerturbs),])
        gamma.m <- postProcessGamma(alpha.v, gamma.m)            
    }
    return(list(alpha.v=alpha.v, beta.m=beta.m, gamma.m=gamma.m,
                mdsine=formatOutput(alpha.v, beta.m, gamma.m, vnames)))
}



#' Function adapted from package "BayesVarSel". Plots are supressed by default
#' @param x BVS object from Bvs function
#' @param ncpu number of CPU used
#' @param n.sim number of simulation per model (default: 5000)
#' @param method method for sampling (default: svd)
#' @param best use the best model (default: FALSE, i.e. choose a model with positive intercept)
BMAcoeff.new <- function (x, ncpu, n.sim = 5000, method = "svd", best=FALSE) {
    registerDoMC(ncpu)    
    if (!is.null(x$lmnull)) {
        if (colnames(x$lmnull$x)[1] == "(Intercept)") {
            colnames(x$lmnull$x)[1] <- "Intercept"
        }
    }
    if (colnames(x$lmfull$x)[1] == "(Intercept)") {
        colnames(x$lmfull$x)[1] <- "Intercept"
    }
    
    n.model <- dim(x$modelsprob)[1]
    name.y <- colnames(x$lmfull$model[1])
    bma.coeffs <- matrix(0, nrow = n.sim, ncol = length(x$lmfull$coefficients))
    colnames(bma.coeffs) <- colnames(x$lmfull$x)
    models <- rep(1:n.model, n.sim)    
    t.models <- table(models)
    cs.tmodels <- cumsum(t.models)
    X <- x$lmfull$x
    
    regulation <- foreach (iter=1:length(t.models), .combine=c) %dopar% {
        rMD <- as.numeric(names(t.models)[iter])
        covsrMD <- names(x$modelsprob[rMD, ])[x$modelsprob[rMD, ] == "*"]
        datarMD <- as.data.frame(cbind(x$lmfull$model[, name.y], X[, covsrMD]))
        colnames(datarMD) <- c(name.y, covsrMD)
        if (!is.null(x$lmnull)) {
            datarMD <- cbind(datarMD, x$lmnull$x)
        }
        colnames(datarMD) <- gsub("`", "", colnames(datarMD))
        formMD <- as.formula(paste(name.y, "~.-1", sep = ""))
        fitrMD <- lm(formula = formMD, data = as.data.frame(datarMD), qr = TRUE)
        prod(fitrMD$coefficients)
    }
    if(best || all(regulation<0)){
        idx <- 1
    }else{
        idx <- which(regulation>0)[1]
    }    
    
    rMD <- as.numeric(names(t.models)[idx])
    howmany <- t.models[idx]
    covsrMD <- names(x$modelsprob[rMD, ])[x$modelsprob[rMD, ] == "*"]
    datarMD <- as.data.frame(cbind(x$lmfull$model[, name.y], X[, covsrMD]))
    colnames(datarMD) <- c(name.y, covsrMD)
    if (!is.null(x$lmnull)) {
        datarMD <- cbind(datarMD, x$lmnull$x)
    }
    colnames(datarMD) <- gsub("`", "", colnames(datarMD))
    formMD <- as.formula(paste(name.y, "~.-1", sep = ""))
    fitrMD <- lm(formula = formMD, data = as.data.frame(datarMD), qr = TRUE)
    Rinv <- qr.solve(qr.R(fitrMD$qr))
    iXtX <- Rinv %*% t(Rinv)
    Sigma <- sum(fitrMD$residuals * datarMD[, name.y]) * iXtX/fitrMD$df
    rcoeff <- rmvt(n = howmany, sigma = Sigma, df = fitrMD$df, 
                   delta = fitrMD$coefficients, type = "shifted", 
                   method = method)
    bma.coeffs[,names(fitrMD$coefficients)] <- rcoeff
    
    class(bma.coeffs) <- "bma.coeffs"
    return(bma.coeffs)
}



#' Function to perform cumulative sum scaling (CSS)
#' @param m matrix of data (variables in columns, measurements in rows)
#' @param p quantile used for normalization (default: 0.5)
cumSumScale <- function(m, p=0.5){
    m[m==0] <- NA
    ## find the quantile in each sample
    quant <- apply(m,1,function(x) quantile(x, p=p, na.rm = TRUE))
    ## calculate normalization factor
    f <- rowSums(m*sweep(m, 1, quant, '<='),na.rm = TRUE)
    nf <- f/exp(mean(log(f)))
    dat.css <- sweep(m,1,nf, '/')
    return(list("normCounts" = dat.css, "normFactors"=nf))
}

#' Function to preprocess input data
#' Input counts are normalized with CSS within each subject across all time points 
#' and then normalized across subjects with method described by David et. al.
#' @param counts input data following MDSINE's OTU table (i.e. variables in rows and samples in columns)
#' @param metadata metadata following MDSINE's metadata format
#' @param rsp reference species used for alr transformation
#' @param dev deviation (dev * mad) from the median to be considered as outliers
#' @param scaling scale the normalized counts (default 5000)
#' @param smooth_data data will be smoothed before initial normalization (default: TRUE)
#' @param forceBreak force to break the trajectory to handle pulsed perturbation (or species invasion) (default: NULL)
#' @param finite_diff use finite difference method to calculate gradients (default FALSE)
#' @param ncpu number of CPUs (default: 10)
preProcess <- function(counts, metadata, rsp, dev=100, scaling=5000, smooth_data=TRUE, ncpu=10, finite_diff=FALSE, forceBreak=NULL){    
    ## get perturbation identifier matrix (mu)
    if(any(metadata[,5]!=0)){
        mu <- model.matrix(~factor(metadata[,5]))[,-1]
        if(is.null(ncol(mu))) mu <- matrix(mu, ncol=1)
    }else{mu = NULL}
    ## normalize counts into relative abundances
    dat.norm <- apply(counts, 2, function(x) x/sum(x, na.rm = TRUE))
    ## compute alr transformed data
    dat.alr <- alr.transformation(t(dat.norm), rsp)
    if(smooth_data && !finite_diff){
        ## data smoothing    
        dat.alr.smoothed <- foreach(i=unique(metadata$subjectID), .combine='rbind') %do%{
            sel <- metadata$subjectID == i
            brk <- unique(sort(c(which(diff(metadata[sel,5]) !=0) + 1,
                                 which(forceBreak[sel] == 1))))
            SmoothGradient(dat.alr[sel, ], metadata[sel, 4], breaks=brk,
                           deriv=0, ncpu=10, method='pspline',smethod=3) 
        }
        dat.alr.smoothed[is.na(dat.alr)] <- NA
        ## transform back
        dat.norm.smoothed <- t(alr.transformation.inv(dat.alr.smoothed, rsp))
        dat.norm <- dat.norm.smoothed
    }
    ## CSS normalization across all the samples
    dat.norm <- t(cumSumScale(t(dat.norm))$normCounts)
    ## compute gradients
    dat.norm.scale <- data.frame(matrix(NA, nrow(counts),0))
    isOutlier <- NULL
    dalr_x_dt <- data.frame(matrix(NA, nrow(counts)-1,0))    
    for(i in unique(metadata$subjectID)){
        sel <- metadata$subjectID==i
        dat.tmp.norm <- dat.norm[, sel]
        ## compute addative log ratio transformation
        alr.tmp <- dat.alr[sel,]
        brk <- unique(sort(c(which(diff(metadata[sel,5]) !=0) + 1,
                             which(forceBreak[sel] == 1))))
        if(finite_diff){
            dalr_x_dt.tmp <- t(rbind(diff(alr.tmp[,-rsp])/diff(metadata[sel,4]),NA))
        }else{
            dalr_x_dt.tmp <- t(SmoothGradient(alr.tmp[,-rsp],
                                              metadata[sel, 4], breaks=brk,
                                              method = "pspline", ncpu=10, smethod=3))
        }
        dalr_x_dt <- cbind(dalr_x_dt, dalr_x_dt.tmp)
        ## detect outliers
        if (!finite_diff){
            outliers.tmp <- t(apply(dalr_x_dt.tmp, 1, function(x)
                tsoutliers.w.brks(x,metadata[sel,4],brk,dev=dev) > 0))
        }else{
            outliers.tmp <- matrix(FALSE, nrow(dalr_x_dt.tmp), ncol(dalr_x_dt.tmp))
        }
        isOutlier <- cbind(isOutlier, outliers.tmp)
        ## scaling
        dat.tmp.norm <- dat.tmp.norm * scaling ##* sf[subj.id==i]
        dat.norm.scale <- cbind(dat.norm.scale, dat.tmp.norm)
    }
    if(finite_diff){
        ## detect outliers for all gradients
        isOutlier <- t(apply(dalr_x_dt,1, function(x)
            abs(x-median(x, na.rm = TRUE))/mad(x, na.rm = TRUE) > dev ))
    }
    list(normData=dat.norm.scale, perturbInd=mu, alrGradient=dalr_x_dt, isOutlier=isOutlier)    
}

#' Function to filter small interactions
#' small interaction is defined relative to the mean of the diagonal of the matrix
#' @param beta interaction matrix estimated with BLASSO
#' @param thre_pos threshold used to remove small positive interactions (default 1e-3 * slef interaction)
#' @param thre_neg threshold used to remove small negative interactions (default 1e-3 * slef interaction)
postProcessBeta <- function(beta, thre_pos=1e-3, thre_neg=1e-3){
    tmp <- abs(diag(beta))
    ##scale <- median(tmp[tmp!=0])
    f <- foreach(r=1:nrow(beta), .combine = 'rbind') %do%{
        cr1 <- beta[r,]>0 & abs(beta[r,]) < (tmp[r] * thre_pos)
        cr2 <- beta[r,]<0 & abs(beta[r,]) < (tmp[r] * thre_neg)
        res <- cr1 | cr2
        res[r] <- FALSE
        res
    }
    beta[f] <- 0
    beta
}

#' Function to filter small perturbation effect
#' small interaction is defined relative to the mean of the diagonal of the matrix
#' @param alpha growth vector estimated with BLASSO
#' @param gamma perterbation effect matrix estimated with BLASSO
#' @param thre threshold used to remove small perturbation effect (default 1e-2 * self interactions)
postProcessGamma <- function(alpha, gamma, thre=1e-2){
    f <- min(abs(alpha)) * thre    
    gamma * sweep(abs(gamma), 1, f , '>')    
}

#' Function to calculate biomass for each sample
#' @param tss matrix of relative abundances (proportions), variables in rows and sample in columns
#' @param gradients matrix of gradients of addative log (alr) transformed abundances
#' @param perturbInd matrix of perturbation indicator generated by preProcess
#' @param meta metadata
#' @param rmSp Speices removed for alr tranformation
#' @param params list(alpha=, beta=) estimated with BLASSO
#' @param ncpu number of CPUs used (default: 10)
#' @param scale scale the median of all samples
#' @param smooth smooth the biomass after normalization
#' @param forceBreak force to break the trajectory to handle pulsed perturbation (or species invasion) (default: NULL)
NORM <- function(tss, gradients, perturbInd, metadata, rmSp, params, ncpu=10, scale=NA, smooth=FALSE, forceBreak=NULL){
    registerDoMC(ncpu)
    tss[is.na(tss)] <- 0    
    if(is.null(params$gamma)){
        perturbOffset <- 0
    }else{
        perturbOffset <- (params$gamma %*% perturbInd)[-rmSp,]
    }    
    gradients[is.na(gradients)] <- 0
    Ys <- gradients - params$alpha[-rmSp] - perturbOffset
    Xs <- params$beta[-rmSp,] %*% tss
    

    comb <- function(x, ...) {
        ## adapted from http://stackoverflow.com/questions/19791609/saving-multiple-outputs-of-foreach-dopar-loop
        ## might have an efficient implementation        
        lapply(seq_along(x),
               function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
    }
    
    w.list <- foreach (t = 1:ncol(tss), .combine='comb', .multicombine=TRUE) %dopar%{        
        Y <- Ys[,t]
        X <- Xs[,t]
        ## remove points that have different signs
        wrong.sign <- sign(X) != sign(Y)        
        outliers <- abs(Y - median(Y)) > 1.5 * IQR(Y) | abs(X - median(X)) > 1.5 * IQR(X) | wrong.sign
        if(all(outliers)) outliers <- wrong.sign
        if(all(wrong.sign)) {
            w.tmp <- scale
            warning(paste0("Sample @", t,' has a negative biomass estimation!'))
        }else{
            ## OLS
            model <- lm(Y[!outliers]~X[!outliers]-1)
            w.tmp <- abs(model$coeff[1])
        }
        mse <- (X*w.tmp-Y)^2 
        names(w.tmp) <- NULL
        list(w.tmp, mse)
    }

    mse.mat <- (matrix(unlist((w.list[[2]])), p-1)) ##* apply(abs(Xs),2,function(x)x/sum(x))
    w <- unlist(w.list[[1]])

    mse <- mean(c(mse.mat), na.rm = TRUE)
    mse.weighted <- mean(apply(mse.mat,2,mean)*1/w, na.rm = TRUE)

    if(smooth){
        w <- foreach(i=unique(metadata$subjectID), .combine="c") %dopar%{
            sel <- metadata$subjectID==i
            brk <- unique(sort(c(which(diff(metadata[sel,5]) !=0) + 1,
                                 which(forceBreak[sel] == 1))))            
            tmp <- abs(SmoothGradient(matrix(w[sel], ncol=1), metadata[sel, 4],
                               breaks = brk, deriv=0, smethod=1))            
        }
    }
    
    if(!is.na(scale)){ 
        w <- w/median(w) * scale
    }
    normCounts <- t(sweep(tss, 2, w, '*'))
    list(biomass=w, normCounts=normCounts, mse=mse, mse.weighted=mse.weighted)
}

#' Function to determine whether to stop the EM
#' @param x a vector of mse trace
#' @param epsilon tolerance of mse
#' @return a boolean value about the termination
testStop <- function(x, epsilon=0.1){
    ## if(which.min(x)==length(x)) {return(FALSE)} ## MSE is still decreasing
    min_mse <- min(x)
    if(x[length(x)]-min_mse> min_mse * epsilon){return(TRUE)}    
    return(FALSE)
}

#' Function to select references for ALR
#' @param dat input data in relative abundances
#' @param meta metadata following MDSINE's metadata format
#' @param scaling library size to scale the data to proportions
suggestRefs <- function(dat, meta, scaling=1){
    ##dat <- apply(dat,2,function(x)x/sum(x))
    sps <- rownames(dat)
    ##1. remove ref with 0 values
    message("The following species are not recommended due to 0 values:")
    fil1 <- rowSums(dat==0)/ncol(dat) > 0.05
    message(paste0(sps[fil1], collapse=', '))
    ##2. filtering two tails 
    tmp <- rowMeans(apply(dat,2, function(x)x/sum(x)))
    fil2 <- tmp<0.05
    fil3 <- tmp==max(tmp)
    message("The following species are not recommended due to their low/high abudances:")
    message(paste0(sps[fil2|fil3], collapse=', '))
    cv <- foreach(s=unique(meta$subjectID), .combine=rbind) %do% {
        apply(dat[, meta$subjectID==s], 1, function(x) sd(x)/mean(x))
    }
    if(!is.null(dim(cv))){
        cv <- apply(cv,2,mean)
    }
    message("The following species is recommended as the reference:")
    message(names(sort(cv[!(fil1 | fil2 | fil3)])[1]))
    data.frame(cv=cv, index=1:length(cv), hasZero=fil1, isTooHigh=fil3, isTooLow=fil2)[order(cv),]
}

#' Function to estimate biomass and parameters simultaneously
#' @param dat input data following MDSINE's OTU table (i.e. variables in rows and samples in columns)
#' @param meta metadata following MDSINE's metadata format
#' @param forceBreak force to break the trajectory to handle pulsed perturbation (or species invasion) (default: NULL)
#' @param useSpline use spline to smooth data and calculate gradients (default: TRUE)
#' @param dev deviation (dev * mad) from the median to be considered as outliers (default:Inf, no filtering)
#' @param verbose print more information (default: TRUE)
#' @param refSp reference OTU for addative log ratio transformation (default: selected by BEEM)
#' @param max_iter maximal number of iterations for the EM algorithm (default: 100)
#' @param min_iter minimal number of iterations for the EM algorithm (default: 30)
#' @param epsilon tolerance threshold to early stop the iterations (default: 0.01)
#' @param ncpu maximal number of CPUs used (default:10)
#' @param seed seed used in BLASSO (default:NULL)
#' @param scaling median total biomass to scale all biomass data (default:1000)
EM <- function(dat, meta, forceBreak=NULL, useSpline=TRUE,
               dev=Inf, verbose=TRUE,
               refSp = NULL,
               min_iter = 30,  max_iter = 100, epsilon=0.01,
               ncpu=10, seed=NULL, scaling=1000){

    if(nrow(dat) < 7){
        warning("There are less than 7 species. This might results in an inaccurate model.")
    }    
    refRank <- suggestRefs(dat, meta)
    if(is.null(refSp)){
        message("No input for reference species, selecting one with the lowest coefficient of variation...")
        refSp <- refRank$index[1]
    }
    if(sum(dat[refSp,]==0)>0){
        message("[!]: The reference species has zero abundance in some samples. This will treated as non-zeros by adding a pseudo count.")
    }
    if(refRank[refSp, 1]>0.9) {
        message("[!]: The reference species has high CV (>90%). Parameter estiamtes might be inaccurate (check the trace of weighted mse for convergence).")
    }
    message(paste0("Reference species: ",  rownames(dat)[refSp]))
    
    p <- nrow(dat)
    dat.tss <- apply(dat, 2, function(x)x/sum(x, na.rm = TRUE))
    if(verbose) message("Preprocessing data ...")
    tmpPreProcessed <- preProcess(dat, meta, rsp=refSp, dev=dev, scaling=scaling, smooth_data = TRUE, finite_diff=!useSpline, ncpu=ncpu, forceBreak=forceBreak)
    gradients.T <- tmpPreProcessed$alrGradient
    gradients <- t(gradients.T)
    perturbInd <- tmpPreProcessed$perturbInd
    isOutlier <- t(tmpPreProcessed$isOutlier)
    ##print(sum(c(isOutlier), na.rm = TRUE))
    ## if we know the species abundance should be 0, we force them to be outliers
    isOutlier[t(is.na(dat))[,-refSp]] <- TRUE
    
    ## initialization
    SpNames <- rownames(dat)
    dat.iter <- t(tmpPreProcessed$normData)
    biomass.traj <- data.frame("0"=rowSums(dat.iter))
    mse.traj <- c(Inf)
    mse.weighted.traj <- c(Inf)

    parm.traj <- apply(expand.grid(SpNames, SpNames), 1,
                       function(x) paste0(x[2], '->', x[1]))
    parm.traj <- c(paste0(NA, '->',SpNames), parm.traj)
    for(i in unique(meta[,5])){
        if(i!=0){
            parm.traj <- c(parm.traj, paste0("perturb" ,i, '->',SpNames))
        }            
    }
    parm.traj <- data.frame(name=parm.traj)

    for(iter in 1:max_iter){
        if(verbose){
            message(paste0("##Iteration ", iter))
            message("####solve alpha and beta (E step)")
        }
        
        parms <- BLASSO(X=dat.iter,P=perturbInd,Ys=gradients,
                        Fs=isOutlier, ncpu=ncpu, rmSp=refSp, vnames=SpNames, seed=seed)
        
        alpha.v = parms$alpha.v
        beta.m <- parms$beta.m
        gamma.m <- parms$gamma.m
        mdsine <- parms$mdsine
        tmp <- mdsine$value
        parm.traj <- cbind(parm.traj, tmp)

        ## test for termination
        if(iter > min_iter && testStop(mse.traj, epsilon)) break
        if(verbose) message( "####normalize (M step)####")
        
        tmp <- NORM(dat.tss, gradients.T, t(perturbInd), meta, refSp, list(alpha=alpha.v, beta=beta.m, gamma=gamma.m), ncpu=ncpu, scale=scaling, smooth=useSpline, forceBreak = forceBreak)
        biomass.traj <- cbind(biomass.traj, tmp$biomass)
        mse.traj <- c(mse.traj, tmp$mse)
        mse.weighted.traj <- c(mse.weighted.traj, tmp$mse.weighted)
        normCounts <- tmp$normCounts
        if(verbose) message(paste0("Weighted mse: ", tmp$mse.weighted))
        dat.iter <- normCounts
        if(verbose) message( "##########################")
        
    }
    
    return(list(
        final.params=mdsine,
        trace.biomass=biomass.traj,
        trace.params=parm.traj,
        trace.mse=mse.traj,
        trace.mse.weighted=mse.weighted.traj,
        min.iter=min_iter,
        max.iter=max_iter,
        epsilon=epsilon,
        refSp=refSp
        ) )
}


#' Independent function for inferring parameters with Bayesian Lasso
#' @param dat input data following MDSINE's OTU table
#' @param metadata metadata following MDSINE's metadata format
#' @param biomass biomass data following MDSINE's biomass data format
#' @param forceBreak force to break the trajectory to handle pulsed perturbation (or species invasion) (default: NULL)
#' @param dev deviation (dev * mad) from the median to be considered as outliers (default:Inf, no filtering)
#' @param ncpu maximal number of CPUs used (default:10)
#' @param infer_flag run inference (default:TRUE)
param.infer <- function(dat, metadata, biomass,
                        forceBreak=NULL, dev=Inf, ncpu=10, infer_flag=TRUE){
    registerDoMC(ncpu)
    log.transform <- function(x){
        tmp <- log(x)
        tmp[!is.finite(tmp)] <- 0
        tmp
    }
    dat[is.na(dat)] <- 0
    counts.tss <- t(apply(dat, 2, function(x)x/sum(x)))
    ## smooth to calculate gradient for log relative abundances
    dlnx_tilde_dt <- foreach(i=unique(metadata$subjectID), .combine='rbind') %do%{
        sel <- metadata$subjectID == i
        brk <- unique(sort(c(which(diff(metadata[sel,5]) !=0) + 1,
                             which(forceBreak[sel] == 1))))            
        SmoothGradient(log.transform(counts.tss)[sel, ], metadata[sel, 4],
                       breaks=brk,method='pspline',smethod=3, deriv=1, ncpu=ncpu)
    }    
    ## smooth to calculate gradient for log biomass
    dlnm_dt <- foreach(i=unique(metadata$subjectID), .combine='c') %do%{
        sel <- metadata$subjectID == i
        brk <- unique(sort(c(which(diff(metadata[sel,5]) !=0) + 1,
                             which(forceBreak[sel] == 1))))            
        SmoothGradient(matrix(log.transform(biomass)[sel],ncol=1), metadata[sel, 4],
                       breaks=brk,method='pspline',smethod=3, deriv=1, ncpu=ncpu)
    }
    Ys <- dlnx_tilde_dt + dlnm_dt   
    Xs <- exp(sweep(log.transform(counts.tss), 1, log.transform(biomass), '+'))
    if(infer_flag == FALSE) return(list(Ys=Ys, Xs=Xs))
    ## detect outliers
    isOutlier <- foreach(i=unique(metadata$subjectID), .combine='rbind') %dopar%{
        sel <- metadata$subjectID==i
        brk <- unique(sort(c(which(diff(metadata[sel,5]) !=0) + 1,
                             which(forceBreak[sel] == 1))))
        apply(Ys[sel,], 2, function(x)
            tsoutliers.w.brks(x,metadata[sel,4],brk,dev=dev) > 0)
    }
    BLASSO(X=Xs, P=NULL, Ys=Ys, Fs=isOutlier, ncpu=ncpu, rmSp=0, vnames=rownames(dat))
}

#' Inferring biomass from BEEM results
#' @param beem.obj BEEM output list
biomassFromEM <- function(beem.obj){
    trace.mse <- beem.obj$trace.mse
    min.mse <- min(trace.mse)
    em.idx <- which((trace.mse-min.mse) < beem.obj$epsilon*min.mse)
    biomass <- apply(beem.obj$trace.biomass[,em.idx],1,median)
    biomass
}

#' Inferring parameters from BEEM results
#' @param beem.obj BEEM output list
#' @param counts counts data following MDSINE's OTU table
#' @param metadata metadata following MDSINE's metadata format
#' @param forceBreak force to break the trajectory to handle pulsed perturbation (or species invasion) (default: NULL)
#' @param ncpu maximal number of CPUs used (default:10)
#' @param enforceLogistic re-estimate the self-interaction parameters (enforce to negative values)
paramFromEM <- function(beem.obj, counts, metadata, forceBreak=NULL, ncpu=10, enforceLogistic=FALSE){
    if(NROW(counts) <7){
        warning('You have less than 7 species. The estimation of parameters might be inaccurate.')
    }
    registerDoMC(ncpu)
    trace.mse <- beem.obj$trace.mse
    min.mse <- min(trace.mse)
    em.idx <- which((trace.mse-min.mse) < beem.obj$epsilon*min.mse)
    ##em.idx <- which((trace.mse-min.mse) < 0.01*min.mse)

    refSp <- beem.obj$refSp
    if(length(em.idx)==1) {
        beem.biomass <- beem.obj$trace.biomass[,em.idx]
        beem.param <- beem.obj$trace.params[,em.idx]
    }else{
        beem.biomass <- apply(beem.obj$trace.biomass[,em.idx],1,median)
        beem.param <- apply(beem.obj$trace.params[,em.idx],1,median)
    }
    p <- nrow(counts)
    ## solve for interaction matrix
    beem.a <- beem.param[1:p]
    beem.b <- matrix(beem.param[-(1:p)],p,p)    
    beem.b.median <- apply(beem.b[-refSp,],2,median)
    beem.beta <- t(t(beem.b) - beem.b.median)

    ## significance
    tmp <- beem.beta
    diag(tmp) <- NA
    beem.sig <- abs(beem.beta)/sd(tmp,na.rm = T)
    ## growth rate for reference species
    tmp <- param.infer(dat=counts, metadata=metadata, biomass=beem.biomass,
                       forceBreak=forceBreak, ncpu=ncpu, infer_flag=FALSE)    
    Xs <- tmp$Xs
    Ys <- tmp$Ys
    beem.alpha <- foreach(idx=1:p, .combine =c) %dopar%{
        tmp <- (Ys[,idx] - (beem.beta %*% t(Xs))[idx,])
        ## tmp <- tmp[(tmp-median(tmp)) <= 2*IQR(tmp)]
        median(tmp[tmp>0])
    }
    if(any(beem.alpha<0) || any(is.na(beem.alpha))) {
        message("Warning: Not enough time points to enforce positive growth rate.")
    }
    ## fix self interaction
    if (enforceLogistic){
        beem.diag <- foreach(idx=1:p, .combine =c) %dopar%{
            beem.tmp <- beem.beta
            diag(beem.tmp) <- 0
            tmp <- (Ys[,idx] - (beem.tmp %*% t(Xs))[idx,] - beem.alpha[idx])/Xs[,idx]
            tmp <- tmp[(tmp-median(tmp)) <= 2*IQR(tmp)]
            median(tmp[tmp<0])
        }
        diag(beem.beta) <- beem.diag
    }

    beem.MDSINE <- beem.obj$final.params
    beem.MDSINE$value <- c(beem.alpha, c(beem.beta))
    beem.MDSINE$significance <- c(rep(10000,p), c(beem.sig))
    if (enforceLogistic){
        beem.MDSINE$significance[beem.MDSINE$parameter_type=='interaction' & beem.MDSINE$source_taxon == beem.MDSINE$target_taxon ] <- rep(10000,p)
    }
    beem.MDSINE
}

