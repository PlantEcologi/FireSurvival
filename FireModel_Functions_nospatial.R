
##### Deviance Information Criterion
UpdateDbar <- function (Y, XBetaW,firerows){
    ## Calculate Dbar
    loglik = CalcLogLik(Y, XBetaW, gid, firerows)
    return(-2*loglik) # return Dbar
}
    
    
###### Function to update Z

UpdateZ<- function(Y, XBetaW, N) { # Y is the 1/0 response data; XBetaW is the regression data matrix * the regression coefficients + the spatial random effects,
                    # and N is the sample size (i.e. the length of Y)                 
    sampleZ={}
    U=runif(N)
    Z=rep(0,N)
    index = (Y==1)
    Z[!index]=qnorm(U[1:sum(!index)]*pnorm(XBetaW[!index]),XBetaW[!index])  
    Z[index]=qnorm(U[(sum(!index)+1):N]*(pnorm(XBetaW[index],lower.tail=F))+pnorm(XBetaW[index]),XBetaW[index])
    sampleZ=Z  
    return(sampleZ)
}

### Fire Process

fire0.update <- function (p, lastnoburn, nGrid, gidrow, leftcensorlength) {
    fire0temp = rep(0, nGrid)
    pmat = as.matrix(p, nrow=nGrid)
    for (i in 1:nGrid) {
        leftcensoredrows = gidrow[i]:(gidrow[i]+leftcensorlength-1)
        pinterval = rev(c(1, cumprod(rev(p[leftcensoredrows]))))
        fire0temp[i] = match(1, rmultinom(1, 1, pinterval))
    }    
    return(fire0temp)
}

#### Functions to calculate likelihood and posterior probability values

UpdateLeftExclude <- function(fire0, nGrid, nSeas) {
    leftexclude = matrix(F, nGrid, nSeas) # matrix of Trues and Falses that indicates which rows are to be excluded since they occurred before the current estimated first fire
    for (i in 1:nGrid) leftexclude[i,1:fire0[i]] = T
    leftexclude = as.vector(t(leftexclude))
    return(leftexclude)
}

CalcLogLik <- function(Y, X, Beta,leftexclude,firerows) {
    logp = pnorm(ifelse(Y == 1, -1, 1) * ((X%*%Beta) ), log.p = TRUE) #+ W[gid]
    logp[firerows+1] = 0  # set likelihood contribution to 0 (i.e p=1) for the seasons immediately following fires
    logp[leftexclude] = 0 # set likelihood contribution to 0 for seasons before the first fire
    loglik = sum(logp)
    return(loglik)
}


CalcPost <- function(Y, X, Beta, leftexclude, adj, num,firerows){
    loglik = CalcLogLik(Y, X, Beta, leftexclude,firerows)
#    Wpriormean = unlist(lapply(adj, fwmean, W))
#    Wpriorsd = sqrt(Sigma2_W/num)
    logpost = loglik + sum(
        ##Priors
        sum(dnorm(Beta, 0, 1000,log=T)),
        sum(dnorm(W, Wpriormean, Wpriorsd,log=T))
    )
    return(logpost)
}

CalcPostNonspatial <- function(Y, X, Beta, leftexclude,firerows){
    loglik = CalcLogLik(Y, X, Beta,leftexclude,firerows)
    logpost = loglik + sum(
        ##Priors
        sum(dnorm(Beta, 0, 1000,log=T))
    )
    return(logpost)
}


#### Generic MCMC functions

metrop.norm <- function(x0, propsd, f, params, this.par) {        # x0 = current value of parameter to be sampled
                                                        # propsd = proposal standard deviation
                                                        # f = function with value proportional to posterior probability function
                                                        # params = list of current values of all parameters supplied to f()
                                                        # this.par = name tag for parameter to be sampled among names(params)                                                  
    this.par = match.arg(this.par, names(params))
    pos = match(this.par, names(params)) # get position of current parameter in argument list
    K = length(x0)
    postcurrent = do.call(f, params)
    paramsnew = params
    propx = rnorm(K, x0, propsd)
    paramsnew[[pos]] = propx
    postnew = do.call(f, paramsnew)
    a = min(1, exp(postnew-postcurrent))
    if (runif(1) < a) { # test for acceptance
        return(propx)
    }
    return(x0)
}

adaptsd <- function(sdpar, Res_par) {
    n = length(Res_par)
    newsdpar = sdpar
    ind = 2:n
    accept = sum(Res_par[ind] != Res_par[ind-1]) / n
    if (accept > 0.3) newsdpar = newsdpar*1.3333
    if (accept < 0.15) newsdpar = newsdpar*0.8
    return(newsdpar)
}


###### W functions ############

### Function to calculate Adjacency Matrix
AdjMat= function(x, celldist)  {
                    #Where x is a data.frame with three columns, x, y, and gid.  Celldist is the distance (in xy units) between centroids.
    
    xd = (celldist)^2; # max squared x distance for neighborhood
    yd = (celldist)^2; # "" y distance ""
    md = xd + yd # max squared XY distance for neighborhood
    
    nc = nrow(x);
    num = vector(mode="numeric", length=nc);
    adj = {}
    adj=vector("list", nc)
    for (i in 1:nc) {
        z = ((x$x[i] - x$x)^2 + (x$y[i]-x$y)^2)
        ind = z<md 
        ind[z==0] = F # which cells are within squared distance md but not the same cell
        nn = sum(ind);
        num[i] = nn;
        adj[[i]] = x$gid[ind]
    }
    
    SumNumNeigh = sum(num);
    weights = vector(mode = "numeric", length = SumNumNeigh);
    weights[1:SumNumNeigh] = 1; # assign equal weights to neighbors
    
    return(list(adj=adj,num=num,SumNumNeigh=SumNumNeigh))
    
}

### Update the W's
UpdateW <- function(Z, XBeta, leftexclude, W, Sigma2_W, num, adj, ncells, nSeas, ecoRF, eco, firerows) {
    nonspatial = XBeta
    nonspatial[leftexclude] = NA # exclude the seasons up to first fire from the nonspatial prediction
    nonspatial[firerows+1] = NA # exclude the seasons immediately after fire from fitting the W, since, these are not in the likelihood either
    nonspatialmat = matrix(nonspatial, nrow = nSeas)
    munonspatial = apply(nonspatialmat, 2, mean, na.rm=T)
    mutemp = Z - munonspatial
    sigtemp= Sigma2_W/(num + Sigma2_W)
    muw={}
    for (i in 1:ncells) 
    { 
        muw[i] = sigtemp[i] * (mutemp[i] + sum(W[adj[[i]]])/Sigma2_W)  
        W[i] = rnorm(1,muw[i],sqrt(sigtemp[i])) 
    }
    W = W-mean(W) 
    return(W)
}

### Functions to update Sigma2_W
fwmean = function(d1,d2){sum(d2[d1])}

UpdateSigma2_W <- function(W, adj, num, a_Sigma2_W, b_Sigma2_W, N) {
    s1 = a_Sigma2_W + N/2
    s2 = b_Sigma2_W + num*t(W)%*%W/2 - t(W)%*%unlist(lapply(adj, fwmean, W))/2
    sampleSigma2_W = 1/rgamma(1, s1, s2)
    return( sampleSigma2_W)
}

### Function to update Sigma2_eco

UpdateSigma2_eco <- function(W, X, Beta ,Z, a_Sigma2_eco, b_Sigma2_eco, N, nSeas, ecoRF,eco) {
    residmat = matrix(Z-(X%*%Beta + W[gid] + ecoRF[eco]), nrow=nSeas)
    residtemp = apply(residmat, 2, mean)
    s1 = a_Sigma2_eco + N/2
    s2 = b_Sigma2_eco + (sum(residtemp)^2)/2
    sampleSigma2_eco = 1/rgamma(1, s1, s2)
    return( sampleSigma2_eco)
}
