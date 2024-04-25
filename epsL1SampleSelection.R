require(sampleSelection)


#############################################################################
# Data generation. We generated the errors from a bivariate normal distribution
# wih mean vector (0,0)and correlation matrix with rho. The correlation
# between variables is corr(xi, xj) = cor^|i-j|. We ensured this via
# "c(rev(cumprod(rep(rho,index-1))),1,cumprod(rep(rho,p-index)))"

##############################################################################
gendatss<-function(seed,n,cor){
  if (match("mvtnorm",.packages(),0)==0) require(mvtnorm)
  if (match("MASS",.packages(),0)==0) require(MASS)
  if (match("sampleSelection",.packages(),0)==0) require(sampleSelection)
  set.seed(seed)
  
  beta <- c(0.5,1, 1, 1.5, 0.2,0.2,0.1,0.1,0.5,1,1,0,0,0,0,0,0,0, 0, 0)
  gamma <- c(1.5,1, 1, 1,0.5,0.1,0.2,0.5,1.5,0.2,0.1,0,0,0,0,0 ,0,0,0,0, 1)
  rho <- 0.5
  p <- length(beta)
  Sig <- matrix(0, ncol=p, nrow=p)
  for(index in 1:p){
    Sig[index,] <-  c(rev(cumprod(rep(rho,index-1))),1,cumprod(rep(rho,p-index)))
  }
  X <- matrix(rnorm(p*n), ncol=p)
  X <- X %*% chol(Sig)
  meane<-c(0,0)
  cove <- matrix(c(1,cor,cor,1),nrow=2)
  e<-mvrnorm(n,meane,cove)
  
  bx <- cbind(1, X[,-20])
  ystar <- (bx%*%beta)+e[,1]
  gw <- cbind(1,X)
  ustar <- (gw%*%gamma)+e[,2]>0
  yobs <- ystar*(ustar>0)
  ymiss <- ifelse(yobs==0,NA,yobs)
  covar <- cov(X)
  dat<-data.frame(yobs,ustar,X,ymiss)
  return(dat)
}


########################################
#Likelihood Functions
########################################

loglik<-function(beta,YS,XS,YO,XO){
  if (match("MASS",.packages(),0)==0) require(MASS)
  options(digit=22)
  NXS <- ncol(XS)
  NXO <- ncol(XO)
  ibetaS <- 1:NXS
  ibetaO <- seq(tail(ibetaS, 1) + 1, length = NXO)
  isigma <- tail(ibetaO, 1) + 1
  irho <- tail(isigma, 1) + 1
  g <- beta[ibetaS]
  b <- beta[ibetaO]
  sigma <- beta[isigma]
  if (sigma < 0){sigma = 0.05}
  rho <- beta[irho]
  if(abs(rho>0.98)){rho = sign(rho)*0.95}
  XS.g = XS %*% g
  XO.b = XO %*% b
  u2<-YO-XO.b
  r <- sqrt(1 - rho^2)
  B <- (XS.g + rho/sigma * u2)/r
  ll <- ifelse(YS == 0, (pnorm(-XS.g, log.p = TRUE)), 
               -1/2 * log(2 * pi) - log(sigma) + (pnorm(B, log.p = TRUE) - 
                                                    0.5 * (u2/sigma)^2))
  return(-sum(ll))
  
}

#############################################
#Hessian and Gradient
#############################################
#Gradient
gradlik <- function(beta,YS,XS,YO,XO)
{
  NXS <- ncol(XS)
  NXO <- ncol(XO)
  ibetaS <- 1:NXS
  ibetaO <- seq(tail(ibetaS, 1) + 1, length = NXO)
  isigma <- tail(ibetaO, 1) + 1
  irho <- tail(isigma, 1) + 1
  nObs <- length(YS)
  nParam <- NXS + NXO + 2
  XS0 <- XS[YS == 0, , drop = FALSE]
  XS1 <- XS[YS == 1, , drop = FALSE]
  YO1 <- YO[YS == 1]
  XO1 <- XO[YS == 1, , drop = FALSE]
  N0 <- sum(YS == 0)
  N1 <- sum(YS == 1)
  w <- rep(1, N0 + N1)
  w0 <- rep(1, N0)
  w1 <- rep(1, N1)
  
  g <- beta[ibetaS]
  b <- beta[ibetaO]
  sigma <- beta[isigma]
  if (sigma < 0) 
    return(matrix(NA, nObs, nParam))
  rho <- beta[irho]
  if ((rho < -1) || (rho > 1)) 
    return(matrix(NA, nObs, nParam))
  XS0.g <- as.numeric(XS0 %*% g)
  XS1.g <- as.numeric(XS1 %*% g)
  XO1.b <- as.numeric(XO1 %*% b)
  u2 <- YO1 - XO1.b
  r <- sqrt(1 - rho^2)
  B <- (XS1.g + rho/sigma * u2)/r
  lambdaB <- exp(dnorm(B, log = TRUE) - pnorm(B, log.p = TRUE))
  gradient <- matrix(0, nObs, nParam)
  gradient[YS == 0, ibetaS] <- -w0 * XS0 * exp(dnorm(-XS0.g, 
                                                     log = TRUE) - pnorm(-XS0.g, log.p = TRUE))
  gradient[YS == 1, ibetaS] <- w1 * XS1 * lambdaB/r
  gradient[YS == 1, ibetaO] <- w1 * XO1 * (u2/sigma^2 - 
                                             lambdaB * rho/sigma/r)
  gradient[YS == 1, isigma] <- w1 * ((u2^2/sigma^3 - lambdaB * 
                                        rho * u2/sigma^2/r) - 1/sigma)
  gradient[YS == 1, irho] <- w1 * (lambdaB * (u2/sigma + 
                                                rho * XS1.g))/r^3
  grad <- apply(gradient,2,sum)
  return(-grad)
  
}



#Hessian
hesslik <- function(beta,YS,XS,YO,XO) 
{
  NXS <- ncol(XS)
  NXO <- ncol(XO)
  ibetaS <- 1:NXS
  ibetaO <- seq(tail(ibetaS, 1) + 1, length = NXO)
  isigma <- tail(ibetaO, 1) + 1
  irho <- tail(isigma, 1) + 1
  nObs <- length(YS)
  nParam <- NXS + NXO + 2
  XS0 <- XS[YS == 0, , drop = FALSE]
  XS1 <- XS[YS == 1, , drop = FALSE]
  YO1 <- YO[YS == 1]
  XO1 <- XO[YS == 1, , drop = FALSE]
  N0 <- sum(YS == 0)
  N1 <- sum(YS == 1)
  w <- rep(1, N0 + N1)
  w0 <- rep(1, N0)
  w1 <- rep(1, N1)
  
  g <- beta[ibetaS]
  b <- beta[ibetaO]
  sigma <- beta[isigma]
  if (sigma < 0) {
    return(matrix(NA, nrow = nParam, ncol = nParam))
  }
  rho <- beta[irho]
  if ((rho < -1) || (rho > 1)) {
    return(matrix(NA, nrow = nParam, ncol = nParam))
  }
  XS0.g <- as.vector(XS0 %*% g)
  XS1.g <- as.vector(XS1 %*% g)
  XO1.b <- as.vector(XO1 %*% b)
  u2 <- YO1 - XO1.b
  r <- sqrt(1 - rho^2)
  B <- (XS1.g + rho/sigma * u2)/r
  lambdaB <- exp(dnorm(B, log = TRUE) - pnorm(B, log.p = TRUE))
  C <- ifelse(B > -500, -exp(dnorm(B, log = TRUE) - pnorm(B, 
                                                          log.p = TRUE)) * B - exp(2 * (dnorm(B, log = TRUE) - 
                                                                                          pnorm(B, log.p = TRUE))), -1)
  hess <- matrix(0, nParam, nParam)
  a <- ifelse(XS0.g < 500, -exp(dnorm(-XS0.g, log = TRUE) - 
                                  pnorm(-XS0.g, log.p = TRUE)) * XS0.g + (exp(dnorm(-XS0.g, 
                                                                                    log = TRUE) - pnorm(-XS0.g, log.p = TRUE)))^2, 1)
  hess[ibetaS, ibetaS] <- -t(XS0) %*% (w0 * XS0 * a) + 
    t(XS1) %*% (w1 * XS1 * C)/r^2
  hess[ibetaS, ibetaO] <- -t(XS1) %*% (w1 * XO1 * C) * 
    rho/r^2/sigma
  hess[ibetaO, ibetaS] <- t(hess[ibetaS, ibetaO])
  hess[ibetaS, isigma] <- -rho/sigma^2/r^2 * t(XS1) %*% 
    (w1 * C * u2)
  hess[isigma, ibetaS] <- t(hess[ibetaS, isigma])
  hess[ibetaS, irho] <- t(XS1) %*% (w1 * (C * (u2/sigma + 
                                                 rho * XS1.g)/r^4 + lambdaB * rho/r^3))
  hess[irho, ibetaS] <- t(hess[ibetaS, irho])
  hess[ibetaO, ibetaO] <- t(XO1) %*% (w1 * (XO1 * ((rho/r)^2 * 
                                                     C - 1)))/sigma^2
  hess[ibetaO, isigma] <- t(XO1) %*% (w1 * (C * rho^2/sigma^3 * 
                                              u2/r^2 + rho/sigma^2 * lambdaB/r - 2 * u2/sigma^3))
  hess[isigma, ibetaO] <- t(hess[ibetaO, isigma])
  hess[ibetaO, irho] <- t(XO1) %*% (w1 * (-C * (u2/sigma + 
                                                  rho * XS1.g)/r^4 * rho - lambdaB/r^3))/sigma
  hess[irho, ibetaO] <- t(hess[ibetaO, irho])
  hess[isigma, isigma] <- sum(w1 * (-3 * u2 * u2/sigma^4 + 
                                      2 * lambdaB * u2/r * rho/sigma^3 + rho^2/sigma^4 * 
                                      u2 * u2/r^2 * C)) + sum(w1)/sigma^2
  hess[isigma, irho] <- hess[irho, isigma] <- -sum(w1 * 
                                                     (C * rho * (u2/sigma + rho * XS1.g)/r + lambdaB) * 
                                                     u2/sigma^2)/r^3
  hess[irho, irho] <- sum(w1 * (C * ((u2/sigma + rho * 
                                        XS1.g)/r^3)^2 + lambdaB * (XS1.g * (1 + 2 * rho^2) + 
                                                                     3 * rho * u2/sigma)/r^5))
  return(-hess)
  
}


############################################################################
#Penalty Approximation Hessian and Gradient
###########################################################################

penalty <- function(beta,lambda,eps,weightt){
  return(sum(lambda*weightt*sqrt(beta*beta+eps)))
}

pen_deriv <- function(beta,lambda,eps,weightt){
  return(weightt*((lambda*beta)/sqrt((beta*beta)+eps)))
}

pen_hess <- function(beta,lambda,eps,weightt){
  hessianvec <- lambda*(eps/(((beta*beta)+eps)^(3/2)))
  return(diag(hessianvec*weightt))
}



##########################################################################
#Line Search Method
#############################################################################


ss_pen_newton <- function(init,lambda,eps,YS,XS,YO,XO,weightt,max_iter=1000){
  beta_old <- init
  tol <- 1e-5
  converged <- FALSE
  iter <- 0
  like <- loglik(beta_old,YS,XS,YO,XO) + penalty(beta_old,lambda,eps,weightt)
  steplength=1
  while(!converged && iter < max_iter){
    hesspen <- hesslik(beta_old,YS=YS,XS=XS,YO=YO,XO=XO)+pen_hess(beta_old,lambda=lambda,eps=eps,weightt=weightt)
    gradpen <- gradlik(beta_old,YS=YS,XS=XS,YO=YO,XO=XO)+pen_deriv(beta_old,lambda=lambda,eps=eps,weightt=weightt)
    descent <- as.vector((solve(hesspen)%*%gradpen))
    steplength=1
    trial <- beta_old - steplength*descent
    while(abs(trial[length(trial)])>=0.98 || trial[length(trial)-1]<=0.02){
      steplength=steplength*0.9
      trial <- beta_old - steplength*descent
    }
    ntrials=1
    c1=10^(-8)
    suffdecrease = c1*sum(gradpen*descent)
    triallik <-loglik(trial,YS,XS,YO,XO) + penalty(trial,lambda,eps,weightt)
    stepfound = FALSE
    ntrials = 0
    while(stepfound == FALSE && ntrials < 20){
      if(triallik> like-suffdecrease*steplength){
        steplength=steplength*0.1
        trial <- beta_old - steplength*descent
        triallik <-loglik(trial,YS,XS,YO,XO) + penalty(trial,lambda,eps,weightt)
      }else{
          stepfound=TRUE
        }
      ntrials = ntrials + 1
    }
    like <- triallik
    beta_new <- beta_old - steplength*descent
    dif <- abs(beta_new-beta_old)
    if(max(abs(descent)) < tol){
      converged <- TRUE 
    }
    beta_old <- beta_new
    iter <- iter + 1
  }
  results = list(theta=beta_old,converged=converged)
  return(results)
}

threshold <- function(theta,eps=10^-8){
  for(i in 1:length(theta)){
    if(abs(theta[i])<(sqrt(eps))){theta[i]=0}
  }
  return(theta)
}
##########################################################
#Trust Region
##########################################################

require("trust")
trust_ss <- function(init,lambda,eps,YS,XS,YO,XO,weightt){
  objfun <- function(beta,...){
    val = loglik(beta,YS,XS,YO,XO) + penalty(beta,lambda,eps,weightt)
    grad <- gradlik(beta,YS=YS,XS=XS,YO=YO,XO=XO)+pen_deriv(beta,lambda=lambda,eps=eps,weightt=weightt)
    hess <- hesslik(beta,YS=YS,XS=XS,YO=YO,XO=XO)+pen_hess(beta,lambda=lambda,eps=eps,weightt=weightt)
    output <- list(value=val,
                   gradient = grad,
                   hessian = hess)
    return(output)
   }
    optim_out <- trust(objfun   = objfun,
                       parinit  = init,
                       rinit    = 1,
                       rmax     = 100000000,
                       iterlim  = 1000,
                       fterm    = 10^-10,
                       mterm    = 10^-10,
                       minimize = TRUE,
                       YS=YS,XS=XS,YO=YO,XO=XO,
                       lambda=lambda,eps=eps,weightt=weightt)
    trust_fit <- optim_out$argument
    for(i in 1:length(trust_fit)){
      if(abs(trust_fit[i])<0.0001){trust_fit[i]=0}
    }
    return(trust_fit)
}

############################################################
#BIC Function
###########################################################


bic_fun <- function(beta,YS,XS,YO,XO){
  n <- length(YO)
  df <- sum(beta!=0)-1
  nvar <- df+1
  
  bic <- loglik(beta,YS,XS,YO,XO)+nvar*log(n)
  return(bic)
}

############################################################
#Select lambda based on Newton-Raphson function using BIC
############################################################

lambda_bic_selection <- function(lambda_grid,init,YS,XS,YO,XO,weightt,method,nparout,nparsel,eps = 0.00000001){
  n <- length(YO)
  models <- matrix(NA,nrow=length(init),ncol=nrow(lambda_grid))
  models_bic <- rep(NA,times=nrow(lambda_grid))
  for(i in 1:nrow(lambda_grid)){
    lambda1 <- lambda_grid[i,1]
    lambda2 <- lambda_grid[i,2]
    lambda=c(rep(lambda1,nparout),rep(lambda2,nparsel),0,0)
    if(method=="Newton"){
      theta = ss_pen_newton(init,lambda,eps,YS,XS,YO,XO,weightt)
    }
    if(method=="Trust"){
      theta <- trust_ss(init,lambda,eps,YS,XS,YO,XO,weightt)
    }
    bic <- bic_fun(theta,YS,XS,YO,XO)
    models[,i]<- theta
    models_bic[i] <- bic
  }
  selected = which(models_bic == min(models_bic))
  result=list(coefficients=models[,selected],bic=models_bic[selected],lambda = lambda_grid[selected,])
  return(result)
}

####################################################################
#DE2D
######################################################################

de_lambda2d <- function(lambda_max,YS,XS,YO,XO,weightt,init,popsize,niter,nparout,nparsel,method="Newton",pcr=0.7,eps = 0.00000001){#
  library("lhs")
  print("Initialising")
  eps <- 0.00000001
  n <- length(YO)
  family <- seq(from=1,to=popsize)
  scale_factor <- 0.5
  lambda_de = lambda_max*maximinLHS(n=popsize,k=2) #initialization
  lambda_de = lambda_de[order(rowSums(lambda_de)^2,decreasing=F),]
  lambda_de = rbind(c(0,0),lambda_de)
  beta_mat <- matrix(NA,nrow=length(init),ncol=popsize+1)
  beta_mat[,1] = init
  bic_de <- rep(NA,popsize+1)
  converged_mat = rep(TRUE,popsize+1)
  bic_de[1]=bic_fun(beta=beta_mat[,1],YS,XS,YO,XO)
  for(m in 1:popsize){
    lambda=c(rep(lambda_de[m+1,1],nparout),rep(lambda_de[m+1,2],nparsel),0,0)
    if(method=="Newton"){
      newtoninit = beta_mat[,m]
      if(abs(newtoninit[length(newtoninit)])>0.95){newtoninit = init}
      if(abs(newtoninit[length(newtoninit)-1])<0.05){newtoninit = init}
      if(converged_mat[m]==FALSE){newtoninit=init}
      fit = ss_pen_newton(newtoninit,lambda,eps,YS,XS,YO,XO,weightt)
      theta_temp = fit$theta
      converged_mat[m+1] = fit$converged
    }
    if(method=="Trust"){
      theta_temp <- trust_ss(init,lambda,eps,YS,XS,YO,XO,weightt)
    }
    bic_de[m+1] <- bic_fun(beta=threshold(theta_temp),YS,XS,YO,XO)
    beta_mat[,m+1] <- theta_temp
  }
  print("Initialised")
  deiter = 1
  convergedde = FALSE
  while(convergedde==FALSE && deiter<=niter){
    print(deiter)
    for(j in 1:popsize+1){
      samp <- sample(family[-j],size=3)
      trial = lambda_de[j,]+(scale_factor*(lambda_de[samp[1],]-lambda_de[samp[2],]))
      prob=runif(2,0,1)
      offspring <- lambda_de[samp[3],]
      for(i in 1:2){
        if(trial[i]<0){trial[i]=0}
        if(trial[i]>lambda_max){trial[i]=lambda_max}
        if(prob[i]<0.7){offspring[i]=trial[i]}
      }
      lambda=c(rep(offspring[1],nparout),rep(offspring[2],nparsel),0,0)
      closest = which(rowMeans(sweep(lambda_de,2,offspring)^2)==min(rowMeans(sweep(lambda_de,2,offspring)^2)))
      if(length(closest)>1){closest=closest[1]}
      if(method=="Newton"){
        newtoninit = beta_mat[,closest]
        if(abs(newtoninit[length(newtoninit)])>0.95){newtoninit = init}
        if(abs(newtoninit[length(newtoninit)-1])<0.05){newtoninit = init}
        if(converged_mat[closest]==FALSE){newtoninit=init}
        fit = ss_pen_newton(newtoninit,lambda,eps,YS,XS,YO,XO,weightt)
        theta_temp = fit$theta
      }
      if(method=="Trust"){
        theta_temp <- trust_ss(init,lambda,eps,YS,XS,YO,XO,weightt)
      }
      bic_trial <- bic_fun(beta=threshold(theta_temp),YS,XS,YO,XO)
      if(bic_trial<bic_de[samp[3]]){
        lambda_de[samp[3],] <- offspring
        bic_de[samp[3]] <- bic_trial
        beta_mat[,samp[3]] <- theta_temp
        converged_mat[samp[3]] <- fit$converged
      }
    }
    deiter = deiter+1
  }
  selected = which(bic_de == min(bic_de))
  if(length(selected)>1){
    selected = selected[1]
  }
  results <- list(lambda=lambda_de[selected,],bic=bic_de[selected],theta=threshold(beta_mat[,selected]))
  return(results)
}

####################################################################
#DE1D
######################################################################

de_lambda1d <- function(lambda_max,YS,XS,YO,XO,weightt,init,popsize,niter,nparout,nparsel,method="Newton",eps = 0.00000001){
  print("Initialising")
  n <- length(YO)
  family <- seq(from=1,to=popsize)
  scale_factor <- 0.5
  lambda_de = seq(from=0,to=lambda_max,length.out=popsize+1)
  beta_mat <- matrix(NA,nrow=length(init),ncol=popsize+1)
  beta_mat[,1] = init
  bic_de <- rep(NA,popsize+1)
  converged_mat = rep(TRUE,popsize+1)
  bic_de[1]=bic_fun(beta=beta_mat[,1],YS,XS,YO,XO)
  for(m in 1:popsize){
    lambda=lambda_de[m+1]
    if(method=="Newton"){
      newtoninit = beta_mat[,m]
      if(abs(newtoninit[length(newtoninit)])>0.95){newtoninit = init}
      if(abs(newtoninit[length(newtoninit)-1])<0.05){newtoninit = init}
      if(converged_mat[m]==FALSE){newtoninit=init}
      fit = ss_pen_newton(newtoninit,lambda,eps,YS,XS,YO,XO,weightt)
      theta_temp = fit$theta
      converged_mat[m+1] = fit$converged
    }
    if(method=="Trust"){
      theta_temp <- trust_ss(init,lambda,eps,YS,XS,YO,XO,weightt)
    }
    bic_de[m+1] <- bic_fun(beta=threshold(theta_temp),YS,XS,YO,XO)
    beta_mat[,m+1] <- theta_temp
  }
  print("Initialised")
  deiter = 1
  convergedde = FALSE
  while(convergedde==FALSE && deiter<=niter){
    print(deiter)
    for(j in 1:popsize+1){
      samp <- sample(family[-j],size=3)
      trial = lambda_de[j]+(scale_factor*(lambda_de[samp[1]]-lambda_de[samp[2]]))
      prob=runif(2,0,1)
      offspring <- trial
      if(trial<0){trial=0}
      if(trial>lambda_max){trial=lambda_max}
      offspring <- trial
      lambda=offspring
      closest = which((abs(lambda_de-lambda))==min(abs(lambda_de-lambda)))
      if(length(closest)>1){closest=closest[1]}
      if(method=="Newton"){
        newtoninit = beta_mat[,closest]
        if(abs(newtoninit[length(newtoninit)])>0.95){newtoninit = init}
        if(abs(newtoninit[length(newtoninit)-1])<0.05){newtoninit = init}
        if(converged_mat[closest]==FALSE){newtoninit=init}
        fit = ss_pen_newton(newtoninit,lambda,eps,YS,XS,YO,XO,weightt)
        theta_temp = fit$theta
      }
      if(method=="Trust"){
        theta_temp <- trust_ss(init,lambda,eps,YS,XS,YO,XO,weightt)
      }
      bic_trial <- bic_fun(beta=threshold(theta_temp),YS,XS,YO,XO)
      if(bic_trial<bic_de[samp[3]]){
        lambda_de[samp[3]] <- offspring
        bic_de[samp[3]] <- bic_trial
        beta_mat[,samp[3]] <- theta_temp
        converged_mat[samp[3]] <- fit$converged
      }
    }
    deiter = deiter+1
  }
  selected = which(bic_de == min(bic_de))
  if(length(selected)>1){
    selected = selected[1]
  }
  results <- list(lambda=lambda_de[selected],bic=bic_de[selected],theta=threshold(beta_mat[,selected]))
  return(results)
}
######################################################################################
#Adjusted Grid Search
#####################################################################################
grid1 <- function(lambda_vec,init,YS,XS,YO,XO,weightt,nparout,nparsel,eps = 0.00000001,method="Trust"){
  n <- length(YO)
  models <- matrix(NA,nrow=length(init),ncol=length(lambda_vec)+1)
  models[,1] = init
  models_bic <- rep(NA,times=length(lambda_vec)+1)
  models_bic[1] <- bic_fun(init,YS,XS,YO,XO)
  for(i in 1:length(lambda_vec)){
    lambda = lambda_vec[i]
    if(method=="Newton"){
      theta_temp = fit$theta
      fit = ss_pen_newton(init,lambda,eps,YS,XS,YO,XO,weightt)
    }
    if(method=="Trust"){
      theta <- trust_ss(init,lambda,eps,YS,XS,YO,XO,weightt)
    }
    bic <- bic_fun(threshold(theta),YS,XS,YO,XO)
    models[,i+1]<- theta
    models_bic[i+1] <- bic
  }
  selected = which(models_bic == min(models_bic))
  if(length(selected)>1){
    selected = selected[1]
  }
  result=list(coefficients=threshold(models[,selected]),bic=models_bic[selected],lambda = lambda_vec[selected])
  return(result)
}

