cbpen <- function(beta,corrmat,lambda){
  penvec = rep(0,length(beta)-1)
  for(i in 1:(length(beta)-1)){
    for(j in (i+1):length(beta)){
      penvec[i] = penvec[i]+(beta[i]-beta[j])^2/(1-corrmat[i,j]) + (beta[i]+beta[j])^2/(1+corrmat[i,j]) 
    }
  }
  penalty = 0.5*lambda*sum(penvec)
  return(penalty)
}

gradcbpen <- function(beta,corrmat,lambda){
  gradcb = rep(NA,length(beta))
  for(i in 1:length(beta)){
    j = (1:length(beta))[-i]
    gradcb[i] = 2*lambda*sum(((1-(corrmat[i,-i])^2)^(-1))*(rep(beta[i],length(beta)-1)-corrmat[i,-i]*beta[-i]))
  }
  return(gradcb)
}



hesscbpen <- function(beta,corrmat,lambda){
  penhess = matrix(NA,nrow=length(beta),ncol=length(beta))
  for(i in 1:length(beta)){
    for(j in 1:length(beta)){
      if(i==j){
        penhess[i,j] = 2*lambda*sum((1-(corrmat[i,-i])^2)^-1)
      }else{
        penhess[i,j] = -2*lambda*(corrmat[i,j])/(1-corrmat[i,j]^2)
      }
    }
  }
  return(penhess)
}

ss_corpen_newton <- function(init,lambda1,lambda2,YS,XS,YO,XO,weightt,ibeta,igamma,eps=1e-8){
  corrxo = cor(XO[,-1])
  corrxs = cor(XS[,-1])
  beta_old <- init
  tol <- 1e-5
  max_iter <- 10000
  converged <- FALSE
  steplength=1
  iter <- 0
  like <- loglik(beta_old,YS,XS,YO,XO) + penalty(beta_old,lambda1,eps,weightt) + cbpen(beta_old[ibeta],corrxo,lambda2) + cbpen(beta_old[igamma],corrxs,lambda2)
  while(!converged && iter<=max_iter && steplength>1e-15){
    betagradvec = rep(0,length(init))
    gammagradvec = rep(0,length(init))
    betagradvec[ibeta]=gradcbpen(beta_old[ibeta],corrxo,lambda2)
    gammagradvec[igamma] = gradcbpen(beta_old[igamma],corrxs,lambda2)
    betahessmat = matrix(0,nrow=length(init),ncol=length(init))
    betahessmat[ibeta,ibeta] = hesscbpen(beta_old[ibeta],corrxo,lambda2)
    gammahessmat = matrix(0,nrow=length(init),ncol=length(init))
    gammahessmat[igamma,igamma] = hesscbpen(beta_old[igamma],corrxo,lambda2)
    hesspen <- hesslik(beta_old,YS=YS,XS=XS,YO=YO,XO=XO)+pen_hess(beta_old,lambda=lambda1,eps=eps,weightt=weightt)+ betahessmat+gammahessmat
    gradpen <- gradlik(beta_old,YS=YS,XS=XS,YO=YO,XO=XO)+pen_deriv(beta_old,lambda=lambda1,eps=eps,weightt=weightt)+  betagradvec+gammagradvec
    descent <- as.vector((solve(hesspen)%*%gradpen))
    steplength=1
    trial <- beta_old - steplength*descent
    while(abs(trial[length(trial)])>=0.98 || trial[length(trial)-1]<=0.02){
      steplength=steplength*0.9
      trial <- beta_old - steplength*descent
    }
    ntrials=1
    c1=10^(-8)
    suffdecrease = c1*t(gradpen)%*%descent
    triallik <-loglik(trial,YS,XS,YO,XO) + penalty(trial,lambda1,eps,weightt) + cbpen(beta_old[ibeta],corrxo,lambda2) + cbpen(beta_old[igamma],corrxs,lambda2)
    stepfound = FALSE
    ntrials = 0
    while(stepfound == FALSE && ntrials <25){
      if(triallik> like-suffdecrease*steplength){
        steplength=steplength*0.1
        trial <- beta_old - steplength*descent
        triallik <-loglik(trial,YS,XS,YO,XO) + penalty(trial,lambda1,eps,weightt) + cbpen(beta_old[ibeta],corrxo,lambda2) + cbpen(beta_old[igamma],corrxs,lambda2)
      }else{
        stepfound=TRUE
      }
      ntrials = ntrials + 1
    }
    like <- triallik
    beta_new <- beta_old - steplength*descent
    if(max(abs(descent)) < tol){
      converged <- TRUE 
    }
    beta_old <- beta_new
    iter <- iter + 1
  }
  results = list(theta=beta_old,converged=converged)
  return(results)
}

de_lambda2dcorr <- function(lambda_max,YS,XS,YO,XO,weightt,init,popsize,niter,ibeta,igamma,method="Newton",pcr=0.7,eps = 0.00000001,ndim=2,lambda_min=0.01){
  library("lhs")
  print("Initialising")
  eps <- 0.00000001
  n <- length(YO)
  family <- seq(from=1,to=popsize)
  scale_factor <- 0.5
  lambda_de = maximinLHS(n=popsize,k=ndim)#initialization
  lambda_de[,1] = (lambda_de[,1]*(lambda_max[1]-lambda_min))+lambda_min
  lambda_de[,2] = (lambda_de[,2]*(lambda_max[2]-lambda_min))+lambda_min
  lambda_de = lambda_de[order(rowSums(lambda_de)^2,decreasing=F),]
  lambda_de = rbind(c(0,0),lambda_de)
  beta_mat <- matrix(NA,nrow=length(init),ncol=popsize+1)
  beta_mat[,1] = init
  bic_de <- rep(NA,popsize+1)
  bic_de[1]=bic_fun(beta=beta_mat[,1],YS,XS,YO,XO)
  converged_mat = rep(TRUE,popsize+1)
  for(m in 1:popsize){
    print(m)
    lambda=lambda_de[m,]
    if(method=="Newton"){
      newtoninit = beta_mat[,m]
      if(abs(newtoninit[length(newtoninit)])>0.8){newtoninit = init}
      if(abs(newtoninit[length(newtoninit)-1])<0.2){newtoninit = init}
      if(converged_mat[m]==FALSE){newtoninit=init}
      fit = ss_corpen_newton(newtoninit,lambda1=lambda_de[m,1],lambda2=lambda_de[m,2],YS,XS,YO,XO,weightt,ibeta,igamma)
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
      prob=runif(ndim,0,1)
      offspring <- lambda_de[samp[3],]
      for(i in 1:ndim){
        if(trial[i]<lambda_min){trial[i]=lambda_min}
        if(trial[i]>lambda_max[i]){trial[i]=lambda_max[i]}
        if(prob[i]<0.7){offspring[i]=trial[i]}
      }
      lambda=offspring
      closest = which(rowMeans(sweep(lambda_de,2,offspring)^2)==min(rowMeans(sweep(lambda_de,2,offspring)^2)))
      if(length(closest)>1){closest=closest[1]}
      if(method=="Newton"){
        newtoninit = beta_mat[,closest]
        if(abs(newtoninit[length(newtoninit)])>0.8){newtoninit = init}
        if(abs(newtoninit[length(newtoninit)-1])<0.2){newtoninit = init}
        if(converged_mat[closest]==FALSE){newtoninit=init}
        fit = ss_corpen_newton(newtoninit,lambda1=lambda[1],lambda2=lambda[2],YS,XS,YO,XO,weightt,ibeta,igamma)
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
