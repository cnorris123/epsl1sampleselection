#################################################################
#Data
###################################################################
rm(list=ls())

gendat <-function(seed,n,true_beta,sigma2){
  set.seed(seed)
  x <- matrix(data=NA,nrow=n,ncol=(length(true_beta)-1))
  for(i in 1:ncol(x)){
    x[,i] = rnorm(n,0,1)
  }
  y <- true_beta[1]+x%*%true_beta[-1]+sqrt(sigma2)*rnorm(n)
  dat<-data.frame(y, x)
  return(dat)
}


##########################################################################
# Write analytical gradient and hessian for linear regression
#############################################################

ols.lf <- function(theta, yo, xo) {
  beta <- theta[-1]
  sigma2 <- theta[1]
  if (sigma2 <= 0) return(NA)
  n <- nrow(xo)
  e <- yo - xo%*%beta                                  # t() = matrix transpose
  logl <- ((-n/2)*log(2*pi)) - ((n/2)*log(sigma2)) - ((t(e)%*%e)/(2*sigma2))
  return(-logl) # since optim() does minimisation by default.
}



ols.score <- function(theta, yo, xo) {
  beta <- theta[-1]
  sigma2 <- theta[1]
  if (sigma2 <= 0) return(rep(NA, length(theta)))
  
  n <- nrow(xo)
  e <- yo - xo %*% beta
  
  # Compute the score vector (negative gradient)
  dL_dbeta <- -t(xo) %*% e / sigma2
  dL_dsigma2 <- ((n / (2 * sigma2)) - (t(e) %*% e) / (2 * sigma2^2))
  
  # Combine the partial derivatives into the score vector
  score <- c(dL_dsigma2, dL_dbeta)
  
  return(score)
}

ols.hessian <- function(theta, yo, xo) {
  beta <- theta[-1]
  sigma2 <- theta[1]
  if (sigma2 <= 0) return(matrix(NA, nrow = length(theta), ncol = length(theta)))
  
  n <- nrow(xo)
  e <- yo - xo %*% beta
  
  # Compute the Hessian matrix
  X_scaled <- xo / sqrt(sigma2)  # Scaled X matrix
  hessian_beta <- t(X_scaled) %*% X_scaled
  hessian_sigma2 <- -((n / (2 * sigma2^2)) + (t(e) %*% e) / (sigma2^3))
  hessian_sigma2beta <- (t(xo) %*% e)/(sigma2^2)
  
  # Combine the Hessian submatrices
  hessian <- matrix(0, nrow = length(theta), ncol = length(theta))
  hessian[1, 1] <- hessian_sigma2
  hessian[1,2:length(theta)] <- hessian_sigma2beta
  hessian[2:length(theta),1] <- hessian_sigma2beta
  hessian[2:length(theta), 2:length(theta)] <- hessian_beta
  
  return(hessian)
}

require(numDeriv)
aa <- lm(y~., data=dat)
sigma2 <- (summary(aa)$sigma)^2
init <- c(sigma2,aa$coefficients)
theta <- init
grad(ols.lf,init,yo=yo,xo=xo)
ols.score(theta=init,xo=xo,yo=yo)

hessian(ols.lf,init,yo=yo,xo=xo)
ols.hessian(theta=init,xo=xo,yo=yo)

#######################################################################
# Penalty approximation functions
#######################################################################

penalty <- function(theta,lambda,eps,weightt){
  beta <- theta[-1]
  return(lambda*sum(weightt[-1]*sqrt(beta*beta+eps)))
}

pen_deriv <- function(theta,lambda,eps,weightt){
  beta <- theta[-1]
  return(weightt*c(0,((lambda*beta)/sqrt((t(beta)*beta)+eps))))
}

pen_hess <- function(theta,lambda,eps,weightt){
  beta<-theta[-1]
  hessianvec <- c(0,(lambda*eps)/((t(beta)*beta+eps)^(3/2)))
  return(diag(hessianvec*weightt))
}


#####################################################################
#Line Search
############################################################################


linear_pen_newton <- function(init,lambda,eps,xo,yo,weightt){
  beta_old <- init
  tol <- 1e-5
  max_iter <- 1000
  converged <- FALSE
  iter <- 0
  steplength <- 1
  like<-ols.lf(theta=beta_old,xo=xo,yo=yo) + penalty(beta_old,lambda,eps,weightt)
  while(!converged && iter < max_iter){
    hesspen <- ols.hessian(beta_old,yo=yo,xo=xo)+pen_hess(beta_old,lambda=lambda,eps=eps,weightt=weightt)
    gradpen <- ols.score(beta_old,yo=yo,xo=xo)+pen_deriv(beta_old,lambda=lambda,eps=eps,weightt=weightt)
    descent <- as.vector((solve(hesspen)%*%gradpen))
    steplength = 1
    trial <- beta_old - steplength*descent
    while(trial[1]<=0){
      steplength = steplength*0.2
      trial <- beta_old - steplength*descent
    }
    triallik <-ols.lf(theta=trial,xo=xo,yo=yo) + penalty(trial,lambda,eps,weightt)
    ntrials=0
    c1=10^(-4)
    suffdecrease = c1*t(gradpen)%*%descent
    stepfound = FALSE
    while(stepfound == FALSE && ntrials<10000){
      trial <- beta_old - steplength*descent
      triallik <-ols.lf(theta=trial,xo=xo,yo=yo) + penalty(trial,lambda,eps,weightt)
      if(triallik> like+suffdecrease*steplength){
        steplength=steplength*0.2
      }else{
        stepfound=TRUE
      }
      ntrials = ntrials + 1
    }
    like <- triallik
    beta_new <- beta_old - steplength*descent
    dif <- abs(beta_new[-1]-beta_old[-1])/steplength
    if(max(dif) < tol){
      converged <- TRUE 
    }
    beta_old <- beta_new
    steplength = min(steplength*2,1)
    iter <- iter + 1
  }
  threshold = (eps*10)^0.5
  for(i in 1:length(beta_new)){
    if(abs(beta_old[i])<threshold){beta_old[i]=0}
  }
  return(beta_old)
}


lambda <- 100
eps <- 0.000000001
beta_test <- linear_pen_newton(init=init,lambda=lambda,eps=eps,xo=xo,yo=yo,weightt=weightt)
beta_test


############################################################
#BIC Function
###########################################################


bic_fun <- function(beta,xo,yo){
  n <- length(yo)
  df <- sum(beta!=0)-1
  nvar <- df+1
  yhat=xo%*%beta
  residuals = (yo- yhat)
  mse <- colMeans(residuals^2)
  bic <- n*log(mse)+nvar*log(n)
  return(bic)
}

############################################################
#Select lambda based on Newton-Raphson function using BIC
############################################################

lambda_grid_search <- function(lambda_vec,xo,yo,weightt){
  lm <- lm(y~.,data=dat)
  sigma2 <- (summary(lm)$sigma)^2
  init <- c(sigma2,lm$coefficients)
  eps <- 0.000000001
  n <- length(yo)
  models <- matrix(NA,nrow=length(init),ncol=length(lambda_vec))
  models_bic <- rep(NA,times=length(lambda_vec))
  for(i in 1:length(lambda_vec)){
    lambda <- lambda_vec[i]
    theta <- linear_pen_newton(init=init,lambda=lambda,eps=eps,xo=xo,yo=yo,weightt=weightt)
    beta <- theta[-1]
    bic <- bic_fun(beta,xo,yo)
    models[,i]<- theta
    models_bic[i] <- bic
  }
  models
  selected = which(models_bic == min(models_bic))
  if(length(selected)>1){
    selected = selected[1]
  }
  result=list(coefficients=models[,selected],bic=models_bic[selected],lambda = lambda_vec[selected])
  return(result)
}

###################################################################
#Differential Evolution
##################################################################


de_lambda <- function(lambda_max,xo,yo,weightt,popsize,niter){
  lm <- lm(y~.,data=dat)
  sigma2 <- (summary(lm)$sigma)^2
  init <- c(sigma2,lm$coefficients) #selecting inital theta
  eps <- 0.000000001
  n <- length(yo)
  family <- seq(from=1,to=popsize)
  scale_factor <- 0.5
  lambda_de <- runif(n=popsize,min=0,max=lambda_max) #initialization
  bic_de <- rep(NA,popsize)
  beta_mat <- matrix(NA,nrow=length(init),ncol=popsize)
  for(m in 1:length(bic_de)){
    theta_temp <- linear_pen_newton(init=init,lambda=lambda_de[m],eps=eps,xo=xo,yo=yo,weightt=weightt)
    beta_temp <- theta_temp[-1]
    bic_de[m] <- bic_fun(beta_temp,xo,yo)
    beta_mat[,m] <- theta_temp
  }
  iterde = 0
  while(iterde < niter){
    for(j in 1:popsize){
      samp <- sample(family[-j],size=2)
      trial = lambda_de[j]+(scale_factor*(lambda_de[samp[1]]-lambda_de[samp[2]]))
      if(trial<0){trial=0}
      if(trial>lambda_max){trial=lambda_max}
      theta_trial <- linear_pen_newton(init=init,lambda=trial,eps=eps,xo=xo,yo=yo,weightt=weightt)
      bic_trial <- bic_fun(theta_trial[-1],xo,yo)
      if(bic_trial<bic_de[j]){
        lambda_de[j] <- trial
        bic_de[j] <- bic_trial
        beta_mat[,j] <- theta_trial
      }
    }
    iterde = iterde + 1
  }
  selected = which(bic_de == min(bic_de))
  if(length(selected)>1){
    selected = selected[1]
  }
  results <- list(lambda=lambda_de[selected],bic=bic_de[selected],theta=beta_mat[,selected])
  return(results)
}




######################################################################################
#DEOptim Package
#######################################################################################
require("DEoptim")
?DEoptim

defun <- function(lambda,xo,yo,weightt){
  lm <- lm(y~.,data=dat)
  sigma2 <- (summary(lm)$sigma)^2
  init <- c(sigma2,lm$coefficients) #selecting inital theta
  eps <- 0.00000001
  theta <- linear_pen_newton(init=init,lambda=lambda,eps=eps,xo=xo,yo=yo,weightt=weightt)
  bic <- bic_fun(theta[-1],xo,yo)
  return(bic)
}

