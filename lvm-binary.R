
# function: latent variable models for binary response --------------------

# usage package
library(MASS) # for sampling from multivariate normal
library(invgamma) # for sampling from inverse gamma

# implementation of Gibbs sampling for latent variable models for binary response
gs_lvm_binary = function(y, # binary response data
                         X, # predictors
                         model="probit", # selection of models, other, "t-link" is available
                         prior="normal", # prior for regression parameter, beta, other, "unif" (uniform prior) is available
                         beta=NULL, # intial value of beta
                         A=NULL, # mean of normal beta prior
                         B=NULL, # variance of normal beta prior
                         u=NULL, # initial value of local scale
                         nu=NULL, # hyperparameter of inverse-gamma prior of u
                         mc=2500, # length of MCMC
                         bn=500, # burn-in
                         seed=1234
){
  set.seed(1234)
  
  # sample size
  n = length(y)
  
  # dimention of predictor
  p = dim(X)[2]
  
  # posterior sample
  post_beta = matrix(NA, mc, p) # for beta
  z = numeric(n) # for latent variable
  
  # initial value
  beta = if(is.null(beta)) beta=rep(0,p) else beta=beta # for beta
  u = if(is.null(u)) u=rep(1,n) else u=u
  
  # prior parameter
  A = if(is.null(A)) A=rep(0,p) else A=A # mean of normal prior beta
  B = if(is.null(B)) B=diag(rep(100,p)) else B=B # for variance of normal prior beta
  nu = if(is.null(nu)) nu=8 else nu=nu
  
  # when the model is probit
  if(model=="probit") {
    
    for (i in 1:mc) {
      # sampling z (sampling from Truncated-Normal using inverse transform)
      lb = ifelse(y==1, 0, -Inf)-X%*%beta
      ub = ifelse(y==1, Inf, 0)-X%*%beta
      z = qnorm( runif(n, pnorm(lb,0,1), pnorm(ub,0,1)), X%*%beta, 1 )
      
      # sampling beta
      if(prior=="unif"){
        beta = mvrnorm(1,
                       mu = solve(t(X)%*%X) %*% t(X)%*%z,
                       Sigma = solve(t(X)%*%X))
      } else if(prior=="normal") {
        beta = mvrnorm(1,
                       mu = solve(t(X)%*%X+solve(B)) %*% {t(X)%*%z+solve(B)%*%A},
                       Sigma = solve(t(X)%*%X+solve(B)))
      }
      post_beta[i,] = beta
    }
    
  }
  
  # when the model is t-link
  if(model=="t-link") {
    
    for (i in 1:mc) {
      # sampling z
      lb = (ifelse(y==1, 0, -Inf)-X%*%beta)/sqrt(u)
      ub = (ifelse(y==1, Inf, 0)-X%*%beta)/sqrt(u)
      z = qnorm( runif(n, pnorm(lb,0,1), pnorm(ub,0,1)), X%*%beta, sqrt(u) )
      
      # sampling beta
      if(prior=="unif"){
        beta = mvrnorm(1,
                       mu = solve(t(X)%*%diag(1/u)%*%X) %*% t(X)%*%z,
                       Sigma = solve(t(X)%*%diag(1/u)%*%X))
      } else if(prior=="normal") {
        beta = mvrnorm(1,
                       mu = solve(t(X)%*%diag(1/u)%*%X+solve(B)) %*% {t(X)%*%diag(1/u)%*%z+solve(B)%*%A},
                       Sigma = solve(t(X)%*%diag(1/u)%*%X+solve(B)))
      }
      post_beta[i,] = beta
      
      # sampling u
      u = rinvgamma(n,
                    shape = {nu+1}/2,
                    rate = {nu+(z-X%*%beta)^2}/2 )
      
    }
    
  }
  
  df_mc = data.frame(post_beta[-(1:bn),])
  colnames(df_mc) = paste0("beta", 0:(p-1))
  
  return(df_mc)
}
