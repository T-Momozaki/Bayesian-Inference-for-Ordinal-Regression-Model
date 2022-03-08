
# Probit model with unknown scale parameter -------------------------------

# usage package
library(MASS) # for sampling from multivariate normal
library(truncnorm) # for sampling from truncated normal
library(invgamma) # for sampling from inverse gamma
library(ggplot2) # for plot of mcmc sampling
library(ggfortify) # for plot of mcmc sampling

# generate the simulation data
set.seed(1234)
n = 2000 # sample size
p = 2 # dimention of predictor
X = cbind(X0 = rep(1, n), 
          X1 = rnorm(n, 0,1)) # predictor
t_beta = c(0.2, 0.6) # coefficient parameter
t_sig = 1 # scale parameter
e = rnorm(n, 0,t_sig) # error term
z = X %*% t_beta + e # latent variable
y = as.numeric(z>0) # response


# Gibbs sampling ----------------------------------------------------------

mc = 2500 # the number of chain
bn = 500 # the number of burn-out

# posterior sample
post_beta = matrix(NA, mc, p) # for beta
post_sig = numeric(mc) # for scale parameter

# initial value
beta = c(0,0)
sig = 0.5

# prior parameter
A = c(0,0); B = diag(c(100,100)) # for beta
a = 3; b = 0.5 # for scale parameter

# implementation of Gibbs sampling
for (i in 1:mc) {
  # sampling z
  z = ifelse(y==1,
             rtruncnorm(n, a = 0, mean = X%*%beta, sd = sig),
             rtruncnorm(n, b = 0, mean = X%*%beta, sd = sig))
  
  # sampling scale parameter
  res = t(z-X%*%beta)%*%(z-X%*%beta)
  sig2 = rinvgamma(1, 
                   shape={(n/2)+a}, 
                   rate={b+(res/2)})
  sig = sqrt(sig2)
  post_sig[i] = sig
  
  # sampling beta
  beta = mvrnorm(1,
                 mu = solve(t(X)%*%X/sig2+solve(B)) %*% {t(X)%*%z/sig2+solve(B)%*%A},
                   Sigma = solve(t(X)%*%X/sig2+solve(B)))
  post_beta[i,] = beta
}

df_mc = data.frame(
  post_beta[-(1:bn),2],
  post_sig[-(1:bn)],
  post_beta[-(1:bn),]/post_sig[-(1:bn)]
)
colnames(df_mc) = c("beta1", "sigma", "beta0/sigma", "beta1/sigma")

ts_mc = as.ts(df_mc)

p1 = autoplot(ts_mc[,1], main = "beta1")
p2 = autoplot(ts_mc[,2], main = "sigma")
p3 = autoplot(ts_mc[,3], main = "beta0/sigma")
p4 = autoplot(ts_mc[,4], main = "beta1/sigma")
gridExtra::grid.arrange(p1,p2,p3,p4)

