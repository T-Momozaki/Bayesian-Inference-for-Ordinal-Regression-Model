
# Numerical experiments: latent variable models for binary response ------------------

# usage package
library(bayesplot) # for plot of mcmc sampling
library(gridExtra) # for plot of mcmc sampling
library(lemon) # for plot of mcmc sampling

# generate the simulation data (normal and t-dist error)
set.seed(1234)
n = 500 # sample size
p = 4 # dimention of predictor
X = cbind(X0 = rep(1, n), 
          X1 = rnorm(n, 0,1),
          X2 = rnorm(n, 0,1),
          X3 = rnorm(n, 0,1)) # predictor
t_beta = c(0.5, 0.3, 0.7, 0.8) # coefficient parameter
e_n = rnorm(n, 0,1) # normal error term
e_t = rt(n, df=4) # t-dist. error term
z_n = X %*% t_beta + e_n # latent variable with normal error
z_t = X %*% t_beta + e_t # latent variable with t-dist error
y_n = as.numeric(z_n>0) # response with normal error
y_t = as.numeric(z_t>0) # response with t-dist error

# read function of implementing the Gibbs sampling for binary response
source("https://raw.githubusercontent.com/T-Momozaki/Bayesian-Inference-for-Ordinal-Regression-Model/binary-probit/lvm-binary.R")

# implementing the Gibbs sampling
df_mc_n_n = gs_lvm_binary(y_n,X,model="probit",prior="unif",mc=5000,bn=2000)
df_mc_n_t = gs_lvm_binary(y_n,X,model="t-link",prior="unif",mc=5000,bn=2000)
df_mc_t_n = gs_lvm_binary(y=y_t,X,model="probit",prior="unif",mc=5000,bn=2000)
df_mc_t_t8 = gs_lvm_binary(y=y_t,X,model="t-link",prior="unif", nu=3,mc=5000,bn=2000)

# values of median and credible intervals
apply(df_mc_n_n, 2, quantile, probs=c(0.025,0.5,0.975))
apply(df_mc_n_t, 2, quantile, probs=c(0.025,0.5,0.975))
apply(df_mc_t_n, 2, quantile, probs=c(0.025,0.5,0.975))
apply(df_mc_t_t, 2, quantile, probs=c(0.025,0.5,0.975))

# c.f. frequency probit and logit model
fit_n_n = glm(y_n ~ .,
              family = binomial(link = "probit"),
              data = as.data.frame(cbind(y_n,X[,-1])))
fit_n_t = glm(y_n ~ .,
              family = binomial(link = "logit"),
              data = as.data.frame(cbind(y_n,X[,-1])))
fit_t_n = glm(y_t ~ .,
              family = binomial(link = "probit"),
              data = as.data.frame(cbind(y_t,X[,-1])))
fit_t_t = glm(y_t ~ .,
              family = binomial(link = "logit"),
              data = as.data.frame(cbind(y_t,X[,-1])))

# plot the credible intervals and true parameter
p1 = mcmc_recover_intervals(df_mc_n_n,t_beta) + 
  scale_y_continuous(breaks=seq(0.2,1.2,by=0.2), limits=c(0.2,1.2)) +
  labs(subtitle = "Probit for normal error")

p2 = mcmc_recover_intervals(df_mc_n_t,t_beta) + 
  scale_y_continuous(breaks=seq(0.2,1.2,by=0.2), limits=c(0.2,1.2)) +
  theme(legend.position = "none") +
  labs(subtitle = "Mixrure for normal error")

p3 = mcmc_recover_intervals(df_mc_t_n,t_beta) + 
  scale_y_continuous(breaks=seq(0.1,0.9,by=0.2), limits=c(0.1,0.9)) +
  theme(legend.position = "none") +
  labs(subtitle = "Probit for t-dist error")

p4 = mcmc_recover_intervals(df_mc_t_t8,t_beta) + 
  scale_y_continuous(breaks=seq(0.1,0.9,by=0.2), limits=c(0.1,0.9)) +
  theme(legend.position = "none") +
  labs(subtitle = "Mixture for t-dist error")

legend = lemon::g_legend(p1)
gridExtra::grid.arrange(p1+theme(legend.position = "none"),p2,p3,p4,
                        legend,
                        layout_matrix=rbind(c(1,1,2,2,5),c(3,3,4,4,5)))

# trace plot
mcmc_trace(df_mc_n_n)
mcmc_trace(df_mc_n_t)
mcmc_trace(df_mc_t_n)
mcmc_trace(df_mc_t_t8)

# plot the posterior density
mcmc_dens(df_mc_n_n)
mcmc_dens(df_mc_n_t)
mcmc_dens(df_mc_t_n)
mcmc_dens(df_mc_t_t8)

# plot the auto correlation
mcmc_acf(df_mc_n_n)
mcmc_acf(df_mc_n_t)
mcmc_acf(df_mc_t_n)
mcmc_acf(df_mc_t_t8)




