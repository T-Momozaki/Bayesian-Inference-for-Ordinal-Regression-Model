
# Numerical experiments: latent variable models for binary response ------------------

# usage package
library(bayesplot) # for plot of mcmc sampling
library(gridExtra) # for plot of mcmc sampling
library(lemon) # for plot of mcmc sampling

# generate the simulation data (normal, t-dist, and logistic error)
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
e_l = rlogis(n, 0,1) # logistic error term
z_n = X %*% t_beta + e_n # latent variable with normal error
z_t = X %*% t_beta + e_t # latent variable with t-dist error
z_l = X %*% t_beta + e_l # latent variable with logistic error
y_n = as.numeric(z_n>0) # response with normal error
y_t = as.numeric(z_t>0) # response with t-dist error
y_l = as.numeric(z_l>0) # response with logistic error

# read function of implementing the Gibbs sampling for binary response
source("https://raw.githubusercontent.com/T-Momozaki/Bayesian-Inference-for-Ordinal-Regression-Model/lvm-binary/lvm-binary.R")

# implementing the Gibbs sampling using the uniform prior of beta for normal error
df_mc_n_n_u = gs_lvm_binary(y_n,X,model="probit",prior="unif",mc=10000,bn=2000,thin=4)
df_mc_n_t8_u = gs_lvm_binary(y_n,X,model="t-link",nu=8,prior="unif",mc=10000,bn=2000,thin=4)
df_mc_n_l_u = gs_lvm_binary(y_n,X,model="logit",prior="unif",mc=10000,bn=2000,thin=4)
# for t-dist. error with 4 degrees of freedom
df_mc_t_n_u = gs_lvm_binary(y_t,X,model="probit",prior="unif",mc=10000,bn=2000,thin=4)
df_mc_t_t8_u = gs_lvm_binary(y_t,X,model="t-link",nu=8,prior="unif",mc=10000,bn=2000,thin=4)
df_mc_t_l_u = gs_lvm_binary(y_t,X,model="logit",prior="unif",mc=10000,bn=2000,thin=4)
# for logistic error
df_mc_l_n_u = gs_lvm_binary(y_l,X,model="probit",prior="unif",mc=10000,bn=2000,thin=4)
df_mc_l_t8_u = gs_lvm_binary(y_l,X,model="t-link",nu=8,prior="unif",mc=10000,bn=2000,thin=4)
df_mc_l_l_u = gs_lvm_binary(y_l,X,model="logit",prior="unif",mc=10000,bn=2000,thin=4)

# implementing the Gibbs sampling using the normal prior of beta for normal error
df_mc_n_n_n = gs_lvm_binary(y_n,X,model="probit",prior="normal",mc=10000,bn=2000,thin=4)
df_mc_n_t8_n = gs_lvm_binary(y_n,X,model="t-link",nu=8,prior="normal",mc=10000,bn=2000,thin=4)
df_mc_n_l_n = gs_lvm_binary(y_n,X,model="logit",prior="normal",mc=10000,bn=2000,thin=4)
# for t-dist. error with 4 degrees of freedom
df_mc_t_n_n = gs_lvm_binary(y_t,X,model="probit",prior="normal",mc=10000,bn=2000,thin=4)
df_mc_t_t8_n = gs_lvm_binary(y_t,X,model="t-link",nu=8,prior="normal",mc=10000,bn=2000,thin=4)
df_mc_t_l_n = gs_lvm_binary(y_t,X,model="logit",prior="normal",mc=10000,bn=2000,thin=4)
# for logistic error
df_mc_l_n_n = gs_lvm_binary(y_l,X,model="probit",prior="normal",mc=10000,bn=2000,thin=4)
df_mc_l_t8_n = gs_lvm_binary(y_l,X,model="t-link",nu=8,prior="normal",mc=10000,bn=2000,thin=4)
df_mc_l_l_n = gs_lvm_binary(y_l,X,model="logit",prior="normal",mc=10000,bn=2000,thin=4)

# values of median and credible intervals
# apply(df_mc_n_n_u, 2, quantile, probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))
# apply(df_mc_n_t8_u, 2, quantile, probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))
# apply(df_mc_n_l_u, 2, quantile, probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))
# 
# apply(df_mc_t_n_u, 2, quantile, probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))
# apply(df_mc_t_t8_u, 2, quantile, probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))
# apply(df_mc_t_l_u, 2, quantile, probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))
# 
# apply(df_mc_l_n_u, 2, quantile, probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))
# apply(df_mc_l_t8_u, 2, quantile, probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))
# apply(df_mc_l_l_u, 2, quantile, probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))
# 
# apply(df_mc_n_n_n, 2, quantile, probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))
# apply(df_mc_n_t8_n, 2, quantile, probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))
# apply(df_mc_n_l_n, 2, quantile, probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))
# 
# apply(df_mc_t_n_n, 2, quantile, probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))
# apply(df_mc_t_t8_n, 2, quantile, probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))
# apply(df_mc_t_l_n, 2, quantile, probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))
# 
# apply(df_mc_l_n_n, 2, quantile, probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))
# apply(df_mc_l_t8_n, 2, quantile, probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))
# apply(df_mc_l_l_n, 2, quantile, probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))

# c.f. frequency probit and logit model
# fit_n_n = glm(y_n ~ .,
#               family = binomial(link = "probit"),
#               data = as.data.frame(cbind(y_n,X[,-1])))
# fit_n_l = glm(y_n ~ .,
#               family = binomial(link = "logit"),
#               data = as.data.frame(cbind(y_n,X[,-1])))
# fit_t_n = glm(y_t ~ .,
#               family = binomial(link = "probit"),
#               data = as.data.frame(cbind(y_t,X[,-1])))
# fit_t_l = glm(y_t ~ .,
#               family = binomial(link = "logit"),
#               data = as.data.frame(cbind(y_t,X[,-1])))
# fit_l_n = glm(y_l ~ .,
#               family = binomial(link = "probit"),
#               data = as.data.frame(cbind(y_l,X[,-1])))
# fit_l_l = glm(y_l ~ .,
#               family = binomial(link = "logit"),
#               data = as.data.frame(cbind(y_l,X[,-1])))

# plot the credible intervals and true parameter
p1_u = mcmc_recover_intervals(df_mc_n_n_u,t_beta) + 
  scale_y_continuous(breaks=seq(0.2,1.2,by=0.2), limits=c(0.2,1.2)) +
  labs(subtitle = "Probit for normal error")
p2_u = mcmc_recover_intervals(df_mc_n_t8_u,t_beta) + 
  scale_y_continuous(breaks=seq(0.2,1.2,by=0.2), limits=c(0.2,1.2)) +
  theme(legend.position = "none") +
  labs(subtitle = "Mixrure for normal error")
p3_u = mcmc_recover_intervals(df_mc_n_l_u,t_beta) + 
  scale_y_continuous(breaks=seq(0.2,1.2,by=0.2), limits=c(0.2,1.2)) +
  theme(legend.position = "none") +
  labs(subtitle = "Logit for normal error")
p4_u = mcmc_recover_intervals(df_mc_t_n_u,t_beta) + 
  scale_y_continuous(breaks=seq(0.2,1.2,by=0.2), limits=c(0.2,1.2)) +
  theme(legend.position = "none") +
  labs(subtitle = "Probit for t-dist. error with df=4")
p5_u = mcmc_recover_intervals(df_mc_t_t8_u,t_beta) + 
  scale_y_continuous(breaks=seq(0.2,1.2,by=0.2), limits=c(0.2,1.2)) +
  theme(legend.position = "none") +
  labs(subtitle = "Mixture for t-dist. error with df=4")
p6_u = mcmc_recover_intervals(df_mc_t_l_u,t_beta) + 
  scale_y_continuous(breaks=seq(0.2,1.2,by=0.2), limits=c(0.2,1.2)) +
  theme(legend.position = "none") +
  labs(subtitle = "Logit for t-dist. error with df=4")
p7_u = mcmc_recover_intervals(df_mc_l_n_u,t_beta) + 
  scale_y_continuous(breaks=seq(0.2,1.2,by=0.2), limits=c(0.2,1.2)) +
  theme(legend.position = "none") +
  labs(subtitle = "Probit for logistic error")
p8_u = mcmc_recover_intervals(df_mc_l_t8_u,t_beta) + 
  scale_y_continuous(breaks=seq(0.2,1.2,by=0.2), limits=c(0.2,1.2)) +
  theme(legend.position = "none") +
  labs(subtitle = "Mixture for logistic error")
p9_u = mcmc_recover_intervals(df_mc_l_l_u,t_beta) + 
  scale_y_continuous(breaks=seq(0.2,1.2,by=0.2), limits=c(0.2,1.2)) +
  theme(legend.position = "none") +
  labs(subtitle = "Logit for logistic error")
legend = lemon::g_legend(p1_u)
gridExtra::grid.arrange(p1_u+theme(legend.position = "none"),
                        p2_u,p3_u,p4_u,p5_u,p6_u,p7_u,p8_u,p9_u,
                        legend,
                        layout_matrix=rbind(c(1,1,2,2,3,3,10),c(4,4,5,5,6,6,10),c(7,7,8,8,9,9,10)))

p1_n = mcmc_recover_intervals(df_mc_n_n_n,t_beta) + 
  scale_y_continuous(breaks=seq(0.2,1.2,by=0.2), limits=c(0.2,1.2)) +
  labs(subtitle = "Probit for normal error")
p2_n = mcmc_recover_intervals(df_mc_n_t8_n,t_beta) + 
  scale_y_continuous(breaks=seq(0.2,1.2,by=0.2), limits=c(0.2,1.2)) +
  theme(legend.position = "none") +
  labs(subtitle = "Mixrure for normal error")
p3_n = mcmc_recover_intervals(df_mc_n_l_n,t_beta) + 
  scale_y_continuous(breaks=seq(0.2,1.2,by=0.2), limits=c(0.2,1.2)) +
  theme(legend.position = "none") +
  labs(subtitle = "Logit for normal error")
p4_n = mcmc_recover_intervals(df_mc_t_n_n,t_beta) + 
  scale_y_continuous(breaks=seq(0.2,1.2,by=0.2), limits=c(0.2,1.2)) +
  theme(legend.position = "none") +
  labs(subtitle = "Probit for t-dist. error with df=4")
p5_n = mcmc_recover_intervals(df_mc_t_t8_n,t_beta) + 
  scale_y_continuous(breaks=seq(0.2,1.2,by=0.2), limits=c(0.2,1.2)) +
  theme(legend.position = "none") +
  labs(subtitle = "Mixture for t-dist. error with df=4")
p6_n = mcmc_recover_intervals(df_mc_t_l_n,t_beta) + 
  scale_y_continuous(breaks=seq(0.2,1.2,by=0.2), limits=c(0.2,1.2)) +
  theme(legend.position = "none") +
  labs(subtitle = "Logit for t-dist. error with df=4")
p7_n = mcmc_recover_intervals(df_mc_l_n_n,t_beta) + 
  scale_y_continuous(breaks=seq(0.2,1.2,by=0.2), limits=c(0.2,1.2)) +
  theme(legend.position = "none") +
  labs(subtitle = "Probit for logistic error")
p8_n = mcmc_recover_intervals(df_mc_l_t8_n,t_beta) + 
  scale_y_continuous(breaks=seq(0.2,1.2,by=0.2), limits=c(0.2,1.2)) +
  theme(legend.position = "none") +
  labs(subtitle = "Mixture for logistic error")
p9_n = mcmc_recover_intervals(df_mc_l_l_n,t_beta) + 
  scale_y_continuous(breaks=seq(0.2,1.2,by=0.2), limits=c(0.2,1.2)) +
  theme(legend.position = "none") +
  labs(subtitle = "Logit for logistic error")
legend = lemon::g_legend(p1_n)
gridExtra::grid.arrange(p1_n+theme(legend.position = "none"),
                        p2_n,p3_n,p4_n,p5_n,p6_n,p7_n,p8_n,p9_n,
                        legend,
                        layout_matrix=rbind(c(1,1,2,2,3,3,10),c(4,4,5,5,6,6,10),c(7,7,8,8,9,9,10)))

# # trace plot
# mcmc_trace(df_mc_n_n_u)
# mcmc_trace(df_mc_n_t8_u)
# mcmc_trace(df_mc_n_l_u)
# 
# mcmc_trace(df_mc_t_n_u)
# mcmc_trace(df_mc_t_t8_u)
# mcmc_trace(df_mc_t_n_u)
# 
# mcmc_trace(df_mc_l_n_u)
# mcmc_trace(df_mc_l_t8_u)
# mcmc_trace(df_mc_l_l_u)
# 
# mcmc_trace(df_mc_n_n_n)
# mcmc_trace(df_mc_n_t8_n)
# mcmc_trace(df_mc_n_l_n)
# 
# mcmc_trace(df_mc_t_n_n)
# mcmc_trace(df_mc_t_t8_n)
# mcmc_trace(df_mc_t_l_n)
# 
# mcmc_trace(df_mc_l_n_n)
# mcmc_trace(df_mc_l_t8_n)
# mcmc_trace(df_mc_l_l_n)
# 
# # plot the posterior density
# mcmc_dens(df_mc_n_l_u)
# mcmc_dens(df_mc_t_l_u)
# mcmc_dens(df_mc_l_l_u)
# 
# mcmc_dens(df_mc_n_l_n)
# mcmc_dens(df_mc_t_l_n)
# mcmc_dens(df_mc_l_l_n)
# 
# # plot the auto correlation
# mcmc_acf(df_mc_n_l_u)
# mcmc_acf(df_mc_t_l_u)
# mcmc_acf(df_mc_l_l_u)
# 
# mcmc_acf(df_mc_n_l_n)
# mcmc_acf(df_mc_t_l_n)
# mcmc_acf(df_mc_l_l_n)




