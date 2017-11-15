setwd("/Users/luzhang/Documents/github/NNGP_STAN")
rm(list = ls())
library(rstan)
library(shinystan)
library(spNNGP)       # Build neighbor index
source("NNmatrix.R")  # Build matrix including nearest neighbor information

#------------------------- generate simulation data ---------------------------#
rmvn <- function(N, mu = 0, V = matrix(1)){
  P <- length(mu)
  if(any(is.na(match(dim(V), P))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(N * P), ncol = P) %*% D + rep(mu, rep(N, P)))
}

set.seed(1234)
N <- 500
coords <- cbind(runif(N), runif(N))
X <- as.matrix(cbind(1, rnorm(N)))
B <- as.matrix(c(1, 5))
sigma.sq <- 2
tau.sq <- 0.1
phi <- 3 / 0.5

D <- as.matrix(dist(coords))
R <- exp(- phi * D)
w <- rmvn(1, rep(0, N), sigma.sq*R)
Y <- rnorm(N, X %*% B + w, sqrt(tau.sq))

#----------------------- Empirival Variogram (for fun) ------------------------#
library(geoR)     # empirical variogram
library(spBayes)  # function iDist 

par(mfrow = c(1, 1))
lm.obj <- lm(Y ~ X[, -1])
d.max <- max(iDist(coords))
v.resid <- variog(coords=coords, data=resid(lm.obj), 
                  uvec=(seq(0, 0.6 * d.max, length = 20))) 
vario.fit <- variofit(v.resid, cov.model="exponential")
plot(v.resid, main = 'Empirical Variogram', cex = 0.2, ylim = c(0, 3.5))
lines(vario.fit, col='red')

#-------------- Build neighbor index on sorted x by "spConjNNGP" --------------#
M = 6 # Number of Nearest Neighbors
sigma.sq.IG <- c(2, 2)
cov.model <- "exponential"
theta.alpha <- c(5, 0.5)
names(theta.alpha) <- c("phi", "alpha")

m.c <- spConjNNGP(Y ~ X[, -1], 
                  coords = coords, 
                  n.neighbors = M,
                  theta.alpha = theta.alpha, k.fold = 1,
                  n.omp.threads = 2, return.neighbors = T,
                  sigma.sq.IG = sigma.sq.IG, cov.model = cov.model)

NN.matrix <- NNMatrix(N, m.c$coords.ord, m.c$n.indx[-1])

#------------------------- Check Neighbors (For fun) --------------------------#

par(mfrow=c(1,1))
Check_Neighbors(m.c$coords.ord, NN = M, NN.matrix, ind = 200)

#-------------------------- Set parameters of priors --------------------------#
P = 1                  # number of regression coefficients
ss = 3 * sqrt(2)       # scale parameters in the normal prior of sigma 
st = 3 * sqrt(0.1)     # scale parameters in the normal prior of tau     
ap = 3/1; bp = 3/0.1   # upper and lower bound of phi 

#--------------------------- NNGP random effects ------------------------------#
options(mc.cores = parallel::detectCores())
data <- list(N = N, M = M, P = P, Y = m.c$y.ord, X = m.c$X.ord, 
             NN_ind = NN.matrix$NN_ind, NN_dist = NN.matrix$NN_dist, 
             NN_distM = NN.matrix$NN_distM, 
             uB = rep(0, P + 1), VB = diag(P + 1)*1000,
             ss = ss, st = st, ap = ap, bp = bp)

myinits <-list(list(beta = c(1, 5), sigma = 1, tau = 0.5, phi = 20, 
                    w_b1 = rep(0, N)), 
               list(beta = c(5, 5), sigma = 1.5, tau = 0.2, phi = 5, 
                    w_b1 = rep(0.1, N)), 
               list(beta = c(0, 0), sigma = 2.5, tau = 0.1, phi = 9 ,
                    w_b1 = rep(0, N)))

parameters <- c("beta", "sigmasq", "tausq", "phi", "w")
samples_w <- stan(
  file = "nngp_random.stan",
  data = data,
  init = myinits,
  pars = parameters,
  iter = 400, 
  chains = 3,
  thin = 1,
  seed = 123
)

print(samples_w)
stan_trace(samples_w)
stan_trace(samples_w, inc_warmup = T)

sampler_params_w <- get_sampler_params(samples_w, inc_warmup = FALSE)
mean_accept_stat_by_chain <- 
  sapply(sampler_params_w, function(x) mean(x[, "accept_stat__"]))
print(mean_accept_stat_by_chain)

max_treedepth_by_chain <- 
  sapply(sampler_params_w, function(x) max(x[, "treedepth__"]))
print(max_treedepth_by_chain)

launch_shinystan(samples_w)

#----------------- NNGP random effects with w centered at b1 ------------------#
options(mc.cores = parallel::detectCores())
data <- list(N = N, M = M, P = P, Y = m.c$y.ord, X = m.c$X.ord, 
             NN_ind = NN.matrix$NN_ind, NN_dist = NN.matrix$NN_dist,
             NN_distM = NN.matrix$NN_distM, 
             uB = rep(0, P + 1), VB = diag(P + 1)*1000,
             ss = ss, st = st, ap = ap, bp = bp)

myinits <-list(list(beta = c(1, 5), sigma = 1, tau = 0.5, phi = 20, 
                    w_b1 = rep(0, N)), 
               list(beta = c(5, 5), sigma = 1.5, tau = 0.2, phi = 5, 
                    w_b1 = rep(0.1, N)), 
               list(beta = c(0, 0), sigma = 2.5, tau = 0.1, phi = 9 ,
                    w_b1 = rep(0, N)))

parameters <- c("beta", "sigmasq", "tausq", "phi", "w_b1")
samples_wb1 <- stan(
  file = "nngp_random_b1.stan",
  data = data,
  init = myinits,
  pars = parameters,
  iter = 400, 
  chains = 3,
  thin = 1,
  seed = 123
)

print(samples_wb1, pars = c("beta", "sigmasq", "tausq", "phi", "w_b1[1]", 
                            "w_b1[2]", "w_b1[3]", "w_b1[4]"))
stan_trace(samples_wb1, pars = c("beta", "sigmasq", "tausq", "phi", "w_b1[1]", 
                                 "w_b1[2]", "w_b1[3]", "w_b1[4]"))
stan_trace(samples_wb1, inc_warmup = T)

sampler_params_wb1 <- get_sampler_params(samples_wb1, inc_warmup = FALSE)
mean_accept_stat_by_chain <- 
  sapply(sampler_params_wb1, function(x) mean(x[, "accept_stat__"]))
print(mean_accept_stat_by_chain)

max_treedepth_by_chain <- 
  sapply(sampler_params_wb1, function(x) max(x[, "treedepth__"]))
print(max_treedepth_by_chain)

launch_shinystan(samples_wb1)

#------------------------------ NNGP response ---------------------------------#

options(mc.cores = parallel::detectCores())
data <- list(N = N, M = M, P = P, Y = m.c$y.ord, X = m.c$X.ord, 
             NN_ind = NN.matrix$NN_ind, NN_dist = NN.matrix$NN_dist, 
             NN_distM = NN.matrix$NN_distM,  
             uB = rep(0, P + 1), VB = diag(P + 1)*1000,
             ss = ss, st = st, ap = ap, bp = bp)

myinits <-list(list(beta = c(1, 5), sigma = 1, tau = 0.4, phi = 20), 
               list(beta = c(5, 5), sigma = 1.5, tau = 0.2, phi = 5), 
               list(beta = c(0, 0), sigma = 2.5, tau = 0.1, phi = 9))

parameters <- c("beta", "sigmasq", "tausq", "phi")
samples <- stan(
  file = "nngp_response.stan",
  data = data,
  init = myinits,
  pars = parameters,
  iter = 400, 
  chains = 3,
  thin = 1,
  seed = 123
)

print(samples)
stan_trace(samples)
stan_trace(samples, inc_warmup = T)

sampler_params <- get_sampler_params(samples, inc_warmup = FALSE)
mean_accept_stat_by_chain <- 
  sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
print(mean_accept_stat_by_chain)

max_treedepth_by_chain <- 
  sapply(sampler_params, function(x) max(x[, "treedepth__"]))
print(max_treedepth_by_chain)

launch_shinystan(samples)

