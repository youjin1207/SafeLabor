rm(list=ls()) # removing everything on work space.
## library
library(statmod)
library(mvQuad)
library(fastGHQuad)
library(Rcpp)
library(RcppParallel)
library(RcppArmadillo)

sourceCpp("Code/Rcpp/rcpp_fun.cpp")
source("Code/sim_splines.R")

## always look over the time we need
ptm<-proc.time()

#DUMZ=1000655 # seed is the number

## name of the seed / fix the number of seed. # we pick up one of the seed
set.seed(DUMZ)

# ------------------------------- setting parameters
n = 500 # n = 1000

sigma = 0.5 # sigma = 1.0
betas = c(1, -2); alpha = 1.5

# baselines of competing-risks
phis = c(0.01, 0.005)
tau = 3/2

# baselines of morbidities
psi = 0.01
omega = 2.5

# ------------------------------- generate data
eta = rnorm(n, 0, sigma)
X = rnorm(n, 0, 1)
W1 = runif(n, 0, 1)

# generate time to event from competing-risks   
d.time = (-log(W1) / (exp(-eta)*exp(betas[1]*X)*(phis[1])^(tau) + exp(-eta)*exp(betas[2]*X)*(phis[2])^(tau)))^(1 / tau) 

#generate competing risks data
csh1 = c(); csh2 = c(); p <- c(); Delta = c()
for(i in 1:n){
  csh1[i] = (d.time[i])^(tau - 1)*exp(-eta[i])*exp(betas[1]*X[i])*(phis[1]^(tau))*tau
  csh2[i] = (d.time[i])^(tau - 1)*exp(-eta[i])*exp(betas[2]*X[i])*(phis[2]^(tau))*tau
  p[i] = csh2[i] / (csh1[i] + csh2[i]) # probability of having C-section (delta = 2)
  Delta[i] = rbinom(1, 1, p[i]) 
  # probability of having C-section at t  = (cause-specific hazards of C-section at t) / (CSH of SVD at t + CSH of C-section at t)
  Delta[i] = Delta[i] + 1
}

W2 = runif(n, 0, 1) 
N = (-log(W2) / (psi^(omega)*exp(alpha*X)*exp(-eta)))^(1/omega) # unobserved onset time of neonatal morbidity 
U = (N <= d.time) # observed current status data.

# make time to event ordere ########  
ind = order(d.time, decreasing = FALSE)
obs.T = d.time[ind]
obs.Delta = Delta[ind]
obs.U = U[ind]
obs.X = X[ind]

obs.lambda1 = (phis[1])^(tau)*tau*(obs.T)^(tau-1)
obs.lambda2 = (phis[2])^(tau)*tau*(obs.T)^(tau-1)
obs.Lambda1 = (phis[1]*obs.T)^(tau); obs.Lambda2 = (phis[2]*obs.T)^(tau)

obs.h = (psi^omega)*omega*(obs.T)^(omega - 1)
obs.H = (psi*obs.T)^(omega)

sample.knot = quantile(obs.T, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
Bsplines = make.Bsplines(obs.T, start = 0, knots = 5,
                         knots_position = sample.knot, degree = 3)
Isplines = make.Isplines(obs.T, start = 0, knots = 5,
                         knots_position = sample.knot, degree = 3)

# ------------------------------- implement optimizations
rectime<-proc.time()-ptm

## implement optimization
init = rep(0,31)
gauss = gaussHermiteData(10)
weights = gauss$w
nodes = gauss$x

optim.result = try(optim(par = init, loglikeli, hessian = FALSE, method = "L-BFGS-B")$par, silent = TRUE)

# if the optim function prints out an error, let us run the rest of the codes 'siltently'..
if(class(optim.result) == "try-error"){
  result = rep(NA,31)
}else{
  result = optim.result
}


if(DUMZ==1000655){
  write.table(t(c(DUMZ,n,rectime[3],result)),file=paste('est_normal05_n1000_Q10_',DUMZ,'.csv',sep=''),sep=',',
              col.names=c('Seed','n','Runtime', 'beta1(0.0)', 'beta2(-2.0)', 'alpha(1.5)', 
                          rep('xi1', 9), rep('xi2', 9), rep('phi', 9), 'sigma(0.5)' ),row.names=F)
}else{
  write.table(t(c(DUMZ,n,rectime[3],result)),file=paste('est_normal05_n1000_Q10_',DUMZ,'.csv',sep=''),sep=',',col.names=F,row.names=F)
}

# ------------------------------- bootstrap
mcbt = matrix(0, nrow = 100, ncol = 31)
orgi.T = obs.T
orgi.Delta = obs.Delta
orgi.U = obs.U
orgi.X = obs.X
orgi.Bsplines = Bsplines
orgi.Isplines = Isplines

nb = 100

for(m in 1:nb){
  
  indbt = sample(1:n, replace=TRUE)
  
  d.time = orgi.T[indbt]
  Delta = orgi.Delta[indbt]
  U = orgi.U[indbt]
  X = orgi.X[indbt]
  
  ind = order(d.time, decreasing = FALSE)
  obs.T = d.time[ind]
  obs.Delta = Delta[ind]
  obs.U = U[ind]
  obs.X = X[ind]
  
  Bsplines = orgi.Bsplines[indbt[ind],]
  Isplines = orgi.Isplines[indbt[ind],]
  
  init = rep(0, 31)
  
  
  optim.result = try(optim(par = init, loglikeli, hessian = FALSE, method = "L-BFGS-B")$par, silent = TRUE)
  
  # if the optim function prints out an error, let us run the rest of the codes 'siltently'..
  if(class(optim.result) == "try-error"){
    result = rep(NA,31)
  }else{
    result = optim.result
  }
  
  
  mcbt[m, ] = result
}


findL = function(x){
  quantile(x,probs=0.025, na.rm = TRUE)
}

findU = function(x){
  quantile(x,probs=0.975, na.rm = TRUE)
}

mcsd = apply(mcbt,2,sd, na.rm = TRUE)
mcL = apply(mcbt,2,findL)
mcU = apply(mcbt,2,findU)
mcavg = apply(mcbt,2,mean, na.rm = TRUE)

rectime<-proc.time()-ptm

if(DUMZ==1000655){
  write.table(t(c(DUMZ,n,rectime[3],mcsd)),file=paste('Qsd_normal05_n1000_Q10', DUMZ,'.csv',sep=''),sep=',',
              col.names=c('Seed','n','Runtime', 'beta1(1.0)', 'beta2(-2.0)', 'alpha(1.5)', 
                          rep('xi1', 9), rep('xi2', 9), rep('phi', 9), 'sigma(0.5)' ),row.names=F)
  write.table(t(c(DUMZ,n,rectime[3],mcL)), file=paste('QciL_normal05_n1000_Q10', DUMZ,'.csv',sep=''),sep=',',
              col.names=c('Seed','n','Runtime', 'beta1(1.0)', 'beta2(-2.0)', 'alpha(1.5)', 
                          rep('xi1', 9), rep('xi2', 9), rep('phi', 9), 'sigma(0.5)' ),row.names=F)
  write.table(t(c(DUMZ,n,rectime[3],mcU)),file=paste('QciU_normal05_n1000_Q10', DUMZ,'.csv',sep=''),sep=',',
              col.names=c('Seed','n','Runtime', 'beta1(1.0)', 'beta2(-2.0)', 'alpha(1.5)',
                          rep('xi1', 9), rep('xi2', 9), rep('phi', 9), 'sigma(0.5)' ),row.names=F)
  write.table(t(c(DUMZ,n,rectime[3],mcavg)),file=paste('Qavg_normal05_n1000_Q10', DUMZ,'.csv',sep=''),sep=',',
              col.names=c('Seed','n','Runtime', 'beta1(1.0)', 'beta2(-2.0)', 'alpha(1.5)',
                          rep('xi1', 9), rep('xi2', 9), rep('phi', 9), 'sigma(0.5)' ),row.names=F)
}else{
  write.table(t(c(DUMZ,n,rectime[3],mcsd)),file=paste('Qsd_normal05_n1000_Q10', DUMZ,'.csv',sep=''),sep=',',
              col.names = F ,row.names=F)
  write.table(t(c(DUMZ,n,rectime[3],mcL)), file=paste('QciL_normal05_n1000_Q10', DUMZ,'.csv',sep=''),sep=',',
              col.names = F ,row.names=F)
  write.table(t(c(DUMZ,n,rectime[3],mcU)),file=paste('QciU_normal05_n1000_Q10', DUMZ,'.csv',sep=''),sep=',',
              col.names = F ,row.names=F)
  write.table(t(c(DUMZ,n,rectime[3],mcavg)),file=paste('Qavg_normal05_n1000_Q10', DUMZ,'.csv',sep=''),sep=',',
              col.names = F ,row.names=F)
}