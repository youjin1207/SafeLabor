rm(list=ls()) # removing everything on work space.
## library
library(statmod)
library(mvQuad)
library(fastGHQuad)
library(Rcpp)
library(RcppParallel)
library(RcppArmadillo)
source("Code/real_splines.R")
#### nulliparous women (age + bmi) ####
ptm<-proc.time()
sourceCpp("Code/Rcpp/rcpp_null.cpp")
DUMZ = 7978997 # random seed
set.seed(DUMZ)
# data = read.csv("Data/nullwithout.csv", header= TRUE, sep = ",") # data in restrictive use
n = nrow(data)
obs.T = data$obs.T
obs.Delta = data$obs.Delta
obs.U = data$obs.U
obs.bmi1 = data$obs.bmi1
obs.bmi2 = data$obs.bmi2
obs.bmi3 = data$obs.bmi3
obs.age = data$obs.age
obs.age2 = (obs.age)^2

obs.T = obs.T*60

sample.knot = quantile(obs.T, probs = c(0.1, 0.25, 0.50, 0.75, 0.9))
Bsplines = make.Bsplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)
Isplines = make.Isplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)

rectime = proc.time()-ptm
init = rep(0,53)
gauss = gauss.quad(10, kind="laguerre")
weights = gauss$weights
nodes = gauss$nodes
optim.result = try(optim(par = init, loglikeli, hessian = FALSE, method = "L-BFGS-B")$par, silent = TRUE)
# if the optim function prints out an error, let us run the rest of the codes 'siltently'..
if(class(optim.result) == "try-error"){
  result = rep(NA,53)
}else{
  result = optim.result
}
rectime<-proc.time()-ptm
write.table(t(c(DUMZ,n,rectime[3],result)),file=paste('Data/Weibull_real/noboot_weibull_null_Q10_',DUMZ,'.csv',sep=''),sep=',',
              col.names=c('Seed','n','Runtime', 'over(SVD)', 'obese(SVD)', 'age(SVD)', 'age2(SVD)',
                          'over(CS)', 'obese(CS)', 'age(CS)', 'age2(CS)', 'over(OVD)', 'obese(OVD)', 'age(OVD)', 'age2(OVD)',
                          rep('xi1', 9), rep('xi2', 9), rep('xi3', 9), 
                          'over(M)', 'obese(M)', 'age(M)', 'age2(M)', rep('phi', 9), 'sigma'),row.names=F)

#### nulliparous women (bmi) ####
ptm<-proc.time()
sourceCpp("rcpp_null_bmi.cpp")
DUMZ = 7978997 # random seed
set.seed(DUMZ)
# data = read.csv("Data/nullwithout.csv", header= TRUE, sep = ",")
n = nrow(data)
obs.T = data$obs.T
obs.Delta = data$obs.Delta
obs.U = data$obs.U
obs.bmi1 = data$obs.bmi1
obs.bmi2 = data$obs.bmi2
obs.bmi3 = data$obs.bmi3
obs.age = data$obs.age
obs.age2 = (obs.age)^2

obs.T = obs.T*60

sample.knot = quantile(obs.T, probs = c(0.1, 0.25, 0.50, 0.75, 0.9))
Bsplines = make.Bsplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)
Isplines = make.Isplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)

rectime = proc.time()-ptm

## implement optimization
init = rep(0,45)
gauss = gauss.quad(10, kind="laguerre")
weights = gauss$weights
nodes = gauss$nodes

optim.result = try(optim(par = init, loglikeli, hessian = FALSE, method = "L-BFGS-B")$par, silent = TRUE)

# if the optim function prints out an error, let us run the rest of the codes 'siltently'..
if(class(optim.result) == "try-error"){
  result = rep(NA,45)
}else{
  result = optim.result
}

rectime<-proc.time()-ptm
write.table(t(c(DUMZ,n,rectime[3],result)),file=paste('Data/Weibull_real/noboot_weibull_null_bmi_Q10_',DUMZ,'.csv',sep=''),sep=',',
              col.names=c('Seed','n','Runtime', 'over(SVD)', 'obese(SVD)', 
                          'over(CS)', 'obese(CS)', 'over(OVD)', 'obese(OVD)', 'over(M)', 'obese(M)',
                          rep('xi1', 9), rep('xi2', 9), rep('xi3', 9), rep('phi', 9), 'sigma'),row.names=F)

#### nulliparous women (age) ####
sourceCpp("Code/Rcpp/rcpp_null_age.cpp")
DUMZ = 7978997 # random seed
set.seed(DUMZ)
ptm = proc.time()
data = read.csv("nullwithout.csv", header= TRUE, sep = ",")
n = nrow(data)
obs.T = data$obs.T
obs.Delta = data$obs.Delta
obs.U = data$obs.U
obs.bmi1 = data$obs.bmi1
obs.bmi2 = data$obs.bmi2
obs.bmi3 = data$obs.bmi3
obs.age = data$obs.age
obs.age2 = (obs.age)^2

obs.T = obs.T*60

sample.knot = quantile(obs.T, probs = c(0.1, 0.25, 0.50, 0.75, 0.9))
Bsplines = make.Bsplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)
Isplines = make.Isplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)

## implement optimization
init = rep(0,45)
gauss = gauss.quad(10, kind="laguerre")
weights = gauss$weights
nodes = gauss$nodes

optim.result = try(optim(par = init, loglikeli, hessian = FALSE, method = "L-BFGS-B")$par, silent = TRUE)

# if the optim function prints out an error, let us run the rest of the codes 'siltently'..
if(class(optim.result) == "try-error"){
  result = rep(NA,45)
}else{
  result = optim.result
}
rectime = proc.time()-ptm
write.table(t(c(DUMZ,n,rectime[3],result)),file=paste('Data/Weibull_real/noboot_weibull_null_age_Q10_',DUMZ,'.csv',sep=''),sep=',',
              col.names=c('Seed','n','Runtime', 'age(SVD)', 'age2(SVD)', 
                          'age(CS)', 'age2(CS)', 'age(OVD)', 'age2(OVD)', 'age(M)', 'age2(M)',
                          rep('xi1', 9), rep('xi2', 9), rep('xi3', 9), rep('phi', 9), 'sigma'),row.names=F)

#### multiparous women (age + bmi) ####
sourceCpp('Code/Rcpp/rcpp_multi.cpp')
DUMZ = 7978997 # random seed
set.seed(DUMZ)
# data = read.csv("Data/multiwithout.csv", header= TRUE, sep = ",") # data in restrictive use
d.time = data$obs.T*60
ind = which(d.time > 10)
obs.T = d.time[ind]
obs.Delta = data$obs.Delta[ind]
obs.Delta = ifelse(obs.Delta == 3, 2, obs.Delta)
obs.U = data$obs.U[ind]
obs.bmi1 = data$obs.bmi1[ind]
obs.bmi2 = data$obs.bmi2[ind]
obs.bmi3 = data$obs.bmi3[ind]
obs.age = data$obs.age[ind]
obs.age2 = (obs.age)^2

n = length(obs.T)
sample.knot = quantile(obs.T, probs = c(0.1, 0.25, 0.50, 0.75, 0.9))
Bsplines = make.Bsplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)
Isplines = make.Isplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)
## implement optimization
init = rep(0,40)
gauss = gauss.quad(10, kind="laguerre")
weights = gauss$weights
nodes = gauss$nodes

optim.result = try(optim(par = init, loglikeli, hessian = FALSE, method = "L-BFGS-B")$par, silent = TRUE)

# if the optim function prints out an error, let us run the rest of the codes 'siltently'..
if(class(optim.result) == "try-error"){
  result = rep(NA,40)
}else{
  result = optim.result
}
rectime = proc.time()-ptm
write.table(t(c(DUMZ,n,rectime[3],result)),file=paste('Data/Weibull_real/noboot_weibull_multi_beyond10_Q10_',DUMZ,'.csv',sep=''),sep=',',
              col.names=c('Seed','n','Runtime', 'over(SVD)', 'obese(SVD)', 'age(SVD)', 'age2(SVD)',
                          'over(CS+OVD)', 'obese(CS+OVD)', 'age(CS+OVD)', 'age2(CS+OVD)', 'over(M)', 'obese(M)', 'age(M)', 'age2(M)',
                          rep('xi1', 9), rep('xi2', 9), rep('phi', 9), 'sigma'),row.names=F)


#### multiparous women (bmi) ####
ptm = proc.time()
DUMZ = 7978997 # random seed
set.seed(DUMZ)
# data = read.csv("Data/multiwithout.csv", header= TRUE, sep = ",")
d.time = data$obs.T*60
ind = which(d.time > 10)
obs.T = d.time[ind]
obs.Delta = data$obs.Delta[ind]
obs.Delta = ifelse(obs.Delta == 3, 2, obs.Delta)
obs.U = data$obs.U[ind]
obs.bmi1 = data$obs.bmi1[ind]
obs.bmi2 = data$obs.bmi2[ind]
obs.bmi3 = data$obs.bmi3[ind]
obs.age = data$obs.age[ind]
obs.age2 = (obs.age)^2

n = length(obs.T)
sample.knot = quantile(obs.T, probs = c(0.1, 0.25, 0.50, 0.75, 0.9))
Bsplines = make.Bsplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)
Isplines = make.Isplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)
## implement optimization
init = rep(0,34)
gauss = gauss.quad(10, kind="laguerre")
weights = gauss$weights
nodes = gauss$nodes
optim.result = try(optim(par = init, loglikeli, hessian = FALSE, method = "L-BFGS-B")$par, silent = TRUE)
# if the optim function prints out an error, let us run the rest of the codes 'siltently'..
if(class(optim.result) == "try-error"){
  result = rep(NA,34)
}else{
  result = optim.result
}
rectime = proc.time()-ptm
write.table(t(c(DUMZ,n,rectime[3],result)),file=paste('Data/Weibull_real/noboot_weibull_multi_beyond10_bmi_Q10_',DUMZ,'.csv',sep=''),sep=',',
              col.names=c('Seed','n','Runtime', 'over(SVD)', 'obese(SVD)', 
                          'over(CS+OVD)', 'obese(CS+OVD)', 'over(M)', 'obese(M)',
                          rep('xi1', 9), rep('xi2', 9), rep('phi', 9), 'sigma'),row.names=F)



#### multiparous women (age) ####
sourceCpp('Code/Rcpp/rcpp_multi_age.cpp')
ptm = proc.time()
DUMZ = 7978997 # random seed
set.seed(DUMZ)
# data = read.csv("Data/multiwithout.csv", header= TRUE, sep = ",") # data in restrictive use
d.time = data$obs.T*60
ind = which(d.time > 10)
obs.T = d.time[ind]
obs.Delta = data$obs.Delta[ind]
obs.Delta = ifelse(obs.Delta == 3, 2, obs.Delta)
obs.U = data$obs.U[ind]
obs.bmi1 = data$obs.bmi1[ind]
obs.bmi2 = data$obs.bmi2[ind]
obs.bmi3 = data$obs.bmi3[ind]
obs.age = data$obs.age[ind]
obs.age2 = (obs.age)^2

n = length(obs.T)
sample.knot = quantile(obs.T, probs = c(0.1, 0.25, 0.50, 0.75, 0.9))
Bsplines = make.Bsplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)
Isplines = make.Isplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)
## implement optimization
init = rep(0,34)
gauss = gauss.quad(10, kind="laguerre")
weights = gauss$weights
nodes = gauss$nodes
optim.result = try(optim(par = init, loglikeli, hessian = FALSE, method = "L-BFGS-B")$par, silent = TRUE)
# if the optim function prints out an error, let us run the rest of the codes 'siltently'..
if(class(optim.result) == "try-error"){
  result = rep(NA,34)
}else{
  result = optim.result
}
rectime = proc.time()-ptm
write.table(t(c(DUMZ,n,rectime[3],result)),file=paste('Data/Weibull_real/noboot_weibull_multi_beyond10_age_Q10_',DUMZ,'.csv',sep=''),sep=',',
              col.names=c('Seed','n','Runtime', 'age(SVD)', 'age2(SVD)', 
                          'age(CS+OVD)', 'age2(CS+OVD)','age(M)', 'age2(M)', 
                          rep('xi1', 9), rep('xi2', 9), rep('phi', 9), 'sigma'),row.names=F)