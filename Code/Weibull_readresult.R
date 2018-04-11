rm(list=ls()) # removing everything on work space.
## library
library(statmod)
library(mvQuad)
library(fastGHQuad)
library(Rcpp)
library(RcppParallel)
library(RcppArmadillo)
source("Code/real_splines.R")

#### nulliparous (bmi + age) ####
nullresult = read.csv("Data/Weibull_real/weibull_null_without_Q10.csv", sep=",", header = FALSE)
var.name = c('Seed','n','Runtime', 'over(SVD)', 'obese(SVD)', 'age(SVD)', 'age2(SVD)',
             'over(CS)', 'obese(CS)', 'age(CS)', 'age2(CS)', 'over(OVD)', 'obese(OVD)', 'age(OVD)', 'age2(OVD)',
             rep('xi1', 9), rep('xi2', 9), rep('xi3', 9), 
             'over(M)', 'obese(M)', 'age(M)', 'age2(M)', rep('phi', 9), 'sigma')
nullresult = na.omit(nullresult)
main.est = as.numeric(read.csv("Data/Weibull_real/noboot_weibull_null_Q10_7978997.csv",
                               sep=",", header = FALSE))
noboot.est =  formatC(c(main.est[c(4:15, 43:46)], mean(exp(main.est[56]))), digits = 4, format = "f")
avg.est =  formatC(c(colMeans(nullresult)[c(4:15, 43:46)], mean(exp(nullresult[,56]))), digits = 4, format = "f")
se.est = formatC(c(apply(nullresult, 2, sd)[c(4:15, 43:46)], (apply(nullresult, 2, sd)[56]) * mean(exp(nullresult[,56]))), digits = 4, format = "f")

nullresult.QciU = formatC(apply(nullresult, 2, findU), digits = 4, format = "f")[c(4:15, 43:46, 56)]
nullresult.QciL = formatC(apply(nullresult, 2, findL), digits = 4, format = "f")[c(4:15, 43:46, 56)]
nullresult.QciU[17] = formatC(exp(findU(nullresult[,56])), digits = 4, format = "f")
nullresult.QciL[17] = formatC(exp(findL(nullresult[,56])), digits = 4, format = "f")

tab = matrix(0, 4, 4)
colnames(tab) = c("Overweight", "Obese", "Age", "Age2")
rownames(tab) = c("SVD", "CS", "OVD", "M")
tab[1,] = c(paste(noboot.est[1], " [", nullresult.QciL[1], ", ", nullresult.QciU[1],"]", sep=""),
            paste(noboot.est[2], " [", nullresult.QciL[2], ", ", nullresult.QciU[2],"]", sep=""),
            paste(noboot.est[3], " [", nullresult.QciL[3], ", ", nullresult.QciU[3],"]", sep=""),
            paste(noboot.est[4], " [", nullresult.QciL[4], ", ", nullresult.QciU[4],"]", sep=""))
tab[2,] = c(paste(noboot.est[5], " [", nullresult.QciL[5], ", ", nullresult.QciU[5],"]", sep=""),
            paste(noboot.est[6], " [", nullresult.QciL[6], ", ", nullresult.QciU[6],"]", sep=""),
            paste(noboot.est[7], " [", nullresult.QciL[7], ", ", nullresult.QciU[7],"]", sep=""),
            paste(noboot.est[8], " [", nullresult.QciL[8], ", ", nullresult.QciU[8],"]", sep=""))
tab[3,] = c(paste(noboot.est[9], " [", nullresult.QciL[9], ", ", nullresult.QciU[9],"]", sep=""),
            paste(noboot.est[10], " [", nullresult.QciL[10], ", ", nullresult.QciU[10],"]", sep=""),
            paste(noboot.est[11], " [", nullresult.QciL[11], ", ", nullresult.QciU[11],"]", sep=""),
            paste(noboot.est[12], " [", nullresult.QciL[12], ", ", nullresult.QciU[12],"]", sep=""))
tab[4,] = c(paste(noboot.est[13], " [", nullresult.QciL[13], ", ", nullresult.QciU[13],"]", sep=""),
            paste(noboot.est[14], " [", nullresult.QciL[14], ", ", nullresult.QciU[14],"]", sep=""),
            paste(noboot.est[15], " [", nullresult.QciL[15], ", ", nullresult.QciU[15],"]", sep=""),
            paste(noboot.est[16], " [", nullresult.QciL[16], ", ", nullresult.QciU[16],"]", sep=""))
print(xtable(tab))

noboot.est[17]
c(nullresult.QciL[17],  nullresult.QciU[17])

### -loglike
#data = read.csv("Data/nullwithout.csv", header= TRUE, sep = ",")
n = nrow(data)
obs.T = data$obs.T; obs.Delta = data$obs.Delta
obs.U = data$obs.U
obs.bmi1 = data$obs.bmi1; obs.bmi2 = data$obs.bmi2; obs.bmi3 = data$obs.bmi3
obs.age = data$obs.age; obs.age2 = (obs.age)^2
obs.T = obs.T*60

sample.knot = quantile(obs.T, probs = c(0.1, 0.25, 0.50, 0.75, 0.9))
Bsplines = make.Bsplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)
Isplines = make.Isplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)
gauss = gauss.quad(30, kind="laguerre")
weights = gauss$weights; nodes = gauss$nodes
sourceCpp("Code/Rcpp/rcpp_null.cpp")
negloglike.weibull.null = loglikeli(main.est[4:length(main.est)])

#### nulliparous (bmi) #####
nullbmi = read.csv("Data/Weibull_real/weibull_null_bmi_Q10.csv", sep=",", header = FALSE)
var.name = c('Seed','n','Runtime', 'over(SVD)', 'obese(SVD)', 
             'over(CS)', 'obese(CS)', 'over(OVD)', 'obese(OVD)', 'over(M)', 'obese(M)',
             rep('xi1', 9), rep('xi2', 9), rep('xi3', 9), rep('phi', 9), 'sigma')
nullbmi = na.omit(nullbmi)
main.est = as.numeric(read.csv("Data/Weibull_real/noboot_weibull_null_bmi_Q10_7978997.csv",
                               sep=",", header = FALSE))
noboot.est =  formatC(c(main.est[c(4:11)], mean(exp(main.est[48]))), digits = 4, format = "f")
avg.est =  formatC(c(colMeans(nullbmi)[c(4:11)], mean(exp(nullbmi[,48]))), digits = 4, format = "f")
se.est = formatC(c(apply(nullbmi, 2, sd)[c(4:11)], (apply(nullbmi, 2, sd)[48]) * mean(exp(nullbmi[,48]))), digits = 4, format = "f")

nullbmi.QciU = formatC(apply(nullbmi, 2, findU), digits = 4, format = "f")[c(4:11, 48)]
nullbmi.QciL = formatC(apply(nullbmi, 2, findL), digits = 4, format = "f")[c(4:11, 48)]
nullbmi.QciU[9] = formatC(exp(findU(nullbmi[,48])), digits = 4, format = "f")
nullbmi.QciL[9] = formatC(exp(findL(nullbmi[,48])), digits = 4, format = "f")

tab = matrix(0, 4, 2)
colnames(tab) = c("Overweight", "Obese")
rownames(tab) = c("SVD", "CS", "OVD", "M")
tab[1,] = c(paste(noboot.est[1], " [", nullbmi.QciL[1], ", ", nullbmi.QciU[1],"]", sep=""),
            paste(noboot.est[2], " [", nullbmi.QciL[2], ", ", nullbmi.QciU[2],"]", sep=""))
tab[2,] = c(paste(noboot.est[3], " [", nullbmi.QciL[3], ", ", nullbmi.QciU[3],"]", sep=""),
            paste(noboot.est[4], " [", nullbmi.QciL[4], ", ", nullbmi.QciU[4],"]", sep=""))
tab[3,] = c(paste(noboot.est[5], " [", nullbmi.QciL[5], ", ", nullbmi.QciU[5],"]", sep=""),
            paste(noboot.est[6], " [", nullbmi.QciL[6], ", ", nullbmi.QciU[6],"]", sep=""))
tab[4,] = c(paste(noboot.est[7], " [", nullbmi.QciL[7], ", ", nullbmi.QciU[7],"]", sep=""),
            paste(noboot.est[8], " [", nullbmi.QciL[8], ", ", nullbmi.QciU[8],"]", sep=""))
print(xtable(tab))

noboot.est[9]
c(nullbmi.QciL[9],  nullbmi.QciU[9])

### -loglike
#data = read.csv("Data/nullwithout.csv", header= TRUE, sep = ",")
n = nrow(data)
obs.T = data$obs.T; obs.Delta = data$obs.Delta
obs.U = data$obs.U
obs.bmi1 = data$obs.bmi1; obs.bmi2 = data$obs.bmi2; obs.bmi3 = data$obs.bmi3
obs.age = data$obs.age;obs.age2 = (obs.age)^2
obs.T = obs.T*60
sample.knot = quantile(obs.T, probs = c(0.1, 0.25, 0.50, 0.75, 0.9))
Bsplines = make.Bsplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)
Isplines = make.Isplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)
gauss = gauss.quad(30, kind="laguerre")
weights = gauss$weights; nodes = gauss$nodes

sourceCpp("Code/Rcpp/rcpp_null_bmi.cpp")
negloglike.weibull.nullbmi = loglikeli(main.est[4:length(main.est)])
##### nulliparous (age) #####
nullage = read.csv("Data/Weibull_real/weibull_null_age_Q10.csv", sep=",", header = FALSE)
var.name = c('Seed','n','Runtime', 'age(SVD)', 'age2(SVD)', 
             'age(CS)', 'age2(CS)', 'age(OVD)', 'age2(OVD)', 'age(M)', 'age2(M)',
             rep('xi1', 9), rep('xi2', 9), rep('xi3', 9), rep('phi', 9), 'sigma')
nullage = na.omit(nullage)
main.est = as.numeric(read.csv("Data/Weibull_real/noboot_weibull_null_age_Q10_7978997.csv",
                               sep=",", header = FALSE))
noboot.est =  formatC(c(main.est[c(4:11)], mean(exp(main.est[48]))), digits = 4, format = "f")
avg.est =  formatC(c(colMeans(nullage)[c(4:11)], mean(exp(nullage[,48]))), digits = 4, format = "f")
se.est = formatC(c(apply(nullage, 2, sd)[c(4:11)], (apply(nullage, 2, sd)[48]) * mean(exp(nullage[,48]))), digits = 4, format = "f")

nullage.QciU = formatC(apply(nullage, 2, findU), digits = 4, format = "f")[c(4:11, 48)]
nullage.QciL = formatC(apply(nullage, 2, findL), digits = 4, format = "f")[c(4:11, 48)]
nullage.QciU[9] = formatC(exp(findU(nullage[,48])), digits = 4, format = "f")
nullage.QciL[9] = formatC(exp(findL(nullage[,48])), digits = 4, format = "f")

tab = matrix(0, 4, 2)
colnames(tab) = c("Age", "Age2")
rownames(tab) = c("SVD", "CS", "OVD", "M")
tab[1,] = c(paste(noboot.est[1], " [", nullage.QciL[1], ", ", nullage.QciU[1],"]", sep=""),
            paste(noboot.est[2], " [", nullage.QciL[2], ", ", nullage.QciU[2],"]", sep=""))
tab[2,] = c(paste(noboot.est[3], " [", nullage.QciL[3], ", ", nullage.QciU[3],"]", sep=""),
            paste(noboot.est[4], " [", nullage.QciL[4], ", ", nullage.QciU[4],"]", sep=""))
tab[3,] = c(paste(noboot.est[5], " [", nullage.QciL[5], ", ", nullage.QciU[5],"]", sep=""),
            paste(noboot.est[6], " [", nullage.QciL[6], ", ", nullage.QciU[6],"]", sep=""))
tab[4,] = c(paste(noboot.est[7], " [", nullage.QciL[7], ", ", nullage.QciU[7],"]", sep=""),
            paste(noboot.est[8], " [", nullage.QciL[8], ", ", nullage.QciU[8],"]", sep=""))
print(xtable(tab))

noboot.est[9]
c(nullage.QciL[9],  nullage.QciU[9])

### -loglike
# data = read.csv("Data/nullwithout.csv", header= TRUE, sep = ",")
n = nrow(data)
obs.T = data$obs.T
obs.Delta = data$obs.Delta
obs.U = data$obs.U
obs.bmi1 = data$obs.bmi1; obs.bmi2 = data$obs.bmi2; obs.bmi3 = data$obs.bmi3
obs.age = data$obs.age; obs.age2 = (obs.age)^2
obs.T = obs.T*60
sample.knot = quantile(obs.T, probs = c(0.1, 0.25, 0.50, 0.75, 0.9))
Bsplines = make.Bsplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)
Isplines = make.Isplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)

gauss = gauss.quad(30, kind="laguerre")
weights = gauss$weights; nodes = gauss$nodes

sourceCpp("Code/Rcpp/rcpp_null_age.cpp")
negloglike.weibull.nullage = loglikeli(main.est[4:length(main.est)])
############################################################
##### multiparous (bmi + age) #####
multiresult10 = read.csv("Data/Weibull_real/weibull_multi_beyond10_Q10.csv", sep=",", header = FALSE)
var.name = c('Seed','n','Runtime', 'over(SVD)', 'obese(SVD)', 'age(SVD)', 'age2(SVD)',
             'over(CS+OVD)', 'obese(CS+OVD)', 'age(CS+OVD)', 'age2(CS+OVD)', 'over(M)', 'obese(M)', 'age(M)', 'age2(M)',
             rep('xi1', 9), rep('xi2', 9), rep('phi', 9), 'sigma')
multiresult10 = na.omit(multiresult10)
main.est = as.numeric(read.csv("Data/Weibull_real/noboot_weibull_multi_beyond10_Q10_7978997.csv",
                               sep=",", header = FALSE))
noboot.est =  formatC(c(main.est[c(4:15)], mean(exp(main.est[43]))), digits = 4, format = "f")
avg.est =  formatC(c(colMeans(multiresult10)[c(4:15)], mean(exp(multiresult10[,43]))), digits = 4, format = "f")
se.est = formatC(c(apply(multiresult10, 2, sd)[c(4:15)], (apply(multiresult10, 2, sd)[43]) * mean(exp(multiresult10[,43]))), digits = 4, format = "f")

multiresult10.QciU = formatC(apply(multiresult10, 2, findU), digits = 4, format = "f")[c(4:15,43)]
multiresult10.QciL = formatC(apply(multiresult10, 2, findL), digits = 4, format = "f")[c(4:15,43)]
multiresult10.QciU[13] = formatC(exp(findU(multiresult10[,43])), digits = 4, format = "f")
multiresult10.QciL[13] = formatC(exp(findL(multiresult10[,43])), digits = 4, format = "f")

tab = matrix(0, 3, 4)
colnames(tab) = c("Overweight", "Obese", "Age", "Age2")
rownames(tab) = c("SVD", "CS+OVD","M")
tab[1,] = c(paste(noboot.est[1], " [", multiresult10.QciL[1], ", ", multiresult10.QciU[1],"]", sep=""),
            paste(noboot.est[2], " [", multiresult10.QciL[2], ", ", multiresult10.QciU[2],"]", sep=""),
            paste(noboot.est[3], " [", multiresult10.QciL[3], ", ", multiresult10.QciU[3],"]", sep=""),
            paste(noboot.est[4], " [", multiresult10.QciL[4], ", ", multiresult10.QciU[4],"]", sep=""))
tab[2,] = c(paste(noboot.est[5], " [", multiresult10.QciL[5], ", ", multiresult10.QciU[5],"]", sep=""),
            paste(noboot.est[6], " [", multiresult10.QciL[6], ", ", multiresult10.QciU[6],"]", sep=""),
            paste(noboot.est[7], " [", multiresult10.QciL[7], ", ", multiresult10.QciU[7],"]", sep=""),
            paste(noboot.est[8], " [", multiresult10.QciL[8], ", ", multiresult10.QciU[8],"]", sep=""))
tab[3,] = c(paste(noboot.est[9], " [", multiresult10.QciL[9], ", ", multiresult10.QciU[9],"]", sep=""),
            paste(noboot.est[10], " [", multiresult10.QciL[10], ", ", multiresult10.QciU[10],"]", sep=""),
            paste(noboot.est[11], " [", multiresult10.QciL[11], ", ", multiresult10.QciU[11],"]", sep=""),
            paste(noboot.est[12], " [", multiresult10.QciL[12], ", ", multiresult10.QciU[12],"]", sep=""))
print(xtable(tab))
noboot.est[13]
c(multiresult10.QciL[13],  multiresult10.QciU[13])

### - loglike
#data = read.csv("Data/multiwithout.csv", header= TRUE, sep = ",")
d.time = data$obs.T*60
ind = which(d.time > 10)
obs.T = d.time[ind]
obs.Delta = data$obs.Delta[ind]; obs.Delta = ifelse(obs.Delta == 3, 2, obs.Delta)
obs.U = data$obs.U[ind]
obs.bmi1 = data$obs.bmi1[ind]; obs.bmi2 = data$obs.bmi2[ind]; obs.bmi3 = data$obs.bmi3[ind]
obs.age = data$obs.age[ind]; obs.age2 = (obs.age)^2

n = length(obs.T)
sample.knot = quantile(obs.T, probs = c(0.1, 0.25, 0.50, 0.75, 0.9))
Bsplines = make.Bsplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)
Isplines = make.Isplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)
gauss = gauss.quad(30, kind="laguerre")
weights = gauss$weights; nodes = gauss$nodes

sourceCpp("Code/Rcpp/rcpp_multi.cpp")
negloglike.weibull.multi = loglikeli(main.est[4:length(main.est)])

##### multiparous (bmi) #####
multiresult10 = read.csv("Data/Weibull_real/weibull_multi_beyond10_bmi_Q10.csv", sep=",", header = FALSE)
var.name = c('Seed','n','Runtime', 'over(SVD)', 'obese(SVD)', 
             'over(CS+OVD)', 'obese(CS+OVD)', 'over(M)', 'obese(M)',
             rep('xi1', 9), rep('xi2', 9), rep('phi', 9), 'sigma')
multiresult10 = na.omit(multiresult10)
main.est = as.numeric(read.csv("Data/Weibull_real/noboot_weibull_multi_beyond10_bmi_Q10_7978997.csv",
                               sep=",", header = FALSE))
noboot.est =  formatC(c(main.est[c(4:9)], mean(exp(main.est[37]))), digits = 4, format = "f")
avg.est =  formatC(c(colMeans(multiresult10)[c(4:9)], mean(exp(multiresult10[,37]))), digits = 4, format = "f")
se.est = formatC(c(apply(multiresult10, 2, sd)[c(4:9)], (apply(multiresult10, 2, sd)[37]) * mean(exp(multiresult10[,37]))), digits = 4, format = "f")

multiresult10.QciU = formatC(apply(multiresult10, 2, findU), digits = 4, format = "f")[c(4:9,37)]
multiresult10.QciL = formatC(apply(multiresult10, 2, findL), digits = 4, format = "f")[c(4:9,37)]
multiresult10.QciU[7] = formatC(exp(findU(multiresult10[,37])), digits = 4, format = "f")
multiresult10.QciL[7] = formatC(exp(findL(multiresult10[,37])), digits = 4, format = "f")

tab = matrix(0, 3, 2)
colnames(tab) = c("Overweight", "Obese")
rownames(tab) = c("SVD", "CS+OVD","M")
tab[1,] = c(paste(noboot.est[1], " [", multiresult10.QciL[1], ", ", multiresult10.QciU[1],"]", sep=""),
            paste(noboot.est[2], " [", multiresult10.QciL[2], ", ", multiresult10.QciU[2],"]", sep=""))

tab[2,] = c(paste(noboot.est[3], " [", multiresult10.QciL[3], ", ", multiresult10.QciU[3],"]", sep=""),
            paste(noboot.est[4], " [", multiresult10.QciL[4], ", ", multiresult10.QciU[4],"]", sep=""))

tab[3,] = c(paste(noboot.est[5], " [", multiresult10.QciL[5], ", ", multiresult10.QciU[5],"]", sep=""),
            paste(noboot.est[6], " [", multiresult10.QciL[6], ", ", multiresult10.QciU[6],"]", sep=""))

print(xtable(tab))
noboot.est[7]
c(multiresult10.QciL[7],  multiresult10.QciU[7])

### - loglike
# data = read.csv("Data/multiwithout.csv", header= TRUE, sep = ",")
d.time = data$obs.T*60
ind = which(d.time > 10)
obs.T = d.time[ind]
obs.Delta = data$obs.Delta[ind]; obs.Delta = ifelse(obs.Delta == 3, 2, obs.Delta)
obs.U = data$obs.U[ind]
obs.bmi1 = data$obs.bmi1[ind]; obs.bmi2 = data$obs.bmi2[ind]; obs.bmi3 = data$obs.bmi3[ind]
obs.age = data$obs.age[ind]; obs.age2 = (obs.age)^2
n = length(obs.T)
sample.knot = quantile(obs.T, probs = c(0.1, 0.25, 0.50, 0.75, 0.9))
Bsplines = make.Bsplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)
Isplines = make.Isplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)
gauss = gauss.quad(30, kind="laguerre")
weights = gauss$weights; nodes = gauss$nodes

sourceCpp("Code/Rcpp/rcpp_multi_bmi.cpp")
negloglike.weibull.multibmi = loglikeli(main.est[4:length(main.est)])
##### multiparous (age) ######
multiresult10 = read.csv("Data/Weibull_real/weibull_multi_beyond10_age_Q10.csv", sep=",", header = FALSE)
var.name = c('Seed','n','Runtime', 'age(SVD)', 'age2(SVD)', 
             'age(CS+OVD)', 'age2(CS+OVD)', 'age(M)', 'age2(M)',
             rep('xi1', 9), rep('xi2', 9), rep('phi', 9), 'sigma')
multiresult10 = na.omit(multiresult10)
avg.est =  formatC(c(colMeans(multiresult10)[c(4:9)], mean(exp(multiresult10[,37]))), digits = 4, format = "f")
se.est = formatC(c(apply(multiresult10, 2, sd)[c(4:9)], (apply(multiresult10, 2, sd)[37]) * mean(exp(multiresult10[,37]))), digits = 4, format = "f")
main.est = as.numeric(read.csv("Data/Weibull_real/noboot_weibull_multi_beyond10_age_Q10_7978997.csv",
                               sep=",", header = FALSE))
noboot.est =  formatC(c(main.est[c(4:9)], mean(exp(main.est[37]))), digits = 4, format = "f")
multiresult10.QciU = formatC(apply(multiresult10, 2, findU), digits = 4, format = "f")[c(4:9,37)]
multiresult10.QciL = formatC(apply(multiresult10, 2, findL), digits = 4, format = "f")[c(4:9,37)]
multiresult10.QciU[7] = formatC(exp(findU(multiresult10[,37])), digits = 4, format = "f")
multiresult10.QciL[7] = formatC(exp(findL(multiresult10[,37])), digits = 4, format = "f")

tab = matrix(0, 3, 2)
colnames(tab) = c("Age", "Age2")
rownames(tab) = c("SVD", "CS+OVD","M")
tab[1,] = c(paste(noboot.est[1], " [", multiresult10.QciL[1], ", ", multiresult10.QciU[1],"]", sep=""),
            paste(noboot.est[2], " [", multiresult10.QciL[2], ", ", multiresult10.QciU[2],"]", sep=""))

tab[2,] = c(paste(noboot.est[3], " [", multiresult10.QciL[3], ", ", multiresult10.QciU[3],"]", sep=""),
            paste(noboot.est[4], " [", multiresult10.QciL[4], ", ", multiresult10.QciU[4],"]", sep=""))

tab[3,] = c(paste(noboot.est[5], " [", multiresult10.QciL[5], ", ", multiresult10.QciU[5],"]", sep=""),
            paste(noboot.est[6], " [", multiresult10.QciL[6], ", ", multiresult10.QciU[6],"]", sep=""))

print(xtable(tab))

noboot.est[7]
c(multiresult10.QciL[7],  multiresult10.QciU[7])
### - loglike
# data = read.csv("data/multiwithout.csv", header= TRUE, sep = ",")
d.time = data$obs.T*60
ind = which(d.time > 10)
obs.T = d.time[ind]
obs.Delta = data$obs.Delta[ind]; obs.Delta = ifelse(obs.Delta == 3, 2, obs.Delta)
obs.U = data$obs.U[ind]
obs.bmi1 = data$obs.bmi1[ind]; obs.bmi2 = data$obs.bmi2[ind]; obs.bmi3 = data$obs.bmi3[ind]
obs.age = data$obs.age[ind]; obs.age2 = (obs.age)^2

n = length(obs.T)
sample.knot = quantile(obs.T, probs = c(0.1, 0.25, 0.50, 0.75, 0.9))
Bsplines = make.Bsplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)
Isplines = make.Isplines(obs.T, start = -0.1, knots = 5,
                         knots_position = sample.knot, degree = 3)
gauss = gauss.quad(30, kind="laguerre")
weights = gauss$weights; nodes = gauss$nodes

sourceCpp("Code/Rcpp/rcpp_multi_age.cpp")
negloglike.weibull.multiage = loglikeli(main.est[4:length(main.est)])

## -loglikelihood table combining lognormal and weibull
tab = matrix(0, 6, 5)
colnames(tab) = c("Parity", "Covariates", "The number of parameters", 
                  "Log-normal", "Weibull")
tab[,1] = c(rep("Nulliparous", 3), rep("Multiparous", 3))
tab[,2] = c("BMI + Age", "BMI", "Age", "BMI + Age", "BMI", "Age")
tab[,3] = c(53, 45, 45, 40, 34, 34)
tab[,4] = formatC(c(negloglike.lognormal.null, negloglike.lognormal.nullbmi, negloglike.lognormal.nullage),2 ,format = "f")
tab[,5] = formatC(c(negloglike.weibull.null, negloglike.weibull.nullbmi, negloglike.weibull.nullage), 2, format = "f")
print(xtable(tab))