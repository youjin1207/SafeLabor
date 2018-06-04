source("Code/real_splines.R")
library(statmod)
library(mvQuad)
library(fastGHQuad)
library(Rcpp)
library(RcppParallel)
library(RcppArmadillo)
findL = function(x){
  quantile(x,probs=0.025, na.rm = TRUE)
}
findU = function(x){
  quantile(x,probs=0.975, na.rm = TRUE)
}
## Nulliparous (log-normal) ##
gauss = gaussHermiteData(30)
weights = gauss$w
nodes = gauss$x
## (1.a) age + BMI
# sourceCpp('Code/Rcpp/rcpp_fun_real.cpp')
# data = read.csv("nullwithout.csv", header= TRUE, sep = ",")
# nullresult = read.csv("Data/lognormal_real/null_without_Q10.csv", sep=",", header = FALSE)
# boot.seed = t(read.table("Data/Seed/seed_lognormal_null.txt"))
## (1.b) BMI
# sourceCpp('Code/Rcpp/rcpp_fun_null_bmi.cpp')
# data = read.csv("nullwithout.csv", header= TRUE, sep = ",")
# nullresult = read.csv("Data/lognormal_real/null_bmi_Q10.csv", sep=",", header = FALSE)
# boot.seed = t(read.table("Data/Seed/seed_lognormal_null_bmi.txt"))
## (1.c) age
# sourceCpp('Code/Rcpp/rcpp_fun_null_age.cpp')
# data = read.csv("nullwithout.csv", header= TRUE, sep = ",")
# nullresult = read.csv("Data/lognormal_real/null_age_Q10.csv", sep=",", header = FALSE)
# boot.seed = t(read.table("Data/Seed/seed_lognormal_null_age.txt"))

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

## (1.a) age + BMI
#main.est = as.numeric(read.csv("Data/lognormal_real/noboot_null_Q10_8169573.csv",
#                               sep=",", header = FALSE))
## (1.b) BMI
#main.est = as.numeric(read.csv("data/lognormal_real/null/noboot_null_bmi_Q10_8169573.csv",
#                               sep=",", header = FALSE))
## (1.c) age 
#main.est = as.numeric(read.csv("Data/lognormal_real/noboot_null_age_Q10_8169573.csv",
#                               sep=",", header = FALSE))
main.est = main.est[-c(1:3)]
frailty.range = seq(0.1, 3, 0.1)
gradient.lognormal.null  = rep(0, length(frailty.range))
for(r in 1:length(frailty.range)){
  print(r)
  for(i in 1:length(obs.T)){
    deno = c()
    for(q in 1:length(weights)){
      deno[q] = joint(nodes[q], main.est, i-1)*weights[q]
    }
    gradient.lognormal.null[r] = gradient.lognormal.null[r] + conditional(-log(frailty.range[r]), main.est, i-1) / max(sum(deno),10^(-20))
  }
}
## (1.a) age + BMI
# write.table(gradient.lognormal.null/n, "Data/lognormal_real/gradient_lognormal_null.txt", sep="\t")
## (1.b) BMI
# write.table(gradient.lognormal.null/n, "Data/lognormal_real/gradient_lognormal_null_bmi.txt", sep="\t")
## (1.c) age 
# write.table(gradient.lognormal.null/n, "Data/lognormal_real/gradient_lognormal_null_age.txt", sep="\t")

orgi.T = obs.T
orgi.Delta = obs.Delta
orgi.U = obs.U
orgi.bmi1 = obs.bmi1
orgi.bmi2 = obs.bmi2
orgi.bmi3 = obs.bmi3
orgi.age = obs.age
orgi.age2 = obs.age2
orgi.Bsplines = Bsplines
orgi.Isplines = Isplines

rectime = c()
for(i in 1:length(boot.seed)){
  rectime[i] = proc.time()-ptm
  main.est = as.numeric(nullresult[which(nullresult[,1] == boot.seed[i]),4:ncol(nullresult)])
  set.seed(boot.seed[i])
  n = nrow(data)
  indbt = sample(1:n, replace=TRUE)
  
  d.time = orgi.T[indbt]
  Delta = orgi.Delta[indbt]
  U = orgi.U[indbt]
  bmi1 = orgi.bmi1[indbt]
  bmi2 = orgi.bmi2[indbt]
  bmi3 = orgi.bmi3[indbt]
  age = orgi.age[indbt]
  age2 = orgi.age2[indbt]
  
  ind = order(d.time, decreasing = FALSE)
  obs.T = d.time[ind]
  obs.Delta = Delta[ind]
  obs.U = U[ind]
  obs.bmi1 = bmi1[ind]
  obs.bmi2 = bmi2[ind]
  obs.bmi3 = bmi3[ind]
  obs.age = age[ind]
  obs.age2 = age2[ind]
  
  Bsplines = orgi.Bsplines[indbt[ind],]
  Isplines = orgi.Isplines[indbt[ind],]
  
  frailty.range = seq(0.1, 3, 0.1)
  
  gradient.lognormal.null  = rep(0, length(frailty.range))
  for(r in 1:length(frailty.range)){
    for(i in 1:length(obs.T)){
      deno = c()
      for(q in 1:length(weights)){
        deno[q] = joint(nodes[q], main.est, i-1)*weights[q]
      }
      gradient.lognormal.null[r] = gradient.lognormal.null[r] + conditional(-log(frailty.range[r]), main.est, i-1) / max(sum(deno),10^(-20))
    }
  }
  
  result[i] = gradient.lognormal.null / n
  rectime[i] = proc.time()-ptm
}

## (1.a) age + BMI
#write.table(cbind(boot.seed, rep(n, length(boot.seed)), rectime, result),
#            file=paste('Data/lognormal_grad_null_without_Q10', '.csv',sep=''),sep=',', row.names=F, col.names = F)
## (1.b) BMI
#write.table(cbind(boot.seed, rep(n, length(boot.seed)), rectime, result),
#            file=paste('Data/lognormal_grad_null_bmi_Q10', '.csv',sep=''),sep=',', row.names=F, col.names = F)
## (1.c) age 
#write.table(cbind(boot.seed, rep(n, length(boot.seed)), rectime, result),
#            file=paste('Data/lognormal_grad_null_age_Q10', '.csv',sep=''),sep=',', row.names=F, col.names = F)


## Multiparous (log-normal) ##
gauss = gaussHermiteData(30)
weights = gauss$w
nodes = gauss$x
## (2.a) age + BMI
# sourceCpp('Code/Rcpp/rcpp_fun_multi.cpp')
# data = read.csv("multiwithout.csv", header= TRUE, sep = ",")
# multiresult = read.csv("Data/lognormal_real/multi_without_beyond10.csv", sep=",", header = FALSE)
# boot.seed = t(read.table("Data/Seed/seed_lognormal_multi_beyond10.txt"))
## (2.b) BMI
# sourceCpp('Code/Rcpp/rcpp_fun_multi_bmi.cpp')
# data = read.csv("multiwithout.csv", header= TRUE, sep = ",")
# multiresult = read.csv("Data/lognormal_real/multi_beyond10_bmi_Q10.csv", sep=",", header = FALSE)
# boot.seed = t(read.table("Data/Seed/seed_lognormal_multi_beyond10_bmi.txt"))
## (2.c) age
# sourceCpp('Code/Rcpp/rcpp_fun_multi_age.cpp.cpp')
# data = read.csv("multiwithout.csv", header= TRUE, sep = ",")
# multiresult = read.csv("Data/lognormal_real/multi_beyond10_age_Q10.csv", sep=",", header = FALSE)
# boot.seed = t(read.table("Data/Seed/seed_lognormal_multi_beyond10_age.txt"))

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
## (2.a) age + BMI
#main.est = as.numeric(read.csv("Data/lognormal_real/noboot_multi_beyond10_Q10_8169573.csv",
#                               sep=",", header = FALSE))
## (2.b) BMI
#main.est = as.numeric(read.csv("Data/lognormal_real/noboot_multi_beyond10_bmi_Q10_8169573.csv",
#                               sep=",", header = FALSE))
## (2.c) age 
#main.est = as.numeric(read.csv("Data/lognormal_real/noboot_multi_beyond10_age_Q10_8169573.csv",
#                               sep=",", header = FALSE))
main.est = main.est[-c(1:3)]
frailty.range = seq(0.1, 3, 0.1)
gradient.lognormal.multi  = rep(0, length(frailty.range))
for(r in 1:length(frailty.range)){
  print(r)
  for(i in 1:length(obs.T)){
    deno = c()
    for(q in 1:length(weights)){
      deno[q] = joint(nodes[q], main.est, i-1)*weights[q]
    }
    gradient.lognormal.multi[r] = gradient.lognormal.multi[r] + conditional(-log(frailty.range[r]), main.est, i-1) / max(sum(deno),10^(-20))
  }
}
## (2.a) age + BMI
# write.table(gradient.lognormal.multi/n, "Data/lognormal_real/gradient_lognormal_multi.txt", sep="\t")
## (2.b) BMI
# write.table(gradient.lognormal.multi/n, "Data/lognormal_real/gradient_lognormal_multi_bmi.txt", sep="\t")
## (2.c) age
# write.table(gradient.lognormal.multi/n, "Data/lognormal_real/gradient_lognormal_multi_age.txt", sep="\t")
orgi.T = obs.T
orgi.Delta = obs.Delta
orgi.U = obs.U
orgi.bmi1 = obs.bmi1
orgi.bmi2 = obs.bmi2
orgi.bmi3 = obs.bmi3
orgi.age = obs.age
orgi.age2 = obs.age2
orgi.Bsplines = Bsplines
orgi.Isplines = Isplines

rectime = c()
for(i in 1:length(boot.seed)){
  rectime[i] = proc.time()-ptm
  main.est = as.numeric(multiresult[which(multiresult[,1] == boot.seed[i]),4:ncol(multiresult)])
  set.seed(boot.seed[i])
  n = nrow(data)
  indbt = sample(1:n, replace=TRUE)
  
  d.time = orgi.T[indbt]
  Delta = orgi.Delta[indbt]
  U = orgi.U[indbt]
  bmi1 = orgi.bmi1[indbt]
  bmi2 = orgi.bmi2[indbt]
  bmi3 = orgi.bmi3[indbt]
  age = orgi.age[indbt]
  age2 = orgi.age2[indbt]
  
  ind = order(d.time, decreasing = FALSE)
  obs.T = d.time[ind]
  obs.Delta = Delta[ind]
  obs.U = U[ind]
  obs.bmi1 = bmi1[ind]
  obs.bmi2 = bmi2[ind]
  obs.bmi3 = bmi3[ind]
  obs.age = age[ind]
  obs.age2 = age2[ind]
  
  Bsplines = orgi.Bsplines[indbt[ind],]
  Isplines = orgi.Isplines[indbt[ind],]
  
  frailty.range = seq(0.1, 3, 0.1)
  
  gradient.lognormal.multi  = rep(0, length(frailty.range))
  for(r in 1:length(frailty.range)){
    for(i in 1:length(obs.T)){
      deno = c()
      for(q in 1:length(weights)){
        deno[q] = joint(nodes[q], main.est, i-1)*weights[q]
      }
      gradient.lognormal.multi[r] = gradient.lognormal.multi[r] + conditional(-log(frailty.range[r]), main.est, i-1) / max(sum(deno),10^(-20))
    }
  }
  
  result[i] = gradient.lognormal.multi / n
  rectime[i] = proc.time()-ptm
}

## (2.a) age + BMI
#write.table(cbind(boot.seed, rep(n, length(boot.seed)), rectime, result),
#            file=paste('Data/lognormal_grad_multi_beyond10_Q10', '.csv',sep=''),sep=',', row.names=F, col.names = F)
## (2.b) BMI
#write.table(cbind(boot.seed, rep(n, length(boot.seed)), rectime, result),
#            file=paste('Data/lognormal_grad_multi_beyond10_bmi_Q10', '.csv',sep=''),sep=',', row.names=F, col.names = F)
## (2.c) age
#write.table(cbind(boot.seed, rep(n, length(boot.seed)), rectime, result),
#            file=paste('Data/lognormal_grad_multi_beyond10_age_Q10', '.csv',sep=''),sep=',', row.names=F, col.names = F)


boot.gradient.lognormal.null = as.matrix(read.csv("Data/lognormal_real/lognormal_grad_null_without_Q10.csv", sep=",", header = FALSE))
boot.gradient.lognormal.null.bmi = as.matrix(read.csv("Data/lognormal_real/lognormal_grad_null_bmi_Q10.csv", sep=",", header = FALSE))
boot.gradient.lognormal.null.age = as.matrix(read.csv("Data/lognormal_real/lognormal_grad_null_age_Q10.csv", sep=",", header = FALSE))

gradient.lognormal.null = read.table("Data/lognormal_real/gradient_lognormal_null.txt", sep="\t")
pdf("Figure/gradient_lognormal_null.pdf", height = 8, width = 14)
par(mfrow = c(1,1),   mar = c(7,12,5,3),  cex.lab = 2.8, 
    cex.main = 3, cex.axis = 2, tcl = 0.5)
plot(frailty.range, as.numeric(t(gradient.lognormal.null)), xlab = expression(sigma),
     ylab = "Gradient function", 
     main = "Nulliparous (BMI + Age), Log-normal",
     ylim = c(0.98, 1.02), pch = 20, cex = 1, type = "b",
     lwd = 4,
     mgp = c(5,1,0))
lines(frailty.range, apply(boot.gradient.lognormal.null, 2, findU)[-c(1:3)],
      lty = 2, lwd = 2)
lines(frailty.range, apply(boot.gradient.lognormal.null, 2, findL)[-c(1:3)],
      lty = 2, lwd = 2)
abline(h = 1, lwd = 2, col = "red")
dev.off()


gradient.lognormal.null.bmi = read.table("Data/lognormal_real/gradient_lognormal_null_bmi.txt", sep="\t")
pdf("Figure/gradient_lognormal_null_bmi.pdf", height = 8, width = 14)
par(mfrow = c(1,1),   mar = c(7,12,5,3),  cex.lab = 2.8, 
    cex.main = 3, cex.axis = 2, tcl = 0.5)
plot(frailty.range, as.numeric(t(gradient.lognormal.null.bmi)), xlab = expression(sigma),
     ylab = "Gradient function", 
     main = "Nulliparous (BMI), Log-normal",
     ylim = c(0.98, 1.02), pch = 20, cex = 1, type = "b",
     lwd = 4,
     mgp = c(5,1,0))
lines(frailty.range, apply(boot.gradient.lognormal.null.bmi, 2, findU)[-c(1:3)],
      lty = 2, lwd = 2)
lines(frailty.range, apply(boot.gradient.lognormal.null.bmi, 2, findL)[-c(1:3)],
      lty = 2, lwd = 2)
abline(h = 1, lwd = 2, col = "red")
dev.off()

gradient.lognormal.null.age = read.table("Data/lognormal_real/gradient_lognormal_null_age.txt", sep="\t")
pdf("Figure/gradient_lognormal_null_age.pdf", height = 8, width = 14)
par(mfrow = c(1,1),   mar = c(7,12,5,3),  cex.lab = 2.8, 
    cex.main = 3, cex.axis = 2, tcl = 0.5)
plot(frailty.range, as.numeric(t(gradient.lognormal.null.age)), xlab = expression(sigma),
     ylab = "Gradient function", 
     main = "Nulliparous (Age), Log-normal",
     ylim = c(0.98, 1.02), pch = 20, cex = 1, type = "b",
     lwd = 4,
     mgp = c(5,1,0))
lines(frailty.range, apply(boot.gradient.lognormal.null.age, 2, findU)[-c(1:3)],
      lty = 2, lwd = 2)
lines(frailty.range, apply(boot.gradient.lognormal.null.age, 2, findL)[-c(1:3)],
      lty = 2, lwd = 2)
abline(h = 1, lwd = 2, col = "red")
dev.off()


###
boot.gradient.lognormal.multi = as.matrix(read.csv("data/lognormal_real/multi/lognormal_grad_multi_beyond10_Q10.csv", sep=",", header = FALSE))
boot.gradient.lognormal.multi.bmi = as.matrix(read.csv("data/lognormal_real/multi/lognormal_grad_multi_beyond10_bmi_Q10.csv", sep=",", header = FALSE))
boot.gradient.lognormal.multi.age = as.matrix(read.csv("data/lognormal_real/multi/lognormal_grad_multi_beyond10_age_Q10.csv", sep=",", header = FALSE))

gradient.lognormal.multi = read.table("Data/lognormal_real/gradient_lognormal_multi.txt", sep="\t")
pdf("Figure/gradient_lognormal_multi.pdf", height = 8, width = 14)
par(mfrow = c(1,1),   mar = c(7,12,5,3),  cex.lab = 2.8, 
    cex.main = 3, cex.axis = 2, tcl = 0.5)
plot(frailty.range, as.numeric(t(gradient.lognormal.multi)), xlab = expression(sigma),
     ylab = "Gradient function", 
     main = "Multiparous (BMI + Age), Log-normal",
     ylim = c(0.98, 1.02), pch = 20, cex = 1, type = "b",
     lwd = 4,
     mgp = c(5,1,0))
lines(frailty.range, apply(boot.gradient.lognormal.multi, 2, findU)[-c(1:3)],
      lty = 2, lwd = 2)
lines(frailty.range, apply(boot.gradient.lognormal.multi, 2, findL)[-c(1:3)],
      lty = 2, lwd = 2)
abline(h = 1, lwd = 2, col = "red")
dev.off()


gradient.lognormal.multi.age = read.table("Data/lognormal_real/gradient_lognormal_multi_age.txt", sep="\t")
pdf("Figure/gradient_lognormal_multi_bmi.pdf", height = 8, width = 14)
par(mfrow = c(1,1),   mar = c(7,12,5,3),  cex.lab = 2.8, 
    cex.main = 3, cex.axis = 2, tcl = 0.5)
plot(frailty.range, as.numeric(t(gradient.lognormal.multi.bmi)), xlab = expression(sigma),
     ylab = "Gradient function", 
     main = "Multiparous (BMI), Log-normal",
     ylim = c(0.98, 1.02), pch = 20, cex = 1, type = "b",
     lwd = 4,
     mgp = c(5,1,0))
lines(frailty.range, apply(boot.gradient.lognormal.multi.bmi, 2, findU)[-c(1:3)],
      lty = 2, lwd = 2)
lines(frailty.range, apply(boot.gradient.lognormal.multi.bmi, 2, findL)[-c(1:3)],
      lty = 2, lwd = 2)
abline(h = 1, lwd = 2, col = "red")
dev.off()

gradient.lognormal.multi.age = read.table("Data/lognormal_real/gradient_lognormal_multi_age.txt", sep="\t")
pdf("Figure/gradient_lognormal_multi_age.pdf", height = 8, width = 14)
par(mfrow = c(1,1),   mar = c(7,12,5,3),  cex.lab = 2.8, 
    cex.main = 3, cex.axis = 2, tcl = 0.5)
plot(frailty.range, as.numeric(t(gradient.lognormal.multi.age)), xlab = expression(sigma),
     ylab = "Gradient function", 
     main = "Multiparous (Age), Log-normal",
     ylim = c(0.98, 1.02), pch = 20, cex = 1, type = "b",
     lwd = 4,
     mgp = c(5,1,0))
lines(frailty.range, apply(boot.gradient.lognormal.multi.age, 2, findU)[-c(1:3)],
      lty = 2, lwd = 2)
lines(frailty.range, apply(boot.gradient.lognormal.multi.age, 2, findL)[-c(1:3)],
      lty = 2, lwd = 2)
abline(h = 1, lwd = 2, col = "red")
dev.off()



## Nulliparous (Weibull) ##
gauss = gauss.quad(30, kind="laguerre")
weights = gauss$weights
nodes = gauss$nodes
## (3.a) age + BMI
# sourceCpp('Code/Rcpp/rcpp_null.cpp')
# data = read.csv("nullwithout.csv", header= TRUE, sep = ",")
# nullresult = read.csv("Data/weibull_real/weibull_null_without_Q10.csv", sep=",", header = FALSE)
# boot.seed = t(read.table("Data/Seed/seed_weibull_null.txt"))
## (3.b) BMI
# sourceCpp('Code/Rcpp/rcpp_null_bmi.cpp')
# data = read.csv("nullwithout.csv", header= TRUE, sep = ",")
# nullresult = read.csv("Data/weibull_real/weibull_null_bmi_Q10.csv", sep=",", header = FALSE)
# boot.seed = t(read.table("Data/Seed/seed_weibull_null_bmi.txt"))
## (3.c) age
# sourceCpp('Code/Rcpp/rcpp_null_age.cpp')
# data = read.csv("nullwithout.csv", header= TRUE, sep = ",")
# nullresult = read.csv("Data/weibull_real/weibull_null_age_Q10.csv", sep=",", header = FALSE)
# boot.seed = t(read.table("Data/Seed/seed_weibull_null_age.txt"))

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

## (3.a) age + BMI
#main.est = as.numeric(read.csv("Data/weibull_real/noboot_weibull_null_Q10_7978997.csv",
#                               sep=",", header = FALSE))
## (3.b) BMI
#main.est = as.numeric(read.csv("Data/weibull_real/noboot_weibull_null_bmi_Q10_7978997.csv",
#sep=",", header = FALSE))
## (3.c) age 
#main.est = as.numeric(read.csv("Data/weibull_real/noboot_weibull_null_age_Q10_7978997.csv",
#sep=",", header = FALSE))
main.est = main.est[-c(1:3)]
frailty.range = seq(0.1, 3, 0.1)
gradient.weibull.null  = rep(0, length(frailty.range))
for(r in 1:length(frailty.range)){
  print(r)
  for(i in 1:length(obs.T)){
    deno = c()
    for(q in 1:length(weights)){
      deno[q] = joint(nodes[q], main.est, i-1)*weights[q]
    }
    gradient.weibull.null[r] = gradient.weibull.null[r] + conditional(frailty.range[r], main.est, i-1) / max(sum(deno),10^(-20))
  }
}

## (3.a) age + BMI
#write.table(gradient.weibull.null/n, "Data/weibull_real/gradient_weibull_null.txt", sep="\t")
## (3.b) BMI
#write.table(gradient.weibull.null/n, "Data/weibull_real/gradient_weibull_null_bmi.txt", sep="\t")
## (3.c) age 
#write.table(gradient.weibull.null/n, "Data/weibull_real/gradient_weibull_null_age.txt", sep="\t")
orgi.T = obs.T
orgi.Delta = obs.Delta
orgi.U = obs.U
orgi.bmi1 = obs.bmi1
orgi.bmi2 = obs.bmi2
orgi.bmi3 = obs.bmi3
orgi.age = obs.age
orgi.age2 = obs.age2
orgi.Bsplines = Bsplines
orgi.Isplines = Isplines

rectime = c()
for(i in 1:length(boot.seed)){
  rectime[i] = proc.time()-ptm
  main.est = as.numeric(nullresult[which(nullresult[,1] == boot.seed[i]),4:ncol(nullresult)])
  set.seed(boot.seed[i])
  n = nrow(data)
  indbt = sample(1:n, replace=TRUE)
  
  d.time = orgi.T[indbt]
  Delta = orgi.Delta[indbt]
  U = orgi.U[indbt]
  bmi1 = orgi.bmi1[indbt]
  bmi2 = orgi.bmi2[indbt]
  bmi3 = orgi.bmi3[indbt]
  age = orgi.age[indbt]
  age2 = orgi.age2[indbt]
  
  ind = order(d.time, decreasing = FALSE)
  obs.T = d.time[ind]
  obs.Delta = Delta[ind]
  obs.U = U[ind]
  obs.bmi1 = bmi1[ind]
  obs.bmi2 = bmi2[ind]
  obs.bmi3 = bmi3[ind]
  obs.age = age[ind]
  obs.age2 = age2[ind]
  
  Bsplines = orgi.Bsplines[indbt[ind],]
  Isplines = orgi.Isplines[indbt[ind],]
  
  frailty.range = seq(0.1, 3, 0.1)
  
  gradient.fun  = rep(0, length(frailty.range))
  for(r in 1:length(frailty.range)){
    
    for(i in 1:length(obs.T)){
      deno = c()
      for(q in 1:length(weights)){
        deno[q] = joint(nodes[q], main, i-1)*weights[q]
      }
      
      gradient.fun[r] = gradient.fun[r] + 
        conditional(frailty.range[r], main, i-1) / max(sum(deno),10^(-20))
    }
  }
  
  result[i] = gradient.fun / n
  rectime[i] = proc.time()-ptm
}

## (3.a) age + BMI
#write.table(cbind(boot.seed, rep(n, length(boot.seed)), rectime, result),
#            file=paste('Data/weibull_real/weibull_grad_null_Q10', '.csv',sep=''),sep=',', row.names=F, col.names = F)
## (3.b) BMI
#write.table(cbind(boot.seed, rep(n, length(boot.seed)), rectime, result),
#            file=paste('Data/weibull_real/weibull_grad_null_bmi_Q10', '.csv',sep=''),sep=',', row.names=F, col.names = F)
## (3.c) age 
#write.table(cbind(boot.seed, rep(n, length(boot.seed)), rectime, result),
#            file=paste('Data/weibull_real/weibull_grad_null_age_Q10', '.csv',sep=''),sep=',', row.names=F, col.names = F)



## Multiparous (Weibull) ## 
gauss = gauss.quad(30, kind="laguerre")
weights = gauss$weights
nodes = gauss$nodes
## (4.a) age + BMI
# sourceCpp('Code/Rcpp/rcpp_multi.cpp')
# data = read.csv("multiwithout.csv", header= TRUE, sep = ",")
# multiresult = read.csv("Data/weibull_real/weibull_multi_beyond10_Q10.csv", sep=",", header = FALSE)
# boot.seed = t(read.table("Data/Seed/seed_weibull_multi_beyond10.txt"))
## (4.b) BMI
# sourceCpp('Code/Rcpp/rcpp_multi_bmi.cpp')
# data = read.csv("multiwithout.csv", header= TRUE, sep = ",")
# multiresult = read.csv("Data/weibull_real/weibull_multi_beyond10_bmi_Q10.csv", sep=",", header = FALSE)
# boot.seed = t(read.table("Data/Seed/seed_weibull_multi_beyond10_bmi.txt"))
## (4.c) age
# sourceCpp('Code/Rcpp/rcpp_multi_age.cpp')
# data = read.csv("multiwithout.csv", header= TRUE, sep = ",")
# multiresult = read.csv("Data/weibull_real/weibull_multi_beyond10_age_Q10.csv", sep=",", header = FALSE)
# boot.seed = t(read.table("Data/Seed/seed_weibull_multi_beyond10_age.txt"))

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
## (4.a) age + BMI
#main.est = as.numeric(read.csv("Data/weibull_real/noboot_weibull_multi_beyond10_Q10_7978997.csv",
#                               sep=",", header = FALSE))
## (4.b) BMI
#main.est = as.numeric(read.csv("data/weibull_real/noboot_weibull_multi_beyond10_bmi_Q10_7978997.csv",
#sep=",", header = FALSE))
## (4.c) age 
#main.est = as.numeric(read.csv("data/weibull_real/noboot_weibull_multi_beyond10_age_Q10_7978997.csv",
#sep=",", header = FALSE))
main.est = main.est[-c(1:3)]
frailty.range = seq(0.1, 3, 0.1)
gradient.weibull.multi  = rep(0, length(frailty.range))
for(r in 1:length(frailty.range)){
  print(r)
  for(i in 1:length(obs.T)){
    deno = c()
    for(q in 1:length(weights)){
      deno[q] = joint(nodes[q], main.est, i-1)*weights[q]
    }
    gradient.weibull.multi[r] = gradient.weibull.multi[r] + conditional(frailty.range[r], main.est, i-1) / max(sum(deno),10^(-20))
  }
}
## (4.a) age + BMI
#write.table(gradient.weibull.multi/n, "Data/weibull_real/gradient_weibull_multi.txt", sep="\t")
## (4.b) BMI
#write.table(gradient.weibull.multi/n, "Data/weibull_real/gradient_weibull_multi_bmi.txt", sep="\t")
## (4.c) age
#write.table(gradient.weibull.multi/n, "Data/weibull_real/gradient_weibull_multi_age.txt", sep="\t")
orgi.T = obs.T
orgi.Delta = obs.Delta
orgi.U = obs.U
orgi.bmi1 = obs.bmi1
orgi.bmi2 = obs.bmi2
orgi.bmi3 = obs.bmi3
orgi.age = obs.age
orgi.age2 = obs.age2
orgi.Bsplines = Bsplines
orgi.Isplines = Isplines

rectime = c()
for(i in 1:length(boot.seed)){
  rectime[i] = proc.time()-ptm
  main.est = as.numeric(multiresult[which(multiresult[,1] == boot.seed[i]),4:ncol(multiresult)])
  set.seed(boot.seed[i])
  n = nrow(data)
  indbt = sample(1:n, replace=TRUE)
  
  d.time = orgi.T[indbt]
  Delta = orgi.Delta[indbt]
  U = orgi.U[indbt]
  bmi1 = orgi.bmi1[indbt]
  bmi2 = orgi.bmi2[indbt]
  bmi3 = orgi.bmi3[indbt]
  age = orgi.age[indbt]
  age2 = orgi.age2[indbt]
  
  ind = order(d.time, decreasing = FALSE)
  obs.T = d.time[ind]
  obs.Delta = Delta[ind]
  obs.U = U[ind]
  obs.bmi1 = bmi1[ind]
  obs.bmi2 = bmi2[ind]
  obs.bmi3 = bmi3[ind]
  obs.age = age[ind]
  obs.age2 = age2[ind]
  
  Bsplines = orgi.Bsplines[indbt[ind],]
  Isplines = orgi.Isplines[indbt[ind],]
  
  frailty.range = seq(0.1, 3, 0.1)
  
  gradient.fun  = rep(0, length(frailty.range))
  for(r in 1:length(frailty.range)){
    
    for(i in 1:length(obs.T)){
      deno = c()
      for(q in 1:length(weights)){
        deno[q] = joint(nodes[q], main, i-1)*weights[q]
      }
      
      gradient.fun[r] = gradient.fun[r] + 
        conditional(frailty.range[r], main, i-1) / max(sum(deno),10^(-20))
    }
  }
  
  result[i] = gradient.fun / n
  rectime[i] = proc.time()-ptm
}
## (4.a) age + BMI
#write.table(cbind(boot.seed, rep(n, length(boot.seed)), rectime, result),
#            file=paste('Data/weibull_real/weibull_grad_multi_Q10', '.csv',sep=''),sep=',', row.names=F, col.names = F)
## (4.b) BMI
#write.table(cbind(boot.seed, rep(n, length(boot.seed)), rectime, result),
#            file=paste('Data/weibull_real/weibull_grad_multi_bmi_Q10', '.csv',sep=''),sep=',', row.names=F, col.names = F)
## (4.c) age 
#write.table(cbind(boot.seed, rep(n, length(boot.seed)), rectime, result),
#            file=paste('Data/weibull_real/weibull_grad_multi_age_Q10', '.csv',sep=''),sep=',', row.names=F, col.names = F)

boot.gradient.weibull.null = as.matrix(read.csv("Data/weibull_real/null/weibull_grad_null_Q10.csv", sep=",", header = FALSE))
boot.gradient.weibull.null.bmi = as.matrix(read.csv("Data/weibull_real/null/weibull_grad_null_bmi_Q10.csv", sep=",", header = FALSE))
boot.gradient.weibull.null.age = as.matrix(read.csv("Data/weibull_real/null/weibull_grad_null_age_Q10.csv", sep=",", header = FALSE))

gradient.weibull.null = read.table("Data/weibull_real/gradient_weibull_null.txt", sep="\t")
pdf("Figure/gradient_weibull_null.pdf", height = 8, width = 14)
par(mfrow = c(1,1),   mar = c(7,12,5,3),  cex.lab = 2.8, 
    cex.main = 3, cex.axis = 2, tcl = 0.5)
plot(frailty.range, as.numeric(t(gradient.weibull.null)), xlab = expression(gamma),
     ylab = "Gradient function", 
     main = "Nulliparous (BMI + Age), Weibull",
     ylim = c(0.95, 1.05), pch = 20, cex = 1, type = "b",
     lwd = 4,
     mgp = c(5,1,0))
lines(frailty.range, apply(boot.gradient.weibull.null, 2, findU)[-c(1:3)],
      lty = 2, lwd = 2)
lines(frailty.range, apply(boot.gradient.weibull.null, 2, findL)[-c(1:3)],
      lty = 2, lwd = 2)
abline(h = 1, lwd = 2, col = "red")
dev.off()

gradient.weibull.null.bmi = read.table("Data/weibull_real/gradient_weibull_null_bmi.txt", sep="\t")
pdf("Figure/gradient_weibull_null_bmi.pdf", height = 8, width = 14)
par(mfrow = c(1,1),   mar = c(7,12,5,3),  cex.lab = 2.8, 
    cex.main = 3, cex.axis = 2, tcl = 0.5)
plot(frailty.range, as.numeric(t(gradient.weibull.null.bmi)), xlab = expression(gamma),
     ylab = "Gradient function", 
     main = "Nulliparous (BMI), Weibull",
     ylim = c(0.95, 1.05), pch = 20, cex = 1, type = "b",
     lwd = 4,
     mgp = c(5,1,0))
lines(frailty.range, apply(boot.gradient.weibull.null.bmi, 2, findU)[-c(1:3)],
      lty = 2, lwd = 2)
lines(frailty.range, apply(boot.gradient.weibull.null.bmi, 2, findL)[-c(1:3)],
      lty = 2, lwd = 2)
abline(h = 1, lwd = 2, col = "red")
dev.off()

gradient.weibull.null.age = read.table("Data/weibull_real/gradient_weibull_null_age.txt", sep="\t")
pdf("Figure/gradient_weibull_null_age.pdf", height = 8, width = 14)
par(mfrow = c(1,1),   mar = c(7,12,5,3),  cex.lab = 2.8, 
    cex.main = 3, cex.axis = 2, tcl = 0.5)
plot(frailty.range, as.numeric(t(gradient.weibull.null.age)), xlab = expression(gamma),
     ylab = "Gradient function", 
     main = "Nulliparous (Age), Weibull",
     ylim = c(0.95, 1.05), pch = 20, cex = 1, type = "b",
     lwd = 4,
     mgp = c(5,1,0))
lines(frailty.range, apply(boot.gradient.weibull.null.age, 2, findU)[-c(1:3)],
      lty = 2, lwd = 2)
lines(frailty.range, apply(boot.gradient.weibull.null.age, 2, findL)[-c(1:3)],
      lty = 2, lwd = 2)
abline(h = 1, lwd = 2, col = "red")
dev.off()


###
boot.gradient.weibull.multi = as.matrix(read.csv("data/weibull_real/multi/weibull_grad_multi_Q10.csv", sep=",", header = FALSE))
boot.gradient.weibull.multi.bmi = as.matrix(read.csv("data/weibull_real/multi/weibull_grad_multi_bmi_Q10.csv", sep=",", header = FALSE))
boot.gradient.weibull.multi.age = as.matrix(read.csv("data/weibull_real/multi/weibull_grad_multi_age_Q10.csv", sep=",", header = FALSE))

gradient.weibull.multi = read.table("Data/weibull_real/gradient_weibull_multi.txt", sep="\t")
pdf("Figure/gradient_weibull_multi.pdf", height = 8, width = 14)
par(mfrow = c(1,1),   mar = c(7,12,5,3),  cex.lab = 2.8, 
    cex.main = 3, cex.axis = 2, tcl = 0.5)
plot(frailty.range, as.numeric(t(gradient.weibull.multi)), xlab = expression(gamma),
     ylab = "Gradient function", 
     main = "Multiparous (BMI + Age), Weibull",
     ylim = c(0.95, 1.05), pch = 20, cex = 1, type = "b",
     lwd = 4,
     mgp = c(5,1,0))
lines(frailty.range, apply(boot.gradient.weibull.multi, 2, findU)[-c(1:3)],
      lty = 2, lwd = 2)
lines(frailty.range, apply(boot.gradient.weibull.multi, 2, findL)[-c(1:3)],
      lty = 2, lwd = 2)
abline(h = 1, lwd = 2, col = "red")
dev.off()


gradient.weibull.multi.bmi = read.table("Data/weibull_real/gradient_weibull_multi_bmi.txt", sep="\t")
pdf("Figure/gradient_weibull_multi_bmi.pdf", height = 8, width = 14)
par(mfrow = c(1,1),   mar = c(7,12,5,3),  cex.lab = 2.8, 
    cex.main = 3, cex.axis = 2, tcl = 0.5)
plot(frailty.range, as.numeric(t(gradient.weibull.multi.bmi)), xlab = expression(gamma),
     ylab = "Gradient function", 
     main = "Multiparous (BMI), Weibull",
     ylim = c(0.95, 1.05), pch = 20, cex = 1, type = "b",
     lwd = 4,
     mgp = c(5,1,0))
lines(frailty.range, apply(boot.gradient.weibull.multi.bmi, 2, findU)[-c(1:3)],
      lty = 2, lwd = 2)
lines(frailty.range, apply(boot.gradient.weibull.multi.bmi, 2, findL)[-c(1:3)],
      lty = 2, lwd = 2)
abline(h = 1, lwd = 2, col = "red")
dev.off()


gradient.weibull.multi.age = read.table("Data/weibull_real/gradient_weibull_multi_age.txt", sep="\t")
pdf("Figure/gradient_weibull_multi_age.pdf", height = 8, width = 14)
par(mfrow = c(1,1),   mar = c(7,12,5,3),  cex.lab = 2.8, 
    cex.main = 3, cex.axis = 2, tcl = 0.5)
plot(frailty.range, as.numeric(t(gradient.weibull.multi.age)), xlab = expression(gamma),
     ylab = "Gradient function", 
     main = "Multiparous (Age), Weibull",
     ylim = c(0.95, 1.05), pch = 20, cex = 1, type = "b",
     lwd = 4,
     mgp = c(5,1,0))
lines(frailty.range, apply(boot.gradient.weibull.multi.age, 2, findU)[-c(1:3)],
      lty = 2, lwd = 2)
lines(frailty.range, apply(boot.gradient.weibull.multi.age, 2, findL)[-c(1:3)],
      lty = 2, lwd = 2)
abline(h = 1, lwd = 2, col = "red")
dev.off()
