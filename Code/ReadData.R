library(sas7bdat)
library(statmod)
### read rawdata
#dat = read.sas7bdat("Data/forndeath1.sas7bdat") # restricted use
dat.null.with = dat[dat$epidural == 2 & dat$cohort == "Nulliparous",] 
dat.multi.with = dat[dat$epidural == 2 & dat$cohort == "Multiparous",] 
dat.null.without = dat[dat$epidural == 1 & dat$cohort == "Nulliparous",] 
dat.multi.without = dat[dat$epidural == 1 & dat$cohort == "Multiparous",] 

### nulliparous women without epidural
d.time = (dat.null.without$stage2)
X1 = dat.null.without$admbmi
X2 = dat.null.without$Momage
Delta = ifelse(dat.null.without$spontvd == 2, 1, 0)
Delta = ifelse(dat.null.without$cesdli == 2, 2, Delta)
Delta = ifelse(dat.null.without$opervd == 2, 3, Delta)
U1 = ifelse(dat.null.without$babyvar1 == 2, 1, 0)
U2 = ifelse(dat.null.without$momvar1 == 2, 1, 0)
U = ifelse( U1 + U2 == 0, 0, 1) # combine neonatal and maternal morbidity

data = cbind(d.time, Delta, U1, U2, U, X1, X2)
data = as.data.frame(na.omit(data))
ind = order(data$d.time, decreasing = FALSE)

obs.T = data$d.time[ind]
obs.Delta = data$Delta[ind]
obs.U = data$U[ind]
obs.X1 = data$X1[ind]
obs.bmi = ifelse(obs.X1 < 25, 1, 3)
obs.bmi = ifelse(25 <= obs.X1 & obs.X1 < 30, 2, obs.bmi)

obs.bmi1 = ifelse(obs.bmi == 1, 1, 0) ## reference group
obs.bmi2 = ifelse(obs.bmi == 2, 1, 0)
obs.bmi3 = ifelse(obs.bmi == 3, 1, 0)

obs.X2 = data$X2[ind] # standardize mother's age
obs.age = (obs.X2 - mean(obs.X2)) / sd(obs.X2)

nullwithout = cbind(obs.T, obs.Delta, obs.U, obs.X1, obs.X2, obs.bmi1, obs.bmi2, obs.bmi3, obs.age)
nullwithout = as.data.frame(null.without)

write.table(nullwithout, file = "Data/nullwithout.csv", col.names = TRUE,
            sep = ",", row.names = FALSE)
#nullwithout = read.csv("Data/nullwithout.csv", header= TRUE, sep = ",") # restricted data

pdf("Figure/NullWithout.pdf", width = 10, height = 5)
par(mfrow = c(1,1), cex.main = 3, cex.lab = 3, cex.axis = 2,
    mar = c(7,7,7,3), tcl = 0.5)
plot(density(obs.T[obs.Delta == 1]), xlab = "Time on secondary stage labor (hour)",
     main = "Nulliparous women without epidural", col = rgb(0,0,0,0.5),ylab = "Density",
     xlim = c(0, 4), ylim = c(0,1),  mgp = c(4,1,0), lwd = 4, yaxt = "n", xaxt = "n")
polygon(density(obs.T[obs.Delta == 1]), col= rgb(0,0,0,0.5) , border=rgb(0,0,0,0.5), density = 20, angle = 90)
polygon(density(obs.T[obs.Delta == 1]), col= rgb(0,0,0,0.5) , border=rgb(0,0,0,0.5), density = 20, angle = 180)
axis(side = 2, at = seq(0, 1, 0.2), 
     labels = seq(0, 1, 0.2), 
     tck = 0.05)
axis(side = 1, at = seq(0, 4, 0.5), 
     labels = seq(0, 4, 0.5), 
     tck = 0.05)
lines(density(obs.T[obs.Delta == 2]), col = rgb(0,0,0,0.7), lwd = 4)
polygon(density(obs.T[obs.Delta == 2]), col= rgb(0,0,0,0.7) , border=rgb(0,0,0,0.3))
lines(density(obs.T[obs.Delta == 3]), col = rgb(0,0,0,0.1), lwd = 4)
polygon(density(obs.T[obs.Delta == 3]), col= rgb(0,0,0,0.1) , border=rgb(0,0,0,0.7))
abline(v = 2, lwd = 2, col = "red")
dev.off()



### multiparous women without epidural
d.time = (dat.multi.without$stage2)
X1 = dat.multi.without$admbmi 
X2 = dat.multi.without$Momage
Delta = ifelse(dat.multi.without$spontvd == 2, 1, 0)
Delta = ifelse(dat.multi.without$cesdli == 2, 2, Delta)
Delta = ifelse(dat.multi.without$opervd == 2, 3, Delta)
U1 = ifelse(dat.multi.without$babyvar1 == 2, 1, 0)
U2 = ifelse(dat.multi.without$momvar1 == 2, 1, 0)
U = ifelse( U1 + U2 == 0, 0, 1)

data = cbind(d.time, Delta, U1, U2, U, X1, X2)
data = as.data.frame(na.omit(data))
ind = order(data$d.time, decreasing = FALSE)

obs.T = data$d.time[ind]
obs.Delta = data$Delta[ind]
obs.U = data$U[ind]
obs.X1 = data$X1[ind]
obs.bmi = ifelse(obs.X1 < 25, 1, 3)
obs.bmi = ifelse(25 <= obs.X1 & obs.X1 < 30, 2, obs.bmi)

obs.bmi1 = ifelse(obs.bmi == 1, 1, 0) ## reference group
obs.bmi2 = ifelse(obs.bmi == 2, 1, 0)
obs.bmi3 = ifelse(obs.bmi == 3, 1, 0)

obs.X2 = data$X2[ind]
obs.age = (obs.X2 - mean(obs.X2))/sd(obs.X2)

multiwithout = cbind(obs.T, obs.Delta, obs.U, obs.X1, obs.X2, obs.bmi1, obs.bmi2, obs.bmi3, obs.age)
multiwithout = as.data.frame(multi.without)

write.table(multiwithout, file = "Data/multiwithout.csv", col.names = TRUE,
            sep = ",", row.names = FALSE)


### Figure : distribution of delivery time stratified by parity
pdf("Figure/MultiWithout.pdf", width = 10, height = 5)
par(mfrow = c(1,1), cex.main = 3, cex.lab = 3, cex.axis = 2,
    mar = c(7,7,7,3), tcl = 0.5)
plot(density(obs.T[obs.Delta == 1]), xlab = "Time on secondary stage labor (hour)",
     main = "Multiparous women without epidural", col = rgb(0,0,0,0.5),ylab = "Density",
     xlim = c(0, 4), ylim = c(0,4),  mgp = c(4,1,0), lwd = 4, yaxt = "n", xaxt = "n")
polygon(density(obs.T[obs.Delta == 1]), col= rgb(0,0,0,0.5) , border=rgb(0,0,0,0.5), density = 20, angle = 90)
polygon(density(obs.T[obs.Delta == 1]), col= rgb(0,0,0,0.5) , border=rgb(0,0,0,0.5), density = 20, angle = 180)
axis(side = 2, at = seq(0, 4, 1), 
     labels = seq(0, 4, 1), 
     tck = 0.05)
axis(side = 1, at = seq(0, 4, 0.5), 
     labels = seq(0, 4, 0.5), 
     tck = 0.05)
lines(density(obs.T[obs.Delta == 2]), col = rgb(0,0,0,0.7), lwd = 4)
polygon(density(obs.T[obs.Delta == 2]), col= rgb(0,0,0,0.7) , border=rgb(0,0,0,0.3))
lines(density(obs.T[obs.Delta == 3]), col = rgb(0,0,0,0.2), lwd = 4)
polygon(density(obs.T[obs.Delta == 3]), col= rgb(0,0,0,0.2) , border=rgb(0,0,0,0.7))
abline(v = 1, lwd = 2, col = "red")
legend("topright", c("SVD", "CD", "OVD"),
       col = c(rgb(0,0,0,0.5), rgb(0,0,0,0.7), rgb(0,0,0,0.2)),
       seg.len = 3,
       pch = c(12, 15, 15),  bty = 'n', cex = 2, xpd = NA)
dev.off()



## table
delta.table = matrix(0, nrow = 2, ncol = 6)
delta.table[,1] = c("Nulliparous", "Multiparous")
delta.table[,2] = c("Wihtout", "Without")

table(null.without$obs.Delta); table(null.without$obs.Delta) / nrow(null.without)
table(multi.without$obs.Delta); table(multi.without$obs.Delta) / nrow(multi.without)

delta.table[1,c(3:6)] = c("5,510 (90.8%)", "110 (1.8%)", "448 (7.4%)", "6,068")
delta.table[2,c(3:6)] = c("12,000 (98.1%)", "37 (0.3%)", "119 (1.6%)", "12,236")

colnames(delta.table) = c("Parous", "Epidural", 
                          "spontaneous vaginal delivery (SVD)",
                          "cesarean delivery (CES)",
                          "operational vaginal delivery (OPV)",
                          "Total")

print(xtable(delta.table), row.names = NULL)

## nulliparous women by BMI classification
null.tab = matrix(0, nrow = 4, ncol = 4)
colnames(null.tab) = c("Underweight/Normal", "Overweight", "Obesity", "Total")
rownames(null.tab) = c("SVD", "CS", "OVD", "Total")
svd.bmi1 = length(nullwithout$obs.X2[which(nullwithout$obs.Delta == 1 & nullwithout$obs.bmi1 == 1)])
svd.bmi2 = length(nullwithout$obs.X2[which(nullwithout$obs.Delta == 1 & nullwithout$obs.bmi2 == 1)])
svd.bmi3 = length(nullwithout$obs.X2[which(nullwithout$obs.Delta == 1 & nullwithout$obs.bmi3 == 1)])
svd.allbmi = length(nullwithout$obs.X2[which(nullwithout$obs.Delta == 1)])
cs.bmi1 = length(nullwithout$obs.X2[which(nullwithout$obs.Delta == 2 & nullwithout$obs.bmi1 == 1)])
cs.bmi2 = length(nullwithout$obs.X2[which(nullwithout$obs.Delta == 2 & nullwithout$obs.bmi2 == 1)])
cs.bmi3 = length(nullwithout$obs.X2[which(nullwithout$obs.Delta == 2 & nullwithout$obs.bmi3 == 1)])
cs.allbmi = length(nullwithout$obs.X2[which(nullwithout$obs.Delta == 2)])
ovd.bmi1 = length(nullwithout$obs.X2[which(nullwithout$obs.Delta == 3 & nullwithout$obs.bmi1 == 1)])
ovd.bmi2 = length(nullwithout$obs.X2[which(nullwithout$obs.Delta == 3 & nullwithout$obs.bmi2 == 1)])
ovd.bmi3 = length(nullwithout$obs.X2[which(nullwithout$obs.Delta == 3 & nullwithout$obs.bmi3 == 1)])
ovd.allbmi = length(nullwithout$obs.X2[which(nullwithout$obs.Delta == 3)])


bmi1.all = length(nullwithout$obs.X2[which(nullwithout$obs.bmi1 == 1)])
bmi2.all = length(nullwithout$obs.X2[which(nullwithout$obs.bmi2 == 1)])
bmi3.all = length(nullwithout$obs.X2[which(nullwithout$obs.bmi3 == 1)])
alls = length(nullwithout$obs.X2)

null.tab[1,] = c(svd.bmi1, svd.bmi2, svd.bmi3, svd.allbmi)
null.tab[2,] = c(cs.bmi1, cs.bmi2, cs.bmi3, cs.allbmi)
null.tab[3,] = c(ovd.bmi1, ovd.bmi2, ovd.bmi3, ovd.allbmi)
null.tab[4,] = c(bmi1.all, bmi2.all, bmi3.all, alls)


## divide multiparous women by delivery time (< 10 min, >= 10 min)
data = read.csv("data/multiwithout.csv", header= TRUE, sep = ",")
ind1 = which(data$obs.T*60 <= 10)
ind2 = which(data$obs.T*60 > 10)
data1 = data[ind1,]; data2 = data[ind2,]
## bmi
data1.tab.bmi = matrix(0, nrow = 5, ncol = 4)
data1.tab.bmi[1,] = c(sum(data1$obs.Delta == 1 & data1$obs.bmi1 == 1),
                      sum(data1$obs.Delta == 1 & data1$obs.bmi2 == 1),
                      sum(data1$obs.Delta == 1 & data1$obs.bmi3 == 1),
                      sum(data1$obs.Delta == 1))
data1.tab.bmi[2,] = c(sum(data1$obs.Delta != 1 & data1$obs.bmi1 == 1),
                      sum(data1$obs.Delta != 1 & data1$obs.bmi2 == 1),
                      sum(data1$obs.Delta != 1 & data1$obs.bmi3 == 1),
                      sum(data1$obs.Delta != 1))
data1.tab.bmi[3,] = c(sum(data1$obs.bmi1 == 1),
                      sum(data1$obs.bmi2 == 1),
                      sum(data1$obs.bmi3 == 1),
                      nrow(data1))
data1.tab.bmi[4,] = c(sum(data1$obs.U == 1 & data1$obs.bmi1 == 1),
                      sum(data1$obs.U == 1 & data1$obs.bmi2 == 1),
                      sum(data1$obs.U == 1 & data1$obs.bmi3 == 1),
                      sum(data1$obs.U == 1))
data1.tab.bmi[5,] = c(sum(data1$obs.U == 0 & data1$obs.bmi1 == 1),
                      sum(data1$obs.U == 0 & data1$obs.bmi2 == 1),
                      sum(data1$obs.U == 0 & data1$obs.bmi3 == 1),
                      sum(data1$obs.U == 0))
print(xtable(formatC(data1.tab.bmi, format = "f", digits = 0)))
### age
data1.tab.age = matrix(0, nrow = 4, ncol = 3)
data1.tab.age[1,] = quantile(data1$obs.X1[data1$obs.Delta == 1] , c(0.25, 0.50, 0.75))
data1.tab.age[2,] = quantile(data1$obs.X1[data1$obs.Delta != 1] , c(0.25, 0.50, 0.75))
data1.tab.age[3,] = quantile(data1$obs.X1[data1$obs.U == 1] , c(0.25, 0.50, 0.75))
data1.tab.age[4,] = quantile(data1$obs.X1[data1$obs.U == 0] , c(0.25, 0.50, 0.75))
print(xtable(formatC(data1.tab.age, format = "f", digits = 2)))

## bmi
data2.tab.bmi = matrix(0, nrow = 5, ncol = 4)
data2.tab.bmi[1,] = c(sum(data2$obs.Delta == 1 & data2$obs.bmi1 == 1),
                      sum(data2$obs.Delta == 1 & data2$obs.bmi2 == 1),
                      sum(data2$obs.Delta == 1 & data2$obs.bmi3 == 1),
                      sum(data2$obs.Delta == 1))
data2.tab.bmi[2,] = c(sum(data2$obs.Delta != 1 & data2$obs.bmi1 == 1),
                      sum(data2$obs.Delta != 1 & data2$obs.bmi2 == 1),
                      sum(data2$obs.Delta != 1 & data2$obs.bmi3 == 1),
                      sum(data2$obs.Delta != 1))
data2.tab.bmi[3,] = c(sum(data2$obs.bmi1 == 1),
                      sum(data2$obs.bmi2 == 1),
                      sum(data2$obs.bmi3 == 1),
                      nrow(data2))
data2.tab.bmi[4,] = c(sum(data2$obs.U == 1 & data2$obs.bmi1 == 1),
                      sum(data2$obs.U == 1 & data2$obs.bmi2 == 1),
                      sum(data2$obs.U == 1 & data2$obs.bmi3 == 1),
                      sum(data2$obs.U == 1))
data2.tab.bmi[5,] = c(sum(data2$obs.U == 0 & data2$obs.bmi1 == 1),
                      sum(data2$obs.U == 0 & data2$obs.bmi2 == 1),
                      sum(data2$obs.U == 0 & data2$obs.bmi3 == 1),
                      sum(data2$obs.U == 0))
print(xtable(formatC(data2.tab.bmi, format = "f", digits = 0)))
### age
data2.tab.age = matrix(0, nrow = 4, ncol = 3)
data2.tab.age[1,] = quantile(data2$obs.X1[data2$obs.Delta == 1] , c(0.25, 0.50, 0.75))
data2.tab.age[2,] = quantile(data2$obs.X1[data2$obs.Delta != 1] , c(0.25, 0.50, 0.75))
data2.tab.age[3,] = quantile(data2$obs.X1[data2$obs.U == 1] , c(0.25, 0.50, 0.75))
data2.tab.age[4,] = quantile(data2$obs.X1[data2$obs.U == 0] , c(0.25, 0.50, 0.75))
print(xtable(formatC(data2.tab.age, format = "f", digits = 2)))

## distribution of age stratified (a) parity and by (b) BMI
pdf("Figure/null_hist.pdf", width = 12, height = 8)
par(mfrow = c(1,3),  cex.lab = 2, oma = c(5, 5, 2, 2),
    cex.main = 2, cex.axis = 1.5, tcl = 0.5, lwd=2,
    mai = c(0.7, 0.3, 0.3, 0.3), oma = c(5, 5, 7, 2))

hist(nullwithout$obs.X2[nullwithout$obs.Delta == 1 & nullwithout$obs.bmi1 == 1],
     col = rgb(0,0,0,0.8),
     main = "SVD",
     xlab = "", ylim = c(0, 400), xlim = c(13, 45))
hist(nullwithout$obs.X2[nullwithout$obs.Delta == 1 & nullwithout$obs.bmi2 == 1],
     col = rgb(0,0,0,0.5), density=10, angle= 90, 
     add = TRUE)
hist(nullwithout$obs.X2[nullwithout$obs.Delta == 1 & nullwithout$obs.bmi2 == 1],
     col = rgb(0,0,0,0.5), density=10, angle= 180, 
     add = TRUE)
hist(nullwithout$obs.X2[nullwithout$obs.Delta == 1 & nullwithout$obs.bmi3 == 1],
     col = rgb(0,0,0,0.3), add = TRUE)

hist(nullwithout$obs.X2[nullwithout$obs.Delta == 2 & nullwithout$obs.bmi1 == 1],
     col = rgb(0,0,0,0.8),
     main = "CS", ylab= "",
     xlab = "", ylim = c(0, 15), xlim = c(13, 45),
     breaks = seq(13, 45, 2))
hist(nullwithout$obs.X2[nullwithout$obs.Delta == 2 & nullwithout$obs.bmi2 == 1],
     col = rgb(0,0,0,0.5), density=10, angle= 90, 
     add = TRUE, breaks = seq(13, 45, 2))
hist(nullwithout$obs.X2[nullwithout$obs.Delta == 2 & nullwithout$obs.bmi2 == 1],
     col = rgb(0,0,0,0.5), density=10, angle= 180, 
     add = TRUE, breaks = seq(13, 45, 2))
hist(nullwithout$obs.X2[nullwithout$obs.Delta == 2 & nullwithout$obs.bmi3 == 1],
     col = rgb(0,0,0,0.3), add = TRUE, breaks = seq(13, 45, 2))


hist(nullwithout$obs.X2[nullwithout$obs.Delta == 3 & nullwithout$obs.bmi1 == 1],
     col = rgb(0,0,0,0.8),
     main = "OVD", ylab= "",
     xlab = "", ylim = c(0, 40), xlim = c(13, 45),
     breaks = seq(13, 45, 2))
hist(nullwithout$obs.X2[nullwithout$obs.Delta == 3 & nullwithout$obs.bmi2 == 1],
     col = rgb(0,0,0,0.5), density=10, angle= 90, 
     add = TRUE, breaks = seq(13, 45, 2))
hist(nullwithout$obs.X2[nullwithout$obs.Delta == 3 & nullwithout$obs.bmi2 == 1],
     col = rgb(0,0,0,0.5), density=10, angle= 180, 
     add = TRUE, breaks = seq(13, 45, 2))
hist(nullwithout$obs.X2[nullwithout$obs.Delta == 3 & nullwithout$obs.bmi3 == 1],
     col = rgb(0,0,0,0.3), add = TRUE,breaks = seq(13, 45, 2))
legend("topright", c("Underweight/Normal", "Overweight", "Obese"), 
       col = c(rgb(0,0,0,0.8), rgb(0,0,0,0.5), rgb(0,0,0,0.3)),
       pch = c(15, 12 , 15), cex=2, bty='n')
mtext("Nulliparous women without epidural", side = 1, line = -53, outer = TRUE, cex= 2)
mtext("Mother's age", side = 1, line = 1, outer = TRUE, cex= 2)
mtext("Frequency", side = 2, line = 1, outer = TRUE, cex= 2)
dev.off()



pdf("Figure/multi_hist.pdf", width = 12, height = 8)
par(mfrow = c(1,2),  cex.lab = 2, oma = c(5, 5, 2, 2),
    cex.main = 2, cex.axis = 1.5, tcl = 0.5, lwd=2,
    mai = c(0.7, 0.3, 0.3, 0.3), oma = c(5, 5, 7, 2))

hist(multiwithout$obs.X2[multiwithout$obs.Delta == 1 & multiwithout$obs.bmi1 == 1],
     col = rgb(0,0,0,0.8),
     main = "SVD",
     xlab = "", ylim = c(0, 800), xlim = c(14, 53))
hist(multiwithout$obs.X2[multiwithout$obs.Delta == 1 & multiwithout$obs.bmi2 == 1],
     col = rgb(0,0,0,0.5), density=10, angle= 90, 
     add = TRUE)
hist(multiwithout$obs.X2[multiwithout$obs.Delta == 1 & multiwithout$obs.bmi2 == 1],
     col = rgb(0,0,0,0.5), density=10, angle= 180, 
     add = TRUE)
hist(multiwithout$obs.X2[multiwithout$obs.Delta == 1 & multiwithout$obs.bmi3 == 1],
     col = rgb(0,0,0,0.3), add = TRUE)

hist(multiwithout$obs.X2[multiwithout$obs.Delta != 1 & multiwithout$obs.bmi1 == 1],
     col = rgb(0,0,0,0.8),
     main = "non-SVD", ylab ="",
     xlab = "", ylim = c(0, 30), xlim = c(14, 53),
     breaks = seq(14, 53, 2))
hist(multiwithout$obs.X2[multiwithout$obs.Delta != 1 & multiwithout$obs.bmi2 == 1],
     col = rgb(0,0,0,0.5), density=10, angle= 90, 
     add = TRUE, breaks = seq(14, 53, 2))
hist(multiwithout$obs.X2[multiwithout$obs.Delta != 1 & multiwithout$obs.bmi2 == 1],
     col = rgb(0,0,0,0.5), density=10, angle= 180, 
     add = TRUE, breaks = seq(14, 53, 2))
hist(multiwithout$obs.X2[multiwithout$obs.Delta != 1 & multiwithout$obs.bmi3 == 1],
     col = rgb(0,0,0,0.3), add = TRUE, breaks = seq(14, 53, 2))
legend("topright", c("Underweight/Normal", "Overweight", "Obese"), 
       col = c(rgb(0,0,0,0.8), rgb(0,0,0,0.5), rgb(0,0,0,0.3)),
       pch = c(15, 12 , 15), cex=1.5, bty = 'n')
mtext("Multiparous women without epidural", side = 1, line = -32, outer = TRUE, cex= 2)
mtext("Mother's age", side = 1, line = 1, outer = TRUE, cex= 2)
mtext("Frequency", side = 2, line = 1.2, outer = TRUE, cex= 2)
dev.off()


pdf("Figure/multi_hist_beyond10.pdf", width = 12, height = 8)
par(mfrow = c(1,2),  cex.lab = 2, oma = c(5, 5, 2, 2),
    cex.main = 2, cex.axis = 1.5, tcl = 0.5, lwd=2,
    mai = c(0.7, 0.3, 0.3, 0.3), oma = c(5, 5, 7, 2))

hist(data2$obs.X2[data2$obs.Delta == 1 & data2$obs.bmi1 == 1],
     col = rgb(0,0,0,0.8),
     main = "SVD",
     xlab = "", ylim = c(0, 400), xlim = c(14, 53))
hist(data2$obs.X2[data2$obs.Delta == 1 & data2$obs.bmi2 == 1],
     col = rgb(0,0,0,0.5), density=10, angle= 90, 
     add = TRUE)
hist(data2$obs.X2[data2$obs.Delta == 1 & data2$obs.bmi2 == 1],
     col = rgb(0,0,0,0.5), density=10, angle= 180, 
     add = TRUE)
hist(data2$obs.X2[data2$obs.Delta == 1 & data2$obs.bmi3 == 1],
     col = rgb(0,0,0,0.3), add = TRUE)

hist(data2$obs.X2[data2$obs.Delta != 1 & data2$obs.bmi1 == 1],
     col = rgb(0,0,0,0.8),
     main = "CS or OVD", ylab ="",
     xlab = "", ylim = c(0, 20), xlim = c(14, 53),
     breaks = seq(14, 53, 2))
hist(data2$obs.X2[data2$obs.Delta != 1 & data2$obs.bmi2 == 1],
     col = rgb(0,0,0,0.5), density=10, angle= 90, 
     add = TRUE, breaks = seq(14, 53, 2))
hist(data2$obs.X2[data2$obs.Delta != 1 & data2$obs.bmi2 == 1],
     col = rgb(0,0,0,0.5), density=10, angle= 180, 
     add = TRUE, breaks = seq(14, 53, 2))
hist(data2$obs.X2[data2$obs.Delta != 1 & data2$obs.bmi3 == 1],
     col = rgb(0,0,0,0.3), add = TRUE, breaks = seq(14, 53, 2))
legend("topright", c("Underweight/Normal", "Overweight", "Obese"), 
       col = c(rgb(0,0,0,0.8), rgb(0,0,0,0.5), rgb(0,0,0,0.3)),
       pch = c(15, 12 , 15), cex=1.5, bty = 'n')
mtext( expression(paste("Multiparous women without epidural ( delivery ", time >= 10, " min)" )), 
       side = 1, line = -32, outer = TRUE, cex= 2)
mtext("Mother's age", side = 1, line = 1, outer = TRUE, cex= 2)
mtext("Frequency", side = 2, line = 1.2, outer = TRUE, cex= 2)
dev.off()