points = seq(0, 5, 0.0)
lognormal05 = dlnorm(points, 0, 0.5)
lognormal10 = dlnorm(points, 0, 1.0)
weibull15 = dweibull(points, shape = 1.5, scale = 1)
weibull20 = dweibull(points, shape = 2.0, scale = 1)
pdf("Figure/frailty.pdf", height = 8, width = 14)
par(mfrow = c(1,1),   mar = c(7,12,5,3),  cex.lab = 2.8, 
    cex.main = 3, cex.axis = 2, tcl = 0.5)
plot(points, lognormal05, xlab = "Frailty (z)",
     ylab = "Distribution of simulated frailty", type = "l", xlim = c(0,5),
     ylim = c(0,1.0), lwd = 4, mgp = c(5,1,0), xaxt = 'n')
axis(1, at = seq(0, 5, 0.5), labels = seq(0, 5, 0.5))
lines(points, lognormal10, lwd = 4, lty = 2)
lines(points, weibull15,  col = "grey", lwd = 4, lty = 1)
lines(points, weibull20,  col = "grey", lwd = 4, lty = 2)
legend("topright", c(expression(paste("Log-normal ", sigma , "=0.5")), 
                     expression(paste("Log-normal ", sigma, "=1.0")),
                     expression(paste("Weibull ", gamma, "=1.5")), 
                     expression(paste("Weibull ", gamma, "=2.0")))
       , lty = c(1,2,1,2),
       col = c("black", "black", "grey", "grey"), bty = 'n', cex = 2, lwd = 3)
dev.off()

points = seq(0, 5, 0.01)
weibull05 = dweibull(points, shape = 0.5, scale = 1)
weibull10 = dweibull(points, shape = 1.0, scale = 1)
weibull15 = dweibull(points, shape = 1.5, scale = 1)
pdf("Figure/weibull.pdf", height = 8, width = 14)
par(mfrow = c(1,1),   mar = c(7,12,5,3),  cex.lab = 2.8, 
    cex.main = 3, cex.axis = 2, tcl = 0.5)
plot(points,  weibull05, xlab = "Frailty (z)",
     ylab = "Distribution of simulated frailty", type = "l", xlim = c(0,5),
     ylim = c(0,1.0), lwd = 4, mgp = c(5,1,0), xaxt = 'n')
axis(1, at = seq(0, 5, 0.5), labels = seq(0, 5, 0.5))
lines(points, weibull10,  col = "grey", lwd = 4, lty = 2)
lines(points, weibull15,  col = "azure4", lwd = 4, lty = 6)
legend("topright", c(expression(paste("Weibull ", gamma, "=0.5")), 
                     expression(paste("Weibull ", gamma, "=1.0")),
                     expression(paste("Weibull ", gamma, "=1.5")))
       , lty = c(1,2,6),
       col = c("black", "grey", "azure4"), bty = 'n', cex = 2, lwd = 3)
dev.off()

##
sigma = 0.5
betas = c(1, -2); alpha = 1.5
phis = c(0.01, 0.005)
tau = 3/2
psi = 0.01
omega = 2.5
intervals = seq(0.1, 160, 0.1);

true.lambda1 = (phis[1])^(tau)*tau*(intervals)^(tau-1)
true.lambda2 = (phis[2])^(tau)*tau*(intervals)^(tau-1)
true.h = (psi^omega)*omega*(intervals)^(omega - 1)

pdf("Figure/baselines.pdf", height = 8, width = 14)
par(mfrow = c(1,1),   mar = c(7,12,5,3),  cex.lab = 2.8, 
    cex.main = 3, cex.axis = 2, tcl = 0.5)
plot(intervals, true.lambda1, xlab = "Time to event (min)",
     ylab = "Cause-specific baseline hazards \nand baseline hazards for morbidity", type = "l", xlim = c(0,160),
     ylim = c(0,0.04), lwd = 4, mgp = c(5,1,0), xaxt = 'n')
axis(1, at = seq(0, 160, 20), labels = seq(0, 160, 20))
lines(intervals, true.lambda2, lwd = 4, lty = 2)
lines(intervals, true.h,  col = "grey", lwd = 4, lty = 5)
legend("topleft", c(expression(lambda[10](t)), expression(lambda[20](t)), 
                    expression(h[0](t))), lty = c(1,2,5),
       col = c("black", "black", "grey"), bty = 'n', cex = 2, lwd = 3)
dev.off()

