options(error=recover, warn=0)
library("tilingArray")
source("scripts/readSegments.R")

n = length(polyA2$"4.+.dat"$x)
J = polyA2$"4.+.seg"$J

## LL = n/2(1-log(2pi)) - n*log(sd)
## J really is -log(estimated standard deviation)

LL = n/2*(1-log(2*pi)+2*J)
x  = 2*(seq(along=LL))
## LL = n/2*log(2*pi*exp(-2*J))

par(mfrow=c(2,1))
plot(x, LL, type="l", xlab="model complexity (n)", ylab="log likelihood")
plot(x[-1], diff(LL)/diff(x), type="l", ylim=c(0,50), xlab="model complexity (n)", ylab="d(log likelihood)/dn")

cols = c("blue", "red")
abline(h=c(2, 2*log(n)), col=cols)
