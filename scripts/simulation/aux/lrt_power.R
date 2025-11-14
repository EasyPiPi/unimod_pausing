# test statistical power
dens1=0.025
dens2=0.03
M=9900
y1=rpois(1000,M*dens1)
y2=rpois(1000,M*dens2)
hist(y1)
hist(y2)
T_stats = y1*log(2*y1/(y1+y2)) + y2*log(2*y2/(y1+y2))
hist(T_stats)
x <- summary(T_stats)
pchisq(2 * x, df = 1, ncp = 0, lower.tail = F, log.p = FALSE)