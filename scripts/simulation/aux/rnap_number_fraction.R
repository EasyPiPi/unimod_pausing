library(tidyverse)
library(truncnorm)

# get parameters for calculating phi
calculate_f_and_kmean <- function(s, k) {
  # sd is set as 25 here
  x <- round(rnorm(1e7, mean = k, sd = 25))
  x <- x[x >= 17 & x <= 200]
  f0 <- mean(x <= s) 
  f1 <- mean((x > s) & (x <= 2 * s))
  f2 <- mean((x > 2 * s) & (x <= 3 * s))
  f3 <- mean(x > 3 * s)
  return(c("f0" = f0, "f1" = f1, "f2" = f2, "f3" = f3))
}

set.seed(20221005)
x <- map2(c("s33" = 33, "s50" = 50, "s70" = 70), 50, calculate_f_and_kmean)

calculate_f_and_kmean_2 <- function(s, k) {
  kmin <- 17
  kmax <- 200
  sd <- 25
  
  p0 <- ptruncnorm(s, a = kmin - 1, b = kmax, mean = k, sd = sd)
  p1 <- ptruncnorm(2 * s, a = kmin - 1, b = kmax, mean = k, sd = sd)
  p2 <- ptruncnorm(3 * s, a = kmin - 1, b = kmax, mean = k, sd = sd)
  
  return(c("f0" = p0, "f1" = p1 - p0, "f2" = p2 - p1, "f3" = 1 - p2))
}

y <- map2(c("s33" = 33, "s50" = 50, "s70" = 70), 50, calculate_f_and_kmean_2)

kmin <- 17
kmax <- 200
sd <- 25
k <- 50
s <- 50

f0 <- integrate(dnorm, lower = kmin, upper = s, mean = k, sd = sd)
f1 <- integrate(dnorm, lower = s, upper = 2 * s, mean = k, sd = sd)
f2 <- integrate(dnorm, lower = 2 * s, upper = 3 * s, mean = k, sd = sd)
f3 <- integrate(dnorm, lower = 3 * s, upper = kmax, mean = k, sd = sd)

total <- f0$value + f1$value + f2$value + f3$value
f0 <- f0$value / total
f1 <- f1$value / total
f2 <- f2$value / total
f3 <- f3$value / total

