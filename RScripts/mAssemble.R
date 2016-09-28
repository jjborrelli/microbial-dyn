library(dplyr)
library(reshape2)
library(ggplot2)
library(deSolve)
library(rootSolve)


ints1 <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/ecomod-ints.csv", row.names = 1)
grow1 <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/ecomod-Growth.csv")

parms <- list(alpha = unlist(grow1), m = as.matrix(ints1))

parms.m <- parms
parms.m$alpha <- parms.m$alpha[1:2]
parms.m$m <- parms.m$m[1:2, 1:2]
m2 <- parms$m


comm <- ode(runif(2), 1:20, parms = parms.m, func = lvmod2, events = list(func = ext1, time = 1:20))
matplot(comm[,-1], typ = "l")

parms.m$m <- rbind(cbind(parms.m$m, rnorm(nrow(parms.m$m), mean(m2), sd(m2))), rnorm(ncol(parms.m$m)+1, mean(m2), sd(m2)))
diag(parms.m$m) <- c(diag(parms.m$m)[-length(diag(parms.m$m))], -abs(rnorm(1, .77, 1.23)))
parms.m$alpha <- c(parms.m$alpha, rbeta(1, 1.5, 2))

comm <- ode(c(comm[20,-1], .1), 1:20, parms = parms.m, func = lvmod2, events = list(func = ext1, time = 1:20))
matplot(comm[,-1], typ = "l")
