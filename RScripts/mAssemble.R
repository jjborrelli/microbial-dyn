library(dplyr)
library(reshape2)
library(ggplot2)
library(deSolve)
library(rootSolve)


ints1 <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/ecomod-ints.csv", row.names = 1)
grow1 <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/ecomod-Growth.csv")

parms <- list(alpha = unlist(grow1), m = as.matrix(ints1))

ti <- 100
iter <- 500
parms.m <- lapply(1:iter, function(x) parms)
parms.m[[1]] <- parms
parms.m[[1]]$alpha <- parms.m[[1]]$alpha[1:2]
parms.m[[1]]$m <- parms.m[[1]]$m[1:2, 1:2]
m2 <- parms$m

comm <- list()
comm[[1]] <- ode(runif(2), 1:ti, parms = parms.m[[1]], func = lvmod2, events = list(func = ext1, time = 1:ti))
matplot(comm[[1]][,-1], typ = "l")

for(i in 319:iter){
  cond <- F
  while(!cond){
    parms.m[[i]]$m <- rbind(cbind(parms.m[[i-1]]$m, rnorm(nrow(parms.m[[i-1]]$m), mean(m2), sd(m2))), rnorm(ncol(parms.m[[i-1]]$m)+1, mean(m2), sd(m2)))
    diag(parms.m[[i]]$m) <- c(diag(parms.m[[i-1]]$m), -abs(rnorm(1, .77, 1.23)))
    parms.m[[i]]$alpha <- c(parms.m[[i-1]]$alpha, rbeta(1, 1.5, 2))
    
    comm[[i]] <- ode(c(comm[[i-1]][ti,-1], .1), 1:ti, parms = parms.m[[i]], func = lvmod2, events = list(func = ext1, time = 1:ti))
    
    cond <- nrow(comm[[i]]) == ti
  }
  matplot(comm[[i]][,-1], typ = "l")
  print(i)
}

comm[[iter]]

comm2 <- lapply(1:iter, function(x){if(x != iter){for(i in (ncol(comm[[x]]) + 1):(iter+2)){comm[[x]] <- cbind(comm[[x]], rep(0, ti))};return(comm[[x]])}else{return(comm[[x]])}}) 

alltraj <- do.call(rbind, comm2)[,-1]
matplot(alltraj, typ = "l")

reltraj <- t(apply(alltraj, 1, function(x) x/sum(x)))
matplot(reltraj, typ = "l")
barplot(t(reltraj))
colnames(reltraj) <- 1:51
ggplot(melt(reltraj), aes(x = Var1, y = value)) + geom_line(aes(col = factor(Var2)))

eqtraj <- sapply(comm2, function(x) x[ti, -1])
barplot(eqtraj)
releq <- apply(eqtraj, 2, function(x) x/sum(x))
barplot(releq)
rownames(releq) <- 1:51
dim(releq)
ggplot(melt(releq), aes(x = Var2, y = value, fill = factor(Var1))) + geom_bar(stat = "identity") + scale_fill_brewer()

numSP <- apply(releq, 2, function(x) sum(x > 0))
plot(numSP, typ = "o")


eqcomm <- apply(releq, 2, function(x) which(x > 0))
dim(parms.m[[iter]]$m)

ity <- lapply(1:500, function(x){itypes(parms.m[[iter]]$m[eqcomm[[x]],eqcomm[[x]]])})
inty <- do.call(rbind, ity)[,1:3]
inty2 <- t(apply(inty, 1, function(x) x/sum(x)))

d <- sapply(1:500, function(x){mean(diag(parms.m[[iter]]$m[eqcomm[[x]],eqcomm[[x]]]))})
gr <- sapply(1:500, function(x){mean(parms.m[[iter]]$alpha[eqcomm[[x]]])})

plot(inty2[,1], numSP)
plot(inty2[,2], numSP)
plot(inty2[,3], numSP)

summary(lm(numSP~inty2[,1:2]+d+gr))
summary(lm(numSP~inty2[,1:2]))
