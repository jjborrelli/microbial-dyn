library(dplyr)
library(reshape2)
library(ggplot2)
library(deSolve)
library(rootSolve)

stOTU <- as.data.frame(t(read.csv("Data/stein-OTUcount.csv", header = F, row.names = 1)))
head(stOTU)
table(stOTU[stOTU$Population == 2,2])



ints1 <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/ecomod-ints.csv", row.names = 1)
grow1 <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/ecomod-Growth.csv")

parms <- list(alpha = unlist(grow1), m = as.matrix(ints1))

matplot(select(stOTU, undefined_genus_of_Enterobacteriaceae:Other), typ = "l")

otuCOUNT <- select(stOTU, undefined_genus_of_Enterobacteriaceae:Other)
barplot(t(otuCOUNT[stOTU$Population == 1 & stOTU$Replicate == 1,]))  



stOTU[stOTU$Population == 1 & stOTU$Replicate == 1,]

res1 <- ode(unlist(otuCOUNT[stOTU$Population == 1 & stOTU$Replicate == 2,][1,]), 1:100, parms = parms, func = lvmod2, events = list(func = ext1, time =  1:100))

times <- stOTU[stOTU$Population == 1 & stOTU$Replicate == 1,4]
barplot(t(res1[res1[,1] %in% times, -1]))

times2 <- rep(times, 11)

ggplot(melt(res1[,-1]), aes(x = times2, y = value, fill = Var2)) + geom_area(position = "stack")


table(stOTU$Population)
table(stOTU$Replicate)
combos <- unique(cbind(stOTU$Population, stOTU$Replicate))
par(mfrow = c(3,3))
for(i in 1:9){
  barplot(t(otuCOUNT[stOTU$Population == combos[i,1] & stOTU$Replicate == combos[i,2],]), main = paste(combos[i,], sep = ","))
}

res1 <- list()
for(i in 1:9){
  res1[[i]] <- ode(unlist(otuCOUNT[stOTU$Population == combos[i,1] & stOTU$Replicate == combos[i,2],][1,]), 0:30, parms = parms, func = lvmod2, events = list(func = ext1, time =  1:30))
  matplot(res1[[i]][,-1], typ = "l", lwd = 2)
}

res2 <- lapply(res1, function(x) x[,-1])
res3 <- lapply(res2, function(x) t(apply(x, 1, function(y) y/sum(y))))
ggplot(melt(res3), aes(x = Var1, y = value, fill = Var2)) + geom_area(position = "stack") + facet_wrap(~L1)

##############################################################################################
##############################################################################################
##############################################################################################
#####   Get equilibrium/steady-state communities
#####   Determine their eigenvalues
#####      1. with C. diff 
#####      2. without C. diff


#B1 <- matrix(runif(11*1000), ncol = 1000, nrow = 11)
B1 <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/intialABUND2.csv", row.names = 1)
B3 <- B2 <- B1
B2[9,] <- runif(1000, 10^-9, 10^-6)
B3[9,] <- 0
strt <- Sys.time()
out1 <- list()
out2 <- list()
out3 <- list()
tte2 <- matrix(nrow = 1000, ncol = 11)
tte3 <- matrix(nrow = 1000, ncol = 11)
for(i in 1:1000){
  out1[[i]] <- ode(B1[,i], 1:1000, parms = parms, func = lvmod2, events = list(func = ext1, time = 1:1000))
  out2[[i]] <- ode(B2[,i], 1:1000, parms = parms, func = lvmod2, events = list(func = ext1, time = 1:1000))
  out3[[i]] <- ode(B3[,i], 1:1000, parms = parms, func = lvmod2, events = list(func = ext1, time = 1:1000))
  
  print(i)
}
end <- Sys.time()
end - strt

eq1 <- t(sapply(out1, function(x) tail(x, 1)[-1]))
unique(apply(eq1, 1, function(x) which(x > 0)))
barplot(t(eq1[999:1000,-1]))
cbind(eq1[1000,][order(eq1[1000,], decreasing = T)],colnames(parms$m)[order(eq1[1000,], decreasing = T)])
unique(apply(eq1, 1, function(x) max(Re(eigen(jacobian.full(x, lvmod2, parms = parms))$values))))
eqabund <- colMeans(eq1[which(eq2[,9] != 0),])
max(eigen(jacobian.full(eqabund[which(eqabund > 0)], lvmod2, 
                        parms = list(alpha = parms$alpha[which(eqabund > 0)], 
                                     m = parms$m[which(eqabund > 0), which(eqabund > 0)])))$values)
colMeans(eq2[which(eq2[,9] != 0),])

eq2 <- t(sapply(out2, function(x) tail(x, 1)[-1]))
unique(apply(eq2, 1, function(x) which(x > 0)))
barplot(t(eq2[999:1000,-1]))
cbind(eq2[1000,][order(eq2[1000,], decreasing = T)],colnames(parms$m)[order(eq2[1000,], decreasing = T)])
unique(apply(eq2, 1, function(x) max(Re(eigen(jacobian.full(x, lvmod2, parms = parms))$values))))
eqabund <- colMeans(eq2)
max(Re(eigen(jacobian.full(eqabund[which(eqabund > 0)], lvmod2, 
                        parms = list(alpha = parms$alpha[which(eqabund > 0)], 
                                     m = parms$m[which(eqabund > 0), which(eqabund > 0)])))$values))
colMeans(eq2)

eq3 <- t(sapply(out3, function(x) tail(x, 1)[-1]))
unique(apply(eq3, 1, function(x) which(x > 0)))
barplot(t(eq3[999:1000,-1]))
cbind(eq3[1000,][order(eq3[1000,], decreasing = T)],colnames(parms$m)[order(eq3[1000,], decreasing = T)])
eqabund2 <- colMeans(eq3)
max(Re(eigen(jacobian.full(eqabund2[which(eqabund2 > 0)], lvmod2, 
                           parms = list(alpha = parms$alpha[which(eqabund2 > 0)], 
                                        m = parms$m[which(eqabund2 > 0),which(eqabund2 > 0)])))$values))


barplot(cbind(eq1[1000,], eq2[1000,], eq3[1000,]))


##############################################################################################
##############################################################################################
##############################################################################################
#####   Perturb equilibrium/steady-state communities
#####   Determine the effect on their eigenvalues
#####      1. with C. diff 
#####      2. without C. diff


## Comm 1 with Cdiff
init.abund <- eqabund[which(eqabund !=0)]
growth <- parms$alpha[which(eqabund !=0)]
imat <- parms$m[which(eqabund !=0),which(eqabund !=0)] 

par1 <- list(alpha = growth, m = imat)
out <- ode(init.abund, 1:1000, parms = par1, func = lvmod2, events = list(func = ext1, time = 1:1000))
matplot(out[,-1], typ = "l")

ia2 <- init.abund
for(i in 1:length(init.abund)){
  ia2[i] <- 0#ia2[i]*(rbeta(1,1,4)*sample(c(1,-1),1,prob = c(.5,.5)))
  out <- ode(ia2, 1:1000, parms = par1, func = lvmod2, events = list(func = ext1, time = 1:1000))
}

imat2 <- imat*-1
diag(imat2) <- diag(imat)
par1 <- list(alpha = growth, m = imat2)
eig <- eigen(jacobian.full(init.abund, lvmod2, parms = par1))
plot(Re(eig$values), Im(eig$values))


e1 <- c()
for(i in 1:length(imat)){
  imat2 <- imat
  imat2[i] <- 0
  par1 <- list(alpha = growth, m = imat2)
  eig <- eigen(jacobian.full(init.abund, lvmod2, parms = par1))
  e1[i] <- max(Re(eig$values))
}

e2 <- c()
cor <- c()
eg1 <- t(combn(1:nrow(imat),2))
for(i in 1:nrow(eg1)){
  imat2 <- imat
  imat2[eg1[i,1], eg1[i,2]] <- 0
  imat2[eg1[i,2], eg1[i,1]] <- 0
  par1 <- list(alpha = growth, m = imat2)
  
  cor[i] <- cor.test(imat2[upper.tri(imat2)], t(imat2)[upper.tri(imat2)])$estimate
  eig <- eigen(jacobian.full(init.abund, lvmod2, parms = par1))
  e2[i] <- max(Re(eig$values))
}

hist(e2)
eg1[which(e2 > 0),]
imat[eg1[,1]]

fun1 <- function(x, mat){
  c(mat[x[1], x[2]],mat[x[2], x[1]])
}

t(apply(eg1, 1, fun1, mat = imat))[which(e2 > 0),]
plot(cor, e2)

e3 <- c()
for(i in 1:1000){
  #imat2 <- matrix(rnorm(length(imat),mean = imat, sd = .01), 10, 10)
  imat2 <- imat
  imat2[unique(eg1[which(e2 > 0),1]),unique(eg1[which(e2 > 0),2])] <- rnorm(length(imat2[unique(eg1[which(e2 > 0),1]),unique(eg1[which(e2 > 0),2])]), imat2[unique(eg1[which(e2 > 0),1]),unique(eg1[which(e2 > 0),2])], .0001)
  imat2[unique(eg1[which(!e2 > 0),1]),unique(eg1[which(!e2 > 0),2])] <- rnorm(length(imat2[unique(eg1[which(!e2 > 0),1]),unique(eg1[which(!e2 > 0),2])]), imat2[unique(eg1[which(!e2 > 0),1]),unique(eg1[which(!e2 > 0),2])], .15)
  par1 <- list(alpha = growth, m = imat2)
  eig <- eigen(jacobian.full(init.abund, lvmod2, parms = par1))
  e3[i] <- max(Re(eig$values))
}

hist(e3)


##############################################################################################
##############################################################################################
##############################################################################################
#####   Look at strong/core microbiome + keep it constant
#####   Vary weakly interacting microbe ints
#####      1 
#####      2

core <- order(eq2[1000,], decreasing = T)[1:5]
periph <- order(eq2[1000,], decreasing = T)[6:11]

par <- parms

## Periphery
par$m[periph, periph] <- rnorm(length(parms$m[periph, periph]), parms$m[periph, periph], abs(parms$m[periph, periph]*.25)) 

par.lst <- lapply(1:1000, function(x){par <- parms; par$m[periph, periph] <- rnorm(length(parms$m[periph, periph]), parms$m[periph, periph], abs(parms$m[periph, periph]));diag(par$m) <- diag(parms$m);return(par)})


strt <- Sys.time()
res1 <- list()
for(i in 1:1000){
  res1[[i]] <- ode(B1[,1], 1:1000, parms = par.lst[[i]], func = lvmod2, events = list(func = ext1, time = 1:1000))
  
  print(i)
}
end <- Sys.time()
end - strt

eqcp <- t(sapply(res1, function(x) tail(x, 1)[-1]))
table(apply(eqcp, 1, function(x) paste0(which(x > 0), collapse = ",")))

eqcomm <- lapply(unique(apply(eqcp, 1, function(x) which(x > 0))), paste0, collapse = ",")
eqc <- apply(eqcp, 1, function(x) paste0(which(x > 0), collapse = ","))

barplot(t(sapply(1:11, function(x) if(sum(eqc %in% eqcomm[[x]]) > 1){colMeans(eqcp[eqc %in% eqcomm[[x]],])}else{eqcp[eqc %in% eqcomm[[x]],]}))) 
mabund <- colMeans(eqcp[eqc %in% eqcomm[[1]],])
plist <- par.lst[eqc %in% eqcomm[[1]]]

j1 <- list()
for(i in 1:length(plist)){
  j1[[i]] <- jacobian.full(mabund, lvmod2, parms = plist[[i]])
}


## Core

par.lst2 <- lapply(1:1000, function(x){par <- parms; par$m[core, core] <- rnorm(length(parms$m[core, core]), parms$m[core, core], abs(parms$m[core, core]*.25));diag(par$m) <- diag(parms$m);return(par)})


strt <- Sys.time()
res1 <- list()
for(i in 1:1000){
  res1[[i]] <- ode(B1[,1], 1:1000, parms = par.lst2[[i]], func = lvmod2, events = list(func = ext1, time = 1:1000))
  
  print(i)
}
end <- Sys.time()
end - strt

eqCORE <- t(sapply(res1, function(x) tail(x, 1)[-1]))
table(apply(eqCORE, 1, function(x) paste0(which(x > 0), collapse = ",")))

eqcommC <- lapply(unique(apply(eqCORE, 1, function(x) which(x > 0))), paste0, collapse = ",")
eqco <- apply(eqCORE, 1, function(x) paste0(which(x > 0), collapse = ","))

barplot(sapply(1:57, function(x) if(sum(eqco %in% eqcommC[[x]]) > 1){colMeans(eqCORE[eqco %in% eqcommC[[x]],])}else{eqCORE[eqco %in% eqcommC[[x]],]})[,-52]) 



######
######

mco <- parms$m[core, core]
mpe <- parms$m[periph, periph]

c(mean(mco[mco < 0]), mean(mco[mco > 0]))
c(mean(mpe[mpe < 0]), mean(mpe[mpe > 0]))

cin <- cbind(mco[upper.tri(mco)], t(mco)[upper.tri(mco)]) # 4 mut, 4 pred, 2 comp
pin <- cbind(mpe[upper.tri(mpe)], t(mpe)[upper.tri(mpe)]) # 3 mut, 6 pred, 6 comp

sum(cin > 0)/length(cin)
sum(pin > 0)/length(pin)
######
######


id <- 1:11
id[core] <- "green"
id[periph] <- "blue"

plot(parms$alpha, diag(parms$m), col = id, pch = 16)
