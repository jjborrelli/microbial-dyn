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

#B1 <- matrix(runif(11*1000), ncol = 1000, nrow = 11)
B1 <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/intialABUND2.csv", row.names = 1)
B2 <- B1
B2[9,] <- 0
strt <- Sys.time()
out2 <- list()
out3 <- list()
tte2 <- matrix(nrow = 1000, ncol = 11)
tte3 <- matrix(nrow = 1000, ncol = 11)
for(i in 1:1000){
  out2[[i]] <- ode(B1[,i], 1:1000, parms = parms, func = lvmod2, events = list(func = ext1, time = 1:1000))
  out3[[i]] <- ode(B2[,i], 1:1000, parms = parms, func = lvmod2, events = list(func = ext1, time = 1:1000))
  tte2[i,] <- apply(out2[[i]][,-1], 2, function(x) 1000 - sum(x == 0))
  tte3[i,] <- apply(out3[[i]][,-1], 2, function(x) 1000 - sum(x == 0))
  
  print(i)
}
end <- Sys.time()
end - strt

eq2 <- t(sapply(out2, function(x) tail(x, 1)[-1]))
unique(apply(eq2, 1, function(x) which(x > 0)))
barplot(t(eq2[999:1000,-1]))
cbind(eq2[1000,][order(eq2[1000,], decreasing = T)],colnames(parms$m)[order(eq2[1000,], decreasing = T)])

eq3 <- t(sapply(out3, function(x) tail(x, 1)[-1]))
unique(apply(eq3, 1, function(x) which(x > 0)))
barplot(t(eq3[999:1000,-1]))
cbind(eq3[1000,][order(eq3[1000,], decreasing = T)],colnames(parms$m)[order(eq3[1000,], decreasing = T)])
