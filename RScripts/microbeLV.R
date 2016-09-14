library(deSolve)
library(rootSolve)


mouse <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/murineMICROBEint.csv", row.names = 1)
mouse2 <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/murineMICROBEint2.csv", row.names = 1)
mouseSE <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/murineMICROBEintSE.csv", row.names = 1)
alphas <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/murineMICROBEa.csv")

mse <- as.matrix(mouseSE)
mse[is.na(mse)] <- 0

m <- as.matrix(mouse)
mALT <- as.matrix(mouse2)
mSTR <- as.matrix(mouse2)
mSTR[abs(mSTR) < 2] <- 0

lvmod <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha * (1 - state/1) + state * parms$m %*% state
    
    list(dB)
  })
}

lvmod2 <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha + state * parms$m %*% state
    
    list(dB)
  })
}

ext1 <- function (times, states, parms){
  with(as.list(states), {
    states[states < 10^-30] <- 0 
    return(c(states))
  })
}


parms <- list(alpha = (alphas$Alpha), m = mALT)
parms <- list(alpha = (alphas$Alpha), m = m)
parms <- list(alpha = (alphas$Alpha), m = mSTR)

system.time(
res2 <- ode(res1[1000,-1], 1:1000, parms = parms, func = lvmod, events = list(func = ext1, time =  1:1000))
)
#res1
matplot(res2[,-1], typ = "l", lwd = 2, ylab = "B", xlab = "Time")


B <- matrix(runif(17 * 1000, .5, 1), nrow = 17, ncol = 1000)

strt <- Sys.time()
out <- list()#matrix(nrow = 1000, ncol = 17)
tte <- matrix(nrow = 1000, ncol = 17)
for(i in 1:1000){
  res <- ode(B[,i], 1:1000, parms = parms, func = lvmod, events = list(func = ext1, time = 1:1000))
  tte[i,] <- apply(res[,-1], 2, function(x) 1000 - sum(x == 0))
  out[[i]] <- res#[1000, -1]
  print(i)
}
end <- Sys.time()
end - strt

n.uniq <- unique(lapply(out, function(x) which(x[1000,-1] > 0)))
size3 <- sapply(lapply(out, function(x) which(x[1000,-1] > 0)), function(y) length(y))
eqcomm <- t(sapply(out, function(x){x[1000,-1]}))
boxplot(eqcomm[size3 == 2,])
boxplot(eqcomm[size3 == 3,])
boxplot(eqcomm[size3 == 5,])


resi <- list()
pert <- c()
for(i in 1:1000){
  p1 <- out[[1]][1000,-1]
  neg <- FALSE
  while(!neg){
    p1[17] <- p1[sample(1:17, 1)] + abs(rnorm(1, 0, .2))
    neg <- sum(p1 < 0) == 0
  }
  pert[i] <- p1[17]
  resi[[i]] <- ode(p1, 1:1000, parms = parms, func = lvmod, events = list(func = ext1, time = 1:1000))
  print(i)
}
plot(round(rowSums(t(sapply(resi, function(x) x[1000,-1]))), 5))

plot(sapply(resi, function(x) sum(x[,18] > 0))~pert)

## TO DO
##  - FINISH WORKING OUT ABOVE METHOD TO TEST PERTURBATIONS OF EQUILIBRIUM COMMUNITIES




persist <- apply(out, 1, function(x){sum(x > 0)})
hist(persist)
p2 <- persist/17
fit1 <- glm(p2~t(B), family = "quasibinomial")
summary(fit1)
plot(fit1$coefficients[-1], apply(out, 2, median))


m2 <- m
m2[which(abs(m) <= 1)] <- 0
parms2 <- list(alpha = (alphas$Alpha)[-1], m = m[-1,-1])

out2 <- matrix(nrow = 100, ncol = 16)
for(i in 1:100){
  out2[i,] <- ode(B[-1,i], 1:1000, parms = parms2, func = lvmod, events = list(func = ext1, time = 1:1000))[1000, -1]
  print(i)
}

apply(out2, 1, function(x){sum(x > 0)})
plot(apply(out2, 1, function(x){sum(x > 0)})~persist[1:100])

dropoutSP <- function(x, iter = 100, t = 1000, B, par, mod = lvmod, ev = ext1){
  par$alpha <- par$alpha[-x]
  par$m <- par$m[-x, -x]
  
  out2 <- matrix(nrow = iter, ncol = nrow(par$m))
  for(j in 1:iter){
    out2[j,] <- ode(B[-x,j], 1:t, parms = par, func = mod, events = list(func = ev, time = 1:t))[t, -1]
  }
  return(out2)
}


dropoutSP2 <- function(x, iter = 100, t = 1000, B, par, mod = lvmod, ev = ext1){
  B[x,] <- 0
  out2 <- list()#matrix(nrow = iter, ncol = nrow(par$m))
  for(j in 1:iter){
    out2[[j]] <- ode(B[,j], 1:t, parms = par, func = mod, events = list(func = ev, time = 1:t))[,-1]
  }
  return(out2)
}


dropoutINT <-  function(x, iter = 100, t = 1000, B, par, mod = lvmod, ev = ext1){
  par$m[m !=0][x] <- 0
  
  out2 <- matrix(nrow = iter, ncol = nrow(par$m))
  for(j in 1:iter){
    out2[j,] <- ode(B[,j], 1:t, parms = par, func = mod, events = list(func = ev, time = 1:t))[t, -1]
    cat(j, " ")
  }
  return(out2)
}


B <- matrix(runif(17 * 1000, .5, 1), nrow = 17, ncol = 1000)
parms <- list(alpha = (alphas$Alpha), m = m)


p1mu <- c()
p1med <- c()
temp <- list()
p1 <- list()
for(i in 1:17){
  temp[[i]] <- dropoutSP(x = i, iter = 500, t = 1000, B = B, par = parms)
  p1[[i]] <- apply(temp[[i]], 1, function(x) sum(x > 0))
  
  p1mu[i] <- mean(p1[[i]])
  p1med[i] <- median(p1[[i]])
  
  print(i)
}


presences <- lapply(lapply(temp, function(x) unlist(apply(x, 1, function(j) which(j != 0)))), table)
1-p1mu/17

nz <- t(sapply(temp2, function(x){apply(x, 2, function(j) sum(j != 0))}))
diag(nz) <- NA
boxplot(nz/200)

strt <- Sys.time()
#p2mu <- c()
#p2med <- c()
temp2 <- lapply(1:17, function(x) lapply(1:500, function(y) matrix(NA, nrow = 1000, ncol = 17)))
#p2 <- list()
for(i in 1:17){
  temp2[[i]] <- dropoutSP2(x = i, iter = 500, t = 1000, B = B, par = parms)
  #p2[[i]] <- apply(temp2[[i]], 1, function(x) sum(x > 0))
  
  #p2mu[i] <- mean(p2[[i]])
  #p2med[i] <- median(p2[[i]])
  
  print(i)
}
end <- Sys.time()
end-strt

length(temp2[[1]])
#get time to extinction
t2e <- lapply(temp2, function(x){t(sapply(x, function(y){apply(y[,-1], 2, function(x) 1000 - sum(x == 0))}))})

loc <- matrix(c((1:17) - .5, rep(-10, 17), (1:17) + .5, rep(1010, 17)), ncol = 4)
par(mfrow = c(6, 3))
for(i in 1:17){
  boxplot(t2e[[i]])
  rect(loc[i,1], loc[i,2], loc[i,3], loc[i,4], col = "grey75")
}

#get equilibrium comms
eq <- lapply(temp2, function(x){t(sapply(x, function(y){tail(y, 1)}))})
eq[[1]]

u2 <- lapply(lapply(eq, function(x){lapply(1:nrow(x), function(y) which(x[y,] > 0))}), unique)
sapply(u2, length)
lapply(u2, function(x) sapply(x, length))
#write.csv(melt(u2), "~/Desktop/GitHub/microbial-dyn/Data/removalCOMM3.csv")



library(rootSolve)
eig <- c()
for(i in 1:500){
  jac <- jacobian.full(B[,i], lvmod, parms = parms)
  eig[i]<- max(Re(eigen(jac)$values))
  
}

hist(eig)
parms$m







# Stochastic

t.max <- 1000
B1 <- matrix(runif(17), nrow = t.max, ncol = 17, byrow = T)

for(i in 2:t.max){
  m.sto <- matrix(rnorm(length(m), m, mse), 17, 17)
  parms <- list(alpha = (alphas$Alpha), m = m.sto)
  
  res1 <- ode(B1[i-1,], 1:2, parms = parms, func = lvmod, events = list(func = ext1, time =  1:2))
  B1[i,] <- res1[2,-1]
}

matplot(B1, typ = "l")
B1[t.max,]


stodyn <- function(t.max, ints = m, intSE = mse, a = alphas$Alpha){
  
  B1 <- matrix(runif(17), nrow = t.max, ncol = 17, byrow = T)
  
  for(i in 2:t.max){
    m.sto <- matrix(rnorm(length(ints), ints, intSE), 17, 17)
    parms <- list(alpha = (a), m = m.sto)
    
    res1 <- ode(B1[i-1,], 1:2, parms = parms, func = lvmod, events = list(func = ext1, time =  1:2))
    B1[i,] <- res1[2,-1]
  }
  return(B1)
}

matplot(stodyn(100), typ = "l")

test <- list()
for(i in 1:100){
  test[[i]] <- stodyn(100)
}

test2 <- lapply(1:17, function(x) sapply(test, function(q) q[,x]))
matplot(sapply(test2, function(x) apply(x, 1, median)))


xi <- runif(17)

mre <- c()
for(i in 1:1000){
  xi <- B[,i]
  mre[i] <- max(Re(eigen(jacobian.full(xi, func = lvmod, parms = parms))$values))
}
plot(mre)

out1 <- ode(xi, 1:1000, parms = parms, func = lvmod, events = list(func = ext1, time = 1:1000))
matplot(out1[,-1], typ = "l")









#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

ints1 <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/ecomod-ints.csv", row.names = 1)
grow1 <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/ecomod-Growth.csv")

parms <- list(alpha = unlist(grow1), m = as.matrix(ints1))
parms$m[abs(parms$m) < quantile(abs(parms$m))[3]] <- 0
diag(parms$m) <- diag(as.matrix(ints1))
parms$m[abs(parms$m) > quantile(abs(parms$m))[4]] <- 0

system.time(
  res1 <- ode(res, 1:13, parms = parms, func = lvmod2, events = list(func = ext1, time =  1:11))
)
res1[1000,-1]
matplot(res1[,-1], typ = "l", lwd = 2)


B1 <- matrix(runif(11*1000), ncol = 1000, nrow = 11)
strt <- Sys.time()
out2 <- list()#matrix(nrow = 1000, ncol = 11)
tte2 <- matrix(nrow = 1000, ncol = 11)
for(i in 1:1000){
  res <- ode(B1[,i], 1:1000, parms = parms, func = lvmod2, events = list(func = ext1, time = 1:1000))
  tte2[i,] <- apply(res[,-1], 2, function(x) 1000 - sum(x == 0))
  out2[[i]] <- res
  print(i)
}
end <- Sys.time()
end - strt

eq <- t(sapply(out2, function(x) tail(x, 1)[-1]))
eq[1,]
unique(apply(eq, 1, function(x) which(x > 0)))


sapply(u2, length)
lapply(u2, function(x) sapply(x, length))
#write.csv(melt(u2), "~/Desktop/GitHub/microbial-dyn/Data/removalCOMM3.csv")


mout <- lapply(out2, function(x){melt(x[,-1])})
mout2 <- lapply(1:1000, function(x) mout[[x]] <- cbind(mout[[x]], iter = x))
library(data.table)
mout3 <- rbindlist(mout2[sample(1:1000, 100)])
#library(ggplot2)
ggplot(mout3, aes(x = Var1, y = value)) + geom_point() + facet_wrap(~Var2, scales = "free_y")

boxplot(tte2)

eq1 <- t(sapply(out2, function(x) x[1000,-1]))
boxplot(eq1)
apply(eq1, 2, median)

B1 <- matrix(runif(11*500), ncol = 500, nrow = 11)
strt <- Sys.time()
dyn <- lapply(1:11, function(x) lapply(1:500, function(y) matrix(NA, nrow = 1000, ncol = 17)))
for(i in 1:11){
  dyn[[i]] <- dropoutSP2(x = i, iter = 500, t = 1000, B = B1, par = parms, mod = lvmod2)
  print(i)
}
end <- Sys.time()
end-strt

length(temp2[[1]])
#get time to extinction
t2e2 <- lapply(dyn, function(x){t(sapply(x, function(y){apply(y, 2, function(x) 1000 - sum(x == 0))}))})

loc <- matrix(c((1:11) - .5, rep(-10, 11), (1:11) + .5, rep(1010, 11)), ncol = 4)
par(mfrow = c(3, 4))
for(i in 1:11){
  boxplot(t2e2[[i]])
  rect(loc[i,1], loc[i,2], loc[i,3], loc[i,4], col = "grey75")
}

#get equilibrium comms
eq <- lapply(dyn, function(y) t(sapply(y, function(x) tail(x, 1))))

u2 <- lapply(lapply(eq, function(x){lapply(1:nrow(x), function(y) which(x[y,] > 0))}), unique)
sapply(u2, length)
melt(lapply(u2, function(x) sapply(x, length)))

u2 <- lapply(lapply(eq, function(x){lapply(1:nrow(x), function(y) which(x[y,] > 0))}), paste0)
lapply(u2, table)

par(mfrow = c(3,4))
for(i in 1:11){barplot(lapply(u2, table)[[i]], main = i)}
#dev.off()

u3 <- lapply(eq, function(x){lapply(1:nrow(x), function(y) x[y,][which(x[y,] > 0)])})




#write.csv(melt(u2), "~/Desktop/GitHub/microbial-dyn/Data/remCOMM1-stein.csv")


sumpos <- c()
sumneg <- c()
for(i in 1:11){
  sumpos[i] <- sum(c(parms$m[i,][parms$m[i,] > 0], parms$m[,i][parms$m[,i] > 0]))
  sumneg[i] <- sum(c(parms$m[i,][parms$m[i,] < 0], parms$m[,i][parms$m[,i] < 0]))
}
sumpos

dim(dyn[[1]][1][1])


m2 <- m
m2[(abs(m) - as.matrix(mouseSE)) < 0] <- 0 

d1 <- cbind(m[upper.tri(m)], t(m)[upper.tri(m)])
d1.2 <- cbind(m2[upper.tri(m2)], t(m2)[upper.tri(m2)])
d1.3 <- cbind(mSTR[upper.tri(mSTR)], t(mSTR)[upper.tri(mSTR)])
d2 <- cbind(as.matrix(ints1)[upper.tri(as.matrix(ints1))], t(as.matrix(ints1))[upper.tri(as.matrix(ints1))])

dx <- d1.2/abs(d1.2)
dx[is.nan(dx)] <- 0
rowSums(dx)

m[1,]
m[,1]
x <- vector(length = 17)


ncomp <- c()
namen <- c()
for(i in 1:17){
  r1 <- ifelse(mALT[i,] < 0, x <- "neg", ifelse(mALT[i,] > 0, x <- "pos", x <- "0"))
  r2 <- ifelse(mALT[,i] < 0, x <- "neg", ifelse(mALT[,i] > 0, x <- "pos", x <- "0"))
  ncomp[i] <- sum(r1 == "neg" & r2 == "neg")
  namen[i] <- sum(r1 == "neg" & r2 == "0" | r2 == "neg" & r1 == "0")
}



lapply(out, function(x) mALT[which(x[1000,-1] > 0), which(x[1000,-1] > 0)])
lapply(out, function(x) which(x[1000,-1] > 0))


rownames(mSTR) <- colnames(mSTR)
rownames(mALT) <- colnames(mALT)
g1 <- melt(mSTR[c(1,3,4,5,8,11,16),c(1,3,4,5,8,11,16)])
g1 <- melt(mALT[c(1,4,6,13,15),c(1,4,6,13,15)])
g2 <- g1[g1$value != 0,]
g3 <- graph.edgelist(as.matrix(g2)[,1:2])
E(g3)$weights <- g2$value
plot(g3, layout = layout.circle, edge.color = factor(g2$value > 0))




test2 <- melt(lapply(dyn, function(x){t(sapply(x, function(y){apply(y, 2, sd)}))})[-4])
ggplot(test2, aes(x = factor(Var2), y = value)) + geom_boxplot(aes(fill = factor(L1))) + facet_wrap(~L1)



library(rootSolve)
library(GGally)
eig <- c()
states <- matrix(runif(11*1000), 11, 1000)
for(i in 1:1000){
  myj <- jacobian.full(y = states[,i], func = lvmod2, parms = parms)
  eig[i] <- max(Re(eigen(myj)$values))
}
hist(eig)
ggpairs(t(states[,which(eig < 0)])) 



#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
# pairwise

com1 <- t(combn(1:11, 5))
par2 <- parms
eigs1 <- c()
for(i in 1:nrow(com1)){
  par2$m <- parms$m[com1[i,],com1[i,]]
  par2$alpha <- parms$alpha[com1[i,]]
  eigs1[i] <- max(Re(eigen(jacobian.full(solve(parms$m[com1[i,], com1[i,]])%*%parms$alpha[com1[i,]], lvmod2, parms = par2))$values))
  print(i)
}
boxplot(eigs1)
