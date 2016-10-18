library(igraph)
library(NetIndices)
library(deSolve)

lvmod <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha + state * parms$m %*% state 
    
    list(dB)
  })
}

ext1 <- function (times, states, parms){
  with(as.list(states), {
    states[states < 10^-30] <- 0 
    if(sum(states >= 100) >= 1){states<-rep(0, length(states))} 
    return(c(states))
  })
}

icor <- function(x){
  i1 <- x[upper.tri(x)]
  i2 <- t(x)[upper.tri(x)] 
  
  ct <- cor.test(i1, i2)
  return(c(ct$statistic, ct$p.value))
}

itypes <- function(x){
  i1 <- x[upper.tri(x)]
  i2 <- t(x)[upper.tri(x)] 
  
  comp <- sum(i1 < 0 & i2 < 0)
  mut <- sum(i1 > 0 & i2 > 0)
  pred <- sum(i1 > 0 & i2 < 0 | i1 < 0 & i2 > 0)
  amens <- sum(i1 < 0 & i2  == 0 | i1 == 0 & i2 < 0)
  comm <- sum(i1 > 0 & i2  == 0 | i1 == 0 & i2 > 0)
  
  return(c(comp = comp, mut = mut, pred = pred, amens = amens, comm = comm))
}

itypes.sp <- function(x){
  mm1 <- matrix(nrow = nrow(x), ncol = 3)
  for(i in 1:nrow(x)){
    i1 <- x[i, -i]
    i2 <- x[-i, i]
    
    comp <- sum(i1 < 0 & i2 < 0)
    mut <- sum(i1 > 0 & i2 > 0)
    pred <- sum(i1 > 0 & i2 < 0 | i1 < 0 & i2 > 0)
    amens <- sum(i1 < 0 & i2  == 0 | i1 == 0 & i2 < 0)
    comm <- sum(i1 > 0 & i2  == 0 | i1 == 0 & i2 > 0)
    
    mm1[i,] <- c(comp = comp, mut = mut, pred = pred)
  }
  return(mm1)
}

S <- 200
growth <- runif(S, .01, 1)
K <- quantile(1:100, rbeta(S, 1, 2))
mats <- get.adjacency(erdos.renyi.game(S, .1, "gnp"), sparse = F)  # emat1[[1]]
nINT <- sum(mats) 

mats[mats != 0] <- runif(nINT, -1, 1)
diag(mats) <- -1

spp <- list()
r2 <- list()
for(i in 1:1000){
  spp[[i]] <- c(sample(71:200, 20), sample(1:70, 30))
  isp <- spp[[i]]
  parms <- list(alpha = growth[isp], m = mats[isp,isp], k = K[isp])
  
  r2[[i]] <- ode(runif(length(isp), .1, 1), 1:1000, parms = parms, func = lvmod, events = list(func = ext1, time =  1:1000))
  matplot(r2[[i]][,-1], typ = "l")
  print(i)
}

use <- sapply(r2, nrow) == 1000 & sapply(r2, function(x) sum(tail(x, 1)[-1] > 0) > 0)
sapply(r2[use], function(x) sum(x[1000,-1] > 0))
eqcomm <- sapply(1:sum(use), function(x) spp[use][[x]][which(r2[use][[x]][1000,-1] > 0)])

eqmat <- matrix(0, nrow = sum(use), ncol = S)
for(i in 1:sum(use)){
  eqmat[i,eqcomm[[i]]] <- 1
}

matuse <- lapply(1:sum(use), function(i) mats[eqcomm[[i]], eqcomm[[i]]])

ity <- t(sapply(matuse, itypes))
plot(ity[,1]/ity[,2])

itySP <- lapply(matuse, itypes.sp)

summary(lm(unlist(lapply(1:sum(use), function(x) r2[use][[x]][1000,-1][r2[use][[x]][1000,-1] > 0]))~do.call(rbind, itySP)))


#############################
#############################
rem <- list()
pers <- c()
x <- 1
spp1 <- eqcomm[[x]]
initial1 <- mats[spp1, spp1]
for(i in 1:nrow(initial1)){
  
  parms <- list(alpha = growth[spp1[-i]], m = initial1[-i,-i])
  states <- r2[use][[x]][1000,-1][r2[use][[x]][1000,-1] > 0] 
  states <- states[-i]
  
  rem[[i]] <- ode(states, 1:1000, parms = parms, func = lvmod, events = list(func = ext1, time =  1:1000))
  matplot(rem[[i]][,-1], typ = "l")
  
  pers[i] <- sum(rem[[i]][1000,-1] > 0)/(nrow(initial1) -1)
  
  print(i)
}

init.biom <- mean(r2[use][[x]][1000,-1][r2[use][[x]][1000,-1] > 0])
delta.biom <- sapply(rem, function(x) mean(x[1000,-1] - init.biom))

lapply(rem, function(x) apply(x[,-1], 2, function(y) sd(y)/mean(y)))

