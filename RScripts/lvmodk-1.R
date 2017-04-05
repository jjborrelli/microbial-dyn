library(igraph)
library(NetIndices)
library(deSolve)
library(ggplot2)
library(rootSolve)
library(reshape2)


#################################################################################################
#################################################################################################
#################################################################################################

# Lotka-Volterra model with evenly distributed K
lvmodK <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha * (1 - state/(parms$K/sum(state > 10^-5))) + state * parms$m %*% state 
    #dB <- state * parms$alpha * (1 - state/(parms$K)) + state * parms$m %*% state 
    
    list(dB)
  })
}


# Lotka-Volterra model with set K
lvmodK2 <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha * (1 - state/parms$K) + state * parms$m %*% state 
    #dB <- state * parms$alpha * (1 - state/(parms$K)) + state * parms$m %*% state 
    
    list(dB)
  })
}

lvmodN <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    imat <- parms$m
    for(i in 1:nrow(parms$m)){
      for(j in 1:ncol(parms$m)){
        if(i == j){next} 
        if(state[i] < 0){state[i] <- 0}
        if(state[j] < 0){state[j] <- 0}
        p1 <- (state[i]/sum(state))*(state[j]/sum(state)); imat[i,j] <- imat[i,j] * sample(c(1,0), 1, prob = c(p1, 1-p1))
      }
      }
    dB <- state * parms$alpha * (1 - state/(parms$K/sum(state > 0))) + state * imat %*% state 
    #dB <- state * parms$alpha * (1 - state/(parms$K)) + state * parms$m %*% state 
    
    list(dB)
  })
}


# Function to detect extinction (prevents negative abundances)
ext1 <- function (times, states, parms){
  with(as.list(states), {
    states[states < 10^-5] <- 0 
    
    return(c(states))
  })
}

# Function to detect extinction and add immigration
ext2 <- function (times, states, parms){
  with(as.list(states), {
    states[states < 10^-5] <- 0 
    #states <- states + (runif(length(states), 10^-5, 10^-2)*sample(c(0,1), length(states), prob = c(.9,.1), replace = T))
    states <- ceiling(states)
    return(c(states))
  })
}


# fill interaction matrices with strengths
fill_mats <- function(mats, sdevp = .5, sdevn = 1){
  t2 <- list()
  for(i in 1:length(mats)){
    t1 <- mats[[i]]
    diag(t1) <- 0  #-rbeta(length(diag(t1)), 1.1, 5)*5
    #t1[t1 == 1] <- abs(rnorm(sum(t1 == 1), -1, sdevp))
    #t1[t1 == -1] <- -abs(rnorm(sum(t1 == -1), -1, sdevn))
    #t2[[i]] <- t1 
    
    t1[t1 == 1] <- runif(sum(t1 == 1), 0, sdevp)
    t1[t1 == -1] <- runif(sum(t1 == -1), sdevn, 0)
    t2[[i]] <- t1     
    
        
  }
  return(t2)
}



#################################################################################################
#################################################################################################
#################################################################################################




multihub <- lapply(1:5, function(x){
  S = 500
  p1 <- runif(1,0,1)
  p2 <- runif(1, p1, 1)
  pow <- rbeta(1, 4, 1)
  mats <- get.adjacency(barabasi.game(S, pow, m = 50, directed = F), sparse = F)
  #tat <- tatoosh*sample(c(1,-1), length(tatoosh), replace = T, prob = c(p1,1-p1))
  tat <- mats*sample(c(-1,1,0), length(mats), replace = T, prob = c(p1,p2-p1,1-(p2)))
  return((tat))
})


multityp <- lapply(1:15, function(x){
  S <- 1000
  p1 <- runif(1,0,1)
  p2 <- runif(1, p1, 1)
  c1 <- runif(1, .1, .3)
  mats <- get.adjacency(erdos.renyi.game(S, c1, "gnp", directed = F), sparse = F)
  tat <- mats*sample(c(-1,1,0), length(mats), replace = T, prob = c(p1,p2-p1,1-(p2)))
  return((tat))
})

multityp.fill <- fill_mats(multityp, sdevn = -2, sdevp = 1)
#multihub.fill <- fill_mats(multihub, sdevn = -2, sdevp = 1)

dfin <- list()
dfin2 <- list()
for(i in 1:15){
  diag(multityp.fill[[i]]) <- (runif(nrow(multityp.fill[[i]]), -5, 0))
  par1 <- list(alpha = runif(nrow(multityp.fill[[i]]), 0,.1), m = multityp.fill[[i]], K = 20)
  dyn <- ode(runif(nrow(multityp.fill[[i]]),.01,.05), times = 1:1000, func = lvmodK, parms = par1, events = list(func = ext1, time =  1:1000))
  matplot(dyn[,-1], typ = "l", main = i)
  if(nrow(dyn) == 1000){dfin[[i]] <- dyn[1000,-1]}else{dfin[[i]] <- NA}
  if(nrow(dyn) == 1000){dfin2[[i]] <- apply(dyn[,-1], 2, mean)}else{dfin2[[i]] <- NA}
}

eqa <- sapply(dfin[!is.na(dfin)], function(x) sort(x[x!=0], decreasing = T))
eqm <- lapply((1:15)[!is.na(dfin)], function(x) multityp.fill[[x]][dfin[[x]] !=0, dfin[[x]]!=0])
eqm1 <- lapply(eqm, itystr)
eqm1 <- lapply(1:length(eqm1), function(x) data.frame(eqm1[[x]], N = nrow(eqm[[x]])))
dstr <- sapply(eqm, function(x) mean(diag(x)))
fzN <- t(sapply(eqa, fzmod))

testdat <- data.frame(prepdat(eqa, eqm1, fzN[,"s"], fzN[,"r2"], dstr), abs = sapply(eqa, sum))
predict(fit.init, testdat)
testdat$sV

dfinA <- list()
dfinA2 <- list()
for(i in 1:5){
  diag(multihub.fill[[i]]) <- (runif(nrow(multihub.fill[[i]]), -5, 0))
  par1 <- list(alpha = runif(nrow(multihub.fill[[i]]), 0, .1), m = multihub.fill[[i]], K = 20)
  dyn <- ode(runif(nrow(multihub.fill[[i]]),0,.04), times = 1:2000, func = lvmodK, parms = par1, events = list(func = ext1, time =  1:2000))
  matplot(dyn[,-1], typ = "l", main = i)
  if(nrow(dyn) == 2000){dfinA[[i]] <- dyn[2000,-1]}else{dfinA[[i]] <- NA}
  if(nrow(dyn) == 2000){dfinA2[[i]] <- apply(dyn[,-1], 2, mean)}else{dfinA2[[i]] <- NA}
}


