### LIBRARIES
###


library(igraph)
library(deSolve)
library(ggplot2)
library(reshape2)

###
### FUNCTIONS
###


# Basic Lotka-Volterra model
lvmod <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha + state * parms$m %*% state 
    
    list(dB)
  })
}

# Function to detect extinction (prevents negative abundances)
ext1 <- function (times, states, parms){
  with(as.list(states), {
    states[states < 10^-5] <- 0 
    #if(sum(states >= 100) >= 1){states<-rep(0, length(states))} 
    return(c(states))
  })
}

make_community <- function(S, C, mu.ints, sd.ints){
  growth <- runif(S, .01, 1)     
  mats <- get.adjacency(erdos.renyi.game(S, C, "gnp", directed = F), sparse = F)  
  # note: altered to directed net, .15 = C
  nINT <- sum(mats)                                                  
  
  INTstr <- rnorm(nINT, mu.ints, sd.ints)
  
  mats[mats != 0] <- INTstr    
  mats[mats != 0] <- mats[mats != 0] * rbinom(nINT, size = 1, prob = .8)
  
  diag(mats) <- -rbeta(nrow(mats), 1.1, 5)*5
  
  return(mats)
}

sim_community <- function(times, state, parms, eq = lvmod, ex = ext1){
  out <- ode(state, times, parms = parms, func = eq, events = list(func = ex, time =  times))
  return(out[,-1])
}


get_eqcomm <- function(S, C, INTs, t1, plot = FALSE){
  cond <- FALSE
  while(!cond){
    c1 <- make_community(S, C, mean(INTs), sd(INTs))
    iA <- runif(nrow(c1))
    iP <- list(alpha = runif(nrow(c1), .1, 1), m = c1)
    
    sc1 <- sim_community(times = t1, state = iA, parms = iP)
    if(nrow(sc1) == max(t1)){cond <- TRUE}
  }
  if(plot){matplot(sc1, typ = "l")}
  return(list(comm.mat = c1, comm.dyn = sc1, init.parms = iP))
}

spp_remove <- function(sc1, iP, t1, track.progress = FALSE){
  scl <- list()
  eqS <- which(tail(sc1, 1) > 0)
  for(i in 1:length(eqS)){
    iA2 <- as.vector(tail(sc1, 1))
    iA2[eqS[i]] <- 0
    
    scl[[i]] <- sim_community(times = t1, state = iA2, parms = iP)
    if(track.progress){cat(i, "\t")}
  }
  return(scl)
}



SteinInt <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/ecomod-ints.csv", row.names = 1)
#SteinInt <- read.csv("C:/Users/jjborrelli/Desktop/GitHub/microbial-dyn/Data/ecomod-ints.csv", row.names = 1)
INTs <- c(SteinInt[upper.tri(SteinInt)],SteinInt[lower.tri(SteinInt)])


e1 <- get_eqcomm(S = 50, C = .3, INTs = INTs, t1 = 1:2000, plot = F)
sr1 <- spp_remove(sc1 = e1$comm.dyn, iP = e1$init.parms, t1 = 1:2000)


# Shuffle strengths
cond <- FALSE
while(!cond){
  ia1 <- runif(50)
  icomm <- make_community(50, .3, mean(INTs), sd(INTs))
  ip1 <- list(alpha = runif(50, .1, 1), m = icomm)
  simc1 <- sim_community(times = 1:2000, state = ia1, parms = ip1)
  
  if(nrow(simc1) == 2000){cond <- TRUE}
}

strt <- Sys.time()
mres <- matrix(nrow = 50, ncol = 50)
allmats <- list()
for(x in 1:50){
  cond <- FALSE
  while(!cond){
    c.mat <- icomm
    ist <- abs(icomm[icomm != 0])
    for(i in 1:500){
      sam1 <- sample(1:length(ist), 2)
      change1 <- runif(1, 0, abs(ist[sam1[1]]))
      ist[sam1[1]] <- ist[sam1[1]] - change1
      ist[sam1[2]] <- ist[sam1[2]] + change1
    }
    c.mat[icomm != 0] <- ist * sign(icomm[icomm != 0])
    ip1$m <- c.mat
    simc2 <- sim_community(times = 1:2000, state = ia1, parms = ip1)
    
    if(nrow(simc2[complete.cases(simc2),]) == 2000){cond <- T}
  }
  allmats[[x]] <- c.mat
  mres[,x] <- simc2[2000,]
}
ends <- Sys.time()
ends - strt

hist(sapply(allmats, mean))

strt <- Sys.time()
mres <- matrix(nrow = 50, ncol = 50)
allmats <- list()
for(x in 1:50){
  cond <- FALSE
  while(!cond){
    c.mat <- icomm
    c.mat[icomm != 0] <- sample(icomm[icomm != 0])
    ip1$m <- c.mat
    simc2 <- sim_community(times = 1:2000, state = ia1, parms = ip1)
    
    if(sum(complete.cases(simc2)) == 1){next}
    if(nrow(simc2[complete.cases(simc2),]) == 2000){cond <- T}
  }
  allmats[[x]] <- c.mat
  mres[,x] <- simc2[2000,]
}
ends <- Sys.time()
ends - strt
