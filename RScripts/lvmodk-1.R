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
    dB <- state * parms$alpha * (1 - state/(parms$K/sum(state > 10^-10))) + state * parms$m %*% state 
    #dB <- state * parms$alpha * (1 - state/(parms$K)) + state * parms$m %*% state 
    
    list(dB)
  })
}


# Lotka-Volterra model with set K
lvmodK2 <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha + state * parms$m %*% state 
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
    states[states < 10^-10] <- 0 
    
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
    t1[t1 == 1] <- abs(rnorm(sum(t1 == 1), .1, sdevp))
    t1[t1 == -1] <- -abs(rnorm(sum(t1 == -1), .1, sdevn))
    t2[[i]] <- t1 
    
    #t1[t1 == 1] <- runif(sum(t1 == 1), 0, sdevp)
    #t1[t1 == -1] <- runif(sum(t1 == -1), sdevn, 0)
    #t2[[i]] <- t1     
    
        
  }
  return(t2)
}



#################################################################################################
#################################################################################################
#################################################################################################

rad_info <- function(dfin, multityp.fill, thres = 5){
  eqa <- lapply(dfin[!is.na(dfin)], function(x) sort(x[x>0], decreasing = T))
  eqa <- eqa[!sapply(eqa, function(x) length(x) == 0)]
  eqm <- lapply((1:length(dfin[!is.na(dfin)])), function(x) multityp.fill[[x]][dfin[!is.na(dfin)][[x]] > 0, dfin[!is.na(dfin)][[x]] > 0])
  eqm <- eqm[!sapply(eqm, function(x) any(is.na(x)))]
  ico <- sapply(eqm, function(x) is.connected(graph.adjacency(abs(sign(x)))))
  eqm1 <- lapply(eqm, itystr)
  eqm1 <- lapply(1:length(eqm1), function(x) data.frame(eqm1[[x]], N = nrow(eqm[[x]])))
  dstr <- sapply(eqm, function(x) mean(diag(x)))
  eqa2 <- lapply(eqa, get_abundvec, 2000)
  fzN <- t(sapply(eqa2[sapply(eqa2, length) > thres], fzmod))
  
  testdat <- data.frame(prepdat(eqa2[sapply(eqa2, length) > thres], eqm1[sapply(eqa2, length) > thres], fzN[,"s"], fzN[,"r2"], dstr[sapply(eqa2, length) > thres]), abs = sapply(eqa2, sum)[sapply(eqa2, length) > thres], ico = ico[sapply(eqa2, length) > thres])
  
  return(testdat)
}


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


multityp <- lapply(1:100, function(x){
  S <- 500
  p1 <- runif(1,0,1)
  p2 <- runif(1, p1, 1)
  c1 <- runif(1, .1, .3)
  mats <- get.adjacency(erdos.renyi.game(S, c1, "gnp", directed = F), sparse = F)
  tat <- mats*sample(c(-1,1,0), length(mats), replace = T, prob = c(p1,p2-p1,1-(p2)))
  return((tat))
})

multityp.fill <- fill_mats(multityp, sdevn = 1, sdevp = .5)
#multityp.fill <- append(multityp.fill, multityp.fill)
#multihub.fill <- fill_mats(multihub, sdevn = -2, sdevp = 1)
s1 <- Sys.time()
dfin <- list()
dfin2 <- list()
dfin3 <- list()
dfin4 <- list()
alphas <- list()
cv <- c()
for(i in 1:100){
  #if(i < 501){
  #  diag(multityp.fill[[i]]) <- -2#(runif(nrow(multityp.fill[[i]]), sample(c(-1,-2,-3),1), 0))
  #}else{
    diag(multityp.fill[[i]]) <- (runif(nrow(multityp.fill[[i]]), -1, -.5))
  #}
  #diag(multityp.fill[[i]]) <- (runif(nrow(multityp.fill[[i]]), sample(c(-1,-2,-3),1), 0))
  par1 <- list(alpha = runif(nrow(multityp.fill[[i]]), -0.01,.1), m = multityp.fill[[i]], K = 100)
  dyn <-(ode(runif(nrow(multityp.fill[[i]]),.01,.2), times = 1:2000, func = lvmodK, parms = par1, events = list(func = ext1, time =  1:2000)))
  
  matplot(dyn[,-1], typ = "l", main = i)
  if(nrow(dyn) == 2000){
    dfin[[i]] <- dyn[2000,-1]
    dfin2[[i]] <- dyn[10,-1]
    dfin3[[i]] <- dyn[50,-1]
    dfin4[[i]] <- dyn[100,-1]
    alphas[[i]] <- par1$alpha
   cv[[i]] <- sd(apply(dyn[1900:2000,-1], 1, function(x) mean(x[x>0])))/mean(apply(dyn[1900:2000,-1], 1, function(x) mean(x[x>0])))
  }else{
    dfin[[i]] <- NA
    dfin2[[i]] <- NA
    dfin3[[i]] <- NA
    dfin4[[i]] <- NA
    cv[i] <- NA
  }
 
}
s2 <- Sys.time()
s2-s1


t2k <- rad_info(dfin, multityp.fill)
t10 <- rad_info(dfin2, multityp.fill)
t50 <- rad_info(dfin3, multityp.fill)
t100 <- rad_info(dfin4, multityp.fill)
pairs(cbind(t2k$sV, t10$sV, t50$sV, t100$sV))

median(t2k$sV)
median(t100$sV)
median(t50$sV)
median(t10$sV)

dfin[[24]]


mt1 <- multityp.fill[!is.na(dfin)][[2]][dfin[!is.na(dfin)][[2]] > 0,dfin[!is.na(dfin)][[2]] > 0]
ab1 <- dfin[!is.na(dfin)][[2]][dfin[!is.na(dfin)][[2]]>0]

library(rootSolve)
mre2 <- c()
for(i in 1:10000){
  jf1 <- jacobian.full(ab1, lvmodK, parms = list(alpha = runif(length(ab1),0,.1), m = mt1, K = 100))
  mre2[i] <- max(Re(eigen(jf1)$values))
}

#################################################################################################
#################################################################################################
#################################################################################################
x <- 2
y <- 4

m1 <- multityp.fill[[x]][dfin[[x]] > 0,dfin[[x]] > 0]
m2 <- multityp.fill[[y]][dfin[[y]] > 0,dfin[[y]] > 0]

m3 <- matrix(0, nrow = nrow(m1), ncol = ncol(m2))
m3[sample(1:length(m3), 10)] <- runif(10, -.1, 1)
m4 <- matrix(0, nrow = nrow(m2), ncol = ncol(m1))
m4[sample(1:length(m4), 10)] <- runif(10, -.1, 1)

m.all <- (cbind(rbind(m1, m4),rbind(m3,m2)))

init.ab <- c(dfin[[x]][dfin[[x]] > 0], dfin[[y]][dfin[[y]] > 0])
alph <- c(alphas[[x]][dfin[[x]] > 0], alphas[[y]][dfin[[y]] > 0])

par1 <- list(alpha = alph, m = m.all, K = 100)
dyn <-(ode(init.ab, times = 1:2000, func = lvmodK, parms = par1, events = list(func = ext1, time =  1:2000)))

matplot(dyn[,-1], typ = "l", main = i)

fzmod(dyn[2000,-1][dyn[2000,-1] > 0])
