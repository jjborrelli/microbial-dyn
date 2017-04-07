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

rad_info <- function(dfin, multityp.fill, thres = 5){
  eqa <- lapply(dfin[!is.na(dfin)], function(x) sort(x[x>0], decreasing = T))
  eqa <- eqa[!sapply(eqa, function(x) length(x) == 0)]
  eqm <- lapply((1:length(dfin[!is.na(dfin)])), function(x) multityp.fill[[x]][dfin[!is.na(dfin)][[x]] > 0, dfin[!is.na(dfin)][[x]] > 0])
  eqm <- eqm[!sapply(eqm, function(x) any(is.na(x)))]
  eqm1 <- lapply(eqm, itystr)
  eqm1 <- lapply(1:length(eqm1), function(x) data.frame(eqm1[[x]], N = nrow(eqm[[x]])))
  dstr <- sapply(eqm, function(x) mean(diag(x)))
  eqa2 <- lapply(eqa, get_abundvec, 100)
  fzN <- t(sapply(eqa2[sapply(eqa2, length) > thres], fzmod))
  
  testdat <- data.frame(prepdat(eqa2[sapply(eqa2, length) > thres], eqm1[sapply(eqa2, length) > thres], fzN[,"s"], fzN[,"r2"], dstr[sapply(eqa2, length) > thres]), abs = sapply(eqa2, sum)[sapply(eqa2, length) > thres])
  
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


multityp <- lapply(1:500, function(x){
  S <- 200
  p1 <- runif(1,0,1)
  p2 <- runif(1, p1, 1)
  c1 <- runif(1, .1, .3)
  mats <- get.adjacency(erdos.renyi.game(S, c1, "gnp", directed = F), sparse = F)
  tat <- mats*sample(c(-1,1,0), length(mats), replace = T, prob = c(p1,p2-p1,1-(p2)))
  return((tat))
})

multityp.fill <- fill_mats(multityp, sdevn = -2, sdevp = 2)
multityp.fill <- append(multityp.fill, multityp.fill)
#multihub.fill <- fill_mats(multihub, sdevn = -2, sdevp = 1)
s1 <- Sys.time()
dfin <- list()
dfin2 <- list()
dfin3 <- list()
dfin4 <- list()
cv <- c()
for(i in 1:100){
  #if(i < 501){
  #  diag(multityp.fill[[i]]) <- -2#(runif(nrow(multityp.fill[[i]]), sample(c(-1,-2,-3),1), 0))
  #}else{
    diag(multityp.fill[[i]]) <- (runif(nrow(multityp.fill[[i]]), sample(c(-1,-2,-3),1), 0))
  #}
  #diag(multityp.fill[[i]]) <- (runif(nrow(multityp.fill[[i]]), sample(c(-1,-2,-3),1), 0))
  par1 <- list(alpha = runif(nrow(multityp.fill[[i]]), 0,2), m = multityp.fill[[i]], K = 20)
  dyn <- ode(runif(nrow(multityp.fill[[i]]),.01,.05), times = 1:2000, func = lvmodK2, parms = par1, events = list(func = ext1, time =  1:2000))
  matplot(dyn[,-1], typ = "l", main = i)
  if(nrow(dyn) == 2000){
    dfin[[i]] <- dyn[2000,-1]
    dfin2[[i]] <- dyn[10,-1]
    dfin3[[i]] <- dyn[50,-1]
    dfin4[[i]] <- dyn[100,-1]
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
dfalt <- dfin
dfin <- dfalt
eqa <- lapply(dfin[!is.na(dfin)], function(x) sort(x[x!=0], decreasing = T))
eqa <- eqa[!sapply(eqa, function(x) length(x) == 0)]
eqm <- lapply((1:length(dfin[!is.na(dfin)])), function(x) multityp.fill[[x]][dfin[!is.na(dfin)][[x]] !=0, dfin[!is.na(dfin)][[x]]!=0])
eqm <- eqm[!sapply(eqm, function(x) any(is.na(x)))]
eqm1 <- lapply(eqm, itystr)
eqm1 <- lapply(1:length(eqm1), function(x) data.frame(eqm1[[x]], N = nrow(eqm[[x]])))
dstr <- sapply(eqm, function(x) mean(diag(x)))
eqa2 <- lapply(eqa, get_abundvec, 100)
fzN <- t(sapply(eqa2[sapply(eqa2, length) > 5], fzmod))

testdat <- data.frame(prepdat(eqa2[sapply(eqa2, length) > 5], eqm1[sapply(eqa2, length) > 5], fzN[,"s"], fzN[,"r2"], dstr[sapply(eqa2, length) > 5]), abs = sapply(eqa2, sum))

t2k <- rad_info(dfin, multityp.fill)
t10 <- rad_info(dfin2, multityp.fill)
t50 <- rad_info(dfin3, multityp.fill)
t100 <- rad_info(dfin4, multityp.fill)
pairs(cbind(t2k$sV, t10$sV, t50$sV, t100$sV))

median(t2k$sV)
median(t100$sV)
median(t50$sV)
median(t10$sV)


plot(testdat$sV~testdat$Nsp)
median(testdat$sV)
#predict(fit.init, testdat)
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


#################################################################################################
#################################################################################################
#################################################################################################

tdat <- (dzipf(1:100, N = 100, s = 1))
fztru <- fitzipf_r(tdat)

news <- c()
news2 <- c()
tmax <- c()
tsec <- c()
for(i in 1:100){
  newtdat.max <- c(runif(1, tdat[2], 1), tdat[-1])
  news[i] <- fitzipf_r(newtdat.max)@coef
  tmax[i] <- newtdat.max[1]
  
  newtdat.sec <- c(tdat[1], runif(1, tdat[3], tdat[1]), tdat[-c(1:2)])
  news2[i] <- fitzipf_r(newtdat.sec)@coef
  tsec[i] <- newtdat.sec[2]
}
diff1 <- tmax-tdat[1]
diff2 <- news-fztru@coef

diffA <- tsec - tdat[2]
diffB <- news2-fztru@coef

par(mfrow = c(1,2))
plot(diff2~diff1)
abline(a = 0, b = 1)

plot(diffB~diffA)
abline(a = 0, b = 1)

sv <- seq(.5, 3, .25)
ns <- seq(1110, 10, -100)
df1 <- list()
for(i in 1:11){
  dz1 <- dzipf(1:ns[i], N = 100, s = sv[i])
  df1[[i]] <- data.frame(dz = dz1, sval = sv[i], x = 1:ns[i])
}

ggplot(do.call(rbind,df1), aes(x = x, y = dz, col = factor(sval))) + geom_point() + facet_grid(~factor(sval), scales = "free_x")
ggplot(do.call(rbind,df1), aes(x = log10(x), y = log10(dz), col = factor(sval))) + geom_point() + facet_grid(~factor(sval))