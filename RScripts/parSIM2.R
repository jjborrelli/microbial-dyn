mats <- readRDS("~/Abundance/mats.rds")
psd7 <- readRDS("~/Abundance/vdat2.rds")

library(deSolve)

# Basic Lotka-Volterra model
lvmodK <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha * (1 - state/(parms$K/sum(state > 0))) + state * parms$m %*% state 
    #dB <- state * parms$alpha * (1 - state/(parms$K)) + state * parms$m %*% state 
    
    list(dB)
  })
}

lvmodK2 <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha * (1 - state/parms$K) + state * parms$m %*% state 
    #dB <- state * parms$alpha * (1 - state/(parms$K)) + state * parms$m %*% state 
    
    list(dB)
  })
}

lvmod <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha + state * parms$m %*% state 
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

fill_mats <- function(mats, sdevp = .5, sdevn = 1){
  t2 <- list()
  for(i in 1:length(mats)){
    t1 <- mats[[i]]
    diag(t1) <- 0  #-rbeta(length(diag(t1)), 1.1, 5)*5
    t1[t1 == 1] <- abs(rnorm(sum(t1 == 1), 0, sdevp))
    #abs(rnorm(sum(t1 == 1), mean(INTs), sd(INTs))) #runif(sum(t1 == 1), 0, 1) 
    t1[t1 == -1] <- -abs(rnorm(sum(t1 == -1), 0, sdevn))
    #-abs(rnorm(sum(t1 == -1), mean(INTs), sd(INTs))) # runif(sum(t1 == -1), -1, 0) 
    t2[[i]] <- t1 
  }
  return(t2)
}

fill_mat <- function(mat, sdevp = .5, sdevn = 1){
  t1 <- mat
  diag(t1) <- 0  #-rbeta(length(diag(t1)), 1.1, 5)*5
  t1[t1 == 1] <- abs(rnorm(sum(t1 == 1), 0, sdevp))
  #abs(rnorm(sum(t1 == 1), mean(INTs), sd(INTs))) #runif(sum(t1 == 1), 0, 1) 
  t1[t1 == -1] <- -abs(rnorm(sum(t1 == -1), 0, sdevn))
  #-abs(rnorm(sum(t1 == -1), mean(INTs), sd(INTs))) # runif(sum(t1 == -1), -1, 0)
  
  return(t1)
}

# Function to get equilibrium communitie
get_eq1 <- function(mat, times, Rmax = 1, Kval = 20, Ki = "val"){
  gr <- runif(nrow(mat), 0, Rmax)
  
  if(Ki == "rand"){
    rK <- rlnorm(nrow(mat))
    Kalt <- (rK/sum(rK))*K.i
    Kval <- Kalt
    parms <- list(alpha = gr, m = mat, K = Kval)
    kvs <- Kval[spp]
    # numerical integration of ODE, simulates dynamics of local community
    test <- ode(runif(nrow(mat), .1, .5), 1:times, parms = parms, func = lvmodK2, events = list(func = ext1, time =  1:times))
    
  }else if(Ki == "val"){
    parms <- list(alpha = gr, m = mat, K = Kval)
    # numerical integration of ODE, simulates dynamics of local community
    test <- ode(runif(nrow(mat), .1, .5), 1:times, parms = parms, func = lvmodK, events = list(func = ext1, time =  1:times))
    kvs <- Kval
  }else if(Ki == "none"){
    parms <- list(alpha = gr, m = mat, K = Kval)
    # numerical integration of ODE, simulates dynamics of local community
    test <- ode(runif(nrow(mat), .1, .5), 1:times, parms = parms, func = lvmod, events = list(func = ext1, time =  1:times))
  }
  
  if(nrow(test) == times){
    spp <- test[times, -1] != 0
    dyn <- test[,-1]
    eqst <- dyn[times, spp]
    gr <- gr[spp]
    
    return(list(spp = spp, eqst = eqst, grs = gr, eqkv = kvs))
  }else{
    return(NA)
  }
}

get_abundvec <- function(abund, N = 10000){
  r.ab <- abund/sum(abund)
  samp2 <- sample(1:length(abund), N, replace = T, prob = r.ab)
  t1 <- table(samp2)
  return(as.numeric(t1))
}

itystr <- function(x){
  if(sum(x) == 0){
    return(data.frame(typ = c("amensalism", "commensalism", "competition", "mutualism", "predation"), 
                      num = c(0,0,0,0,0), m1 = c(0,0,0,0,0), m2 = c(0,0,0,0,0)))
  }
  
  i1 <- x[upper.tri(x)]
  i2 <- t(x)[upper.tri(x)] 
  
  nonI <- i1 == 0 & i2 == 0
  
  ints <- cbind(i1 = apply(cbind(i1, i2), 1, max), i2 = apply(cbind(i1, i2), 1, min))#[!nonI,]
  
  inty <- vector(length = nrow(ints))
  
  inty[ints[,1] < 0 & ints[,2] < 0] <- "competition"
  inty[ints[,1] > 0 & ints[,2] > 0] <- "mutualism"
  inty[ints[,1] > 0 & ints[,2] < 0 | ints[,1] < 0 & ints[,2] > 0] <- "predation"
  inty[ints[,1] < 0 & ints[,2]  == 0 | ints[,1] == 0 & ints[,2] < 0] <- "amensalism"
  inty[ints[,1] > 0 & ints[,2]  == 0 | ints[,1] == 0 & ints[,2] > 0] <- "commensalism"
  inty[ints[,1] == 0 & ints[,2] == 0] <- "none"
  
  df1 <- data.frame(ints, inty)[inty != "none",]
  
  num <- vector(length = 5, mode = "numeric")
  umu <- vector(length = 5, mode = "numeric")
  lmu <- vector(length = 5, mode = "numeric")
  
  
  iL <- c("amensalism", "commensalism", "competition", "mutualism", "predation")
  num[iL %in% aggregate(df1$inty, list(df1$inty), length)$Group.1] <- aggregate(df1$inty, list(df1$inty), length)$x
  umu[iL %in% aggregate(df1$inty, list(df1$inty), length)$Group.1] <- aggregate(df1$i2, list(df1$inty), mean)$x
  lmu[iL %in% aggregate(df1$inty, list(df1$inty), length)$Group.1] <- aggregate(df1$i1, list(df1$inty), mean)$x
  
  
  df2 <- data.frame(typ = iL, num = num, m1 = umu, m2 = lmu)
  
  return(df2)
}
####################################################################################################################
####################################################################################################################
library(parallel)
library(doSNOW)
####################################################################################################################
####################################################################################################################
for(I in 1:length(mats)){
  
  a.i <- psd7$eqa[[I]][order(as.numeric(names(psd7$eqa[[I]])))]
  m.i <- mats[[I]]
  k.i <- psd7$kv[[I]]
  r.i <- psd7$rs[[I]]
  
  par1 <- list(alpha = r.i, m = m.i, K = k.i)
  out <- ode(y = a.i, times = 1:1000, func = lvmodK, parms = par1, events = list(func = ext1, time =  1:1000))
  matplot(out[,-1], typ = "l")
  
  jf <- jacobian.full(a.i, lvmodK, parms = par1)
  print(max(Re(eigen(jf)$values)))
}

#25264

cl <- makeCluster(n.cores, type = "FORK")
clusterExport(cl, c("lvmod", "lvmodK2", "ext1", "psd7", "mats"))
registerDoSNOW(cl)

RESULTS <- foreach(x = 1:length(mats), .packages = c("deSolve", "R.utils", "igraph")) %do% {
  a.i <- psd7$eqa[[x]][order(as.numeric(names(psd7$eqa[[x]])))]
  m.i <- mats[[x]]
  k.i <- psd7$kv[[x]]
  r.i <- psd7$rs[[x]]
  
  par1 <- list(alpha = r.i, m = m.i, K = k.i)
  out <- ode(y = a.i, times = 1:1800, func = lvmodK, parms = par1, events = list(func = ext1, time =  1:1800))
  
  spp <- out[,-1] > 10^-10
  ity <- itystr(m.i[spp,spp])
  
  ico <- is.connected(graph.adjacency(abs(sign(m.i[spp,spp]))))
  Con <- sum(abs(sign(m.i[spp,spp])))/(sum(spp)*sum(spp))
  
  return(list(spp = spp, eqst = out[,-1][spp], grs = r.i[spp], eqkv = k.i, eqm = ity, ico = ico))
}

save(RESULTS, "~/Abundance/dat2ALT.Rdata")

