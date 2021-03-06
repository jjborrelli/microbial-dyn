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

# Function to get equilibrium communities
get_eq <- function(mats, times, INTs, Rmax = 1, Kval = 20, Ki = "val"){
  dyn <- list()
  mtmats <- list()
  grs <- list()
  kvs <- list()
  K.i <- Kval
  for(i in 1:length(mats)){
    gr <- runif(nrow(t1), -.1, Rmax)
    
    
    if(Ki == "rand"){
      rK <- rlnorm(nrow(t1))
      Kalt <- (rK/sum(rK))*K.i
      Kval <- Kalt
      parms <- list(alpha = gr, m = t1, K = Kval)
      
      # numerical integration of ODE, simulates dynamics of local community
      test <- ode(runif(nrow(t1), .1, .5), 1:times, parms = parms, func = lvmodK2, events = list(func = ext1, time =  1:times))
      
    }else if(Ki == "set"){
      dat1 <- apply(apply(mg1, 1, sort, decreasing = T), 1, median)[apply(apply(mg1, 1, sort, decreasing = T), 1, median)!=0]
      sam1 <- sample(dat1, nrow(t1), replace = T)
      #sam2 <- abs(rnorm(length(sam1), sam1, .1))
      Kval <- (sam1)/sum(sam1)*K.i 
      parms <- list(alpha = gr, m = t1, K = Kval)
      # numerical integration of ODE, simulates dynamics of local community
      test <- ode(runif(nrow(t1), .1, .5), 1:times, parms = parms, func = lvmodK2, events = list(func = ext1, time =  1:times))
    }else if(Ki == "val"){
      parms <- list(alpha = gr, m = t1, K = Kval)
      # numerical integration of ODE, simulates dynamics of local community
      test <- ode(runif(nrow(t1), .1, .5), 1:times, parms = parms, func = lvmodK, events = list(func = ext1, time =  1:times))
    }else if(Ki == "none"){
      parms <- list(alpha = gr, m = t1, K = Kval)
      # numerical integration of ODE, simulates dynamics of local community
      test <- ode(runif(nrow(t1), .1, .5), 1:times, parms = parms, func = lvmod, events = list(func = ext1, time =  1:times))
    }
    
    if(nrow(test) == times){
      dyn[[i]] <- test[,-1]
      mtmats[[i]] <- t1
      grs[[i]] <- gr
      kvs[[i]] <- Kval
    }else{
      dyn[[i]] <- NA
      mtmats[[i]] <- NA
      grs[[i]] <- NA
      kvs[[i]] <- NA
    }
    
    #matplot(test[,-1], type = "l", main = i)
    #print(length(dyn[!is.na(dyn)]))
  }
  
  ncomm <- sum(!is.na(dyn))
  wcomm <- which(!is.na(dyn))
  ## initial interaction matrices for communities that worked
  mtmats <- mtmats[!is.na(dyn)]
  ## growth rates of spp for communities that worked
  grs <- grs[!is.na(dyn)]#[sapply(dyn[!is.na(dyn)], function(x) sum(is.na(x)) == 0)]
  ## dynamics of communities that worked
  dyn <- dyn[!is.na(dyn)]#[sapply(dyn[!is.na(dyn)], function(x) sum(is.na(x)) == 0)]
  ## species with positive biomass
  spp1 <- lapply(dyn, function(x) (x[times,] != 0))
  ## equilibrium matrices
  #eqmat <- lapply(1:ncomm, function(x) mtmats[[x]][spp1[[x]], spp1[[x]]])
  eqgrs <- lapply(1:ncomm, function(x) grs[[x]][spp1[[x]]])
  eqst <- lapply(1:ncomm, function(x) dyn[[x]][times, spp1[[x]]])
  eqkv <- lapply(1:ncomm, function(x){
    if(length(kvs[!is.na(kvs)][[x]]) == 1){
      return(kvs[!is.na(kvs)][[x]]/sum(spp1[[x]]))
    }else{
      return(kvs[!is.na(kvs)][[x]][spp1[[x]]])
    }
  })  
  
  #return(list(spp = spp1, eqm = eqmat, eqgr = eqgrs, eqst = eqst, eqkv = eqkv, wrk = wcomm))
  return(list(spp = spp1, eqgr = eqgrs, eqst = eqst, eqkv = eqkv, wrk = wcomm))
}

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
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
############################# Generate Networks
library(igraph)
library(parallel)
library(doSNOW)


t0 <- Sys.time()

multityp <- lapply(1:1000, function(x){
  S <- sample(seq(500, 1000, 100), 1)
  p1 <- runif(1,0,1)
  p2 <- runif(1, p1, 1)
  c1 <- runif(1, .1, .3)
  mats <- get.adjacency(erdos.renyi.game(S, c1, "gnp", directed = F), sparse = F)
  tat <- mats*sample(c(-1,1,0), length(mats), replace = T, prob = c(p1,p2-p1,1-(p2)))
  return((tat))
})


#test <- lapply(1:5, function(x) multityp[[1]])
#sd1 <- c(.5, .5, .5, 1, 1)
#sd2  <- c(.5, 1, 2, 1, .5)
#for(i in 1:length(test)){
#  test[[i]] <- fill_mat(test[[i]], sdevp = sd1[i], sdevn = sd2[i])
#}
multityp.fill <- fill_mats(multityp, sdevn = 2, sdevp = .5)

t1 <- Sys.time()
t1-t0

#hist(multityp.fill[[1]][multityp.fill[[1]]!=0], breaks = 10)

t2 <- Sys.time()
#filepath1 <- "~/Documents/Data/"
filepath1 <- "D:/jjborrelli/parSADnoK/"

cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("filepath1", "lvmod", "lvmodK2", "ext1", "get_eq1", "fill_mat"))
registerDoSNOW(cl)

foreach(x = 15001:17000, .packages = c("deSolve", "R.utils", "igraph")) %dopar% {
  S <- sample(seq(500, 1000, 100), 1)
  p1 <- runif(1,0,1)
  p2 <- runif(1, p1, 1)
  c1 <- runif(1, .1, .2)
  mats <- get.adjacency(erdos.renyi.game(S, c1, "gnp", directed = F), sparse = F)
  multityp <- mats*sample(c(-1,1,0), length(mats), replace = T, prob = c(p1,p2-p1,1-(p2)))
  multityp.fill <- fill_mat(multityp, sdevp = 1, sdevn = 1)
  diag(multityp.fill) <- runif(length(diag(multityp.fill)), -2, 0)
  geq1 <- evalWithTimeout(get_eq1(multityp.fill, times = 2000, Ki = "val", Kval = 100, Rmax = 1), timeout = 600, onTimeout = "warning")
  #geq1 <- get_eq1(multityp.fill, times = 1000, Ki = "val", Kval = 20, Rmax = 1)
  if(is.character(geq1)){geq1 <- NA}
  saveRDS(geq1, file = paste(filepath1, "ge", x, ".rds", sep = ""))
  saveRDS(multityp.fill, file = paste(filepath1, "mat", x, ".rds", sep = "")) 
  # return(geq1)
}

stopCluster(cl)

t3 <- Sys.time()
t3 - t2

# 1:1500, diag = -5, r = 0.1 (all), sdp = .5, sdn = 2
# 1501:2500, diag = -5, r = 0:0.1, sdp = .5, sdn = 2
# 2501:3500, diag = -2, r = 0:0.1, sdp = .5, sdn = 2
# 3501:4500, diag = -2, r = 0:0.1, sdp = 1, sdn = 1


t4 <- Sys.time()
#filepath1 <- "~/Documents/Data/"
filepath2 <- "D:/jjborrelli/parSADtatoosh/"
tatoosh <- as.matrix(read.csv("C:/Users/jjborrelli/Desktop/GitHub/rKeystone/tatoosh.csv", header = F))

cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("tatoosh", "filepath2", "lvmodK", "lvmodK2", "ext1", "get_eq1", "fill_mat"))
registerDoSNOW(cl)

foreach(x = 6001:10000, .packages = c("deSolve", "R.utils", "igraph")) %dopar% {
  multityp <- tatoosh
  multityp.fill <- fill_mat(multityp, sdevp = 1, sdevn = 1)
  diag(multityp.fill) <- runif(length(diag(multityp.fill)), -2, 0)
  geq1 <- evalWithTimeout(get_eq1(multityp.fill, times = 1000, Ki = "val", Kval = 100, Rmax = 1), timeout = 360, onTimeout = "warning")
  #geq1 <- get_eq1(multityp.fill, times = 1000, Ki = "val", Kval = 20, Rmax = 1)
  if(is.character(geq1)){geq1 <- NA}
  saveRDS(geq1, file = paste(filepath2, "TATge", x, ".rds", sep = ""))
  saveRDS(multityp.fill, file = paste(filepath2, "TATmat", x, ".rds", sep = "")) 
  # return(geq1)
}

stopCluster(cl)

t5 <- Sys.time()
t5 - t4


################################################################################
################################################################################

t4 <- Sys.time()
#filepath1 <- "~/Documents/Data/"
filepath2 <- "D:/jjborrelli/parSADtatoosh/"
tatoosh <- as.matrix(read.csv("C:/Users/jjborrelli/Desktop/GitHub/rKeystone/tatoosh.csv", header = F))

cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("tatoosh", "filepath2", "lvmodK", "lvmodK2", "ext1", "get_eq1", "fill_mat"))
registerDoSNOW(cl)

foreach(x = 1:3000, .packages = c("deSolve", "R.utils", "igraph")) %dopar% {
  S <- sample(seq(500, 1000, 100), 1)
  p1 <- runif(1,0,1)
  p2 <- runif(1, p1, 1)
  pow <- rbeta(1, 4, 1)
  c1 <- runif(1, .1, .3)
  mats <- get.adjacency(barabasi.game(S, pow, m = S/10, directed = F), sparse = F)
  multityp <- mats*sample(c(-1,1,0), length(mats), replace = T, prob = c(p1,p2-p1,1-(p2)))
  multityp.fill <- fill_mat(multityp, sdevp = .5, sdevn = 2)
  diag(multityp.fill) <- runif(length(diag(multityp.fill)), -2, 0)
  geq1 <- evalWithTimeout(get_eq1(multityp.fill, times = 1000, Ki = "val", Kval = 20, Rmax = .1), timeout = 360, onTimeout = "warning")
  #geq1 <- get_eq1(multityp.fill, times = 1000, Ki = "val", Kval = 20, Rmax = 1)
  if(is.character(geq1)){geq1 <- NA}
  saveRDS(geq1, file = paste(filepath2, "HUBge", x, ".rds", sep = ""))
  saveRDS(multityp.fill, file = paste(filepath2, "HUBmat", x, ".rds", sep = "")) 
  # return(geq1)
}

stopCluster(cl)

t5 <- Sys.time()
t5 - t4


fna <- c()
for(x in 1:2000){
  fna[x] <- paste("HUBge", x, ".rds", sep = "")
}
which(!fna %in% list.files(filepath2))
