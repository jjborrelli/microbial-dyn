#save.image("C:/Users/jjborrelli/Desktop/GitHub/microbial-dyn/Data/sim.Rdata")
mg1 <- read.csv("C:/Users/jjborrelli/Desktop/GitHub/microbial-dyn/Data/m3sppdat.csv")

library(igraph)
library(NetIndices)
library(deSolve)
library(ggplot2)
library(rootSolve)
library(reshape2)
library(boot)
library(data.table)

#################################################################################################
#################################################################################################
#################################################################################################

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

# Function to detect extinction (prevents negative abundances)
ext1 <- function (times, states, parms){
  with(as.list(states), {
    states[states < 10^-5] <- 0 
    
    return(c(states))
  })
}

# Function to get equilibrium communities
get_eq <- function(mats, times, INTs, Rmax = 1, Kval = 20, Ki = FALSE){
  dyn <- list()
  mtmats <- list()
  grs <- list()
  kvs <- list()
  K.i <- Kval
  for(i in 1:length(mats)){
    t1 <- mats[[i]]
    diag(t1) <- 0  #-rbeta(length(diag(t1)), 1.1, 5)*5
    t1[t1 == 1] <- abs(rnorm(sum(t1 == 1), 0, .5))
    #abs(rnorm(sum(t1 == 1), mean(INTs), sd(INTs))) #runif(sum(t1 == 1), 0, 1) 
    t1[t1 == -1] <- -abs(rnorm(sum(t1 == -1), 0, .5))
    #-abs(rnorm(sum(t1 == -1), mean(INTs), sd(INTs))) # runif(sum(t1 == -1), -1, 0) 
    
    gr <- runif(nrow(t1), .1, Rmax)
    
   
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
    
    matplot(test[,-1], type = "l", main = i)
    print(length(dyn[!is.na(dyn)]))
  }
  
  ncomm <- sum(!is.na(dyn))
  wcomm <- which(!is.na(dyn))
  ## initial interaction matrices for communities that worked
  mtmats <- mtmats[!is.na(dyn)]#[sapply(dyn[!is.na(dyn)], function(x) sum(is.na(x)) == 0) & !sapply(dyn, function(x) sum(x < 0) > 0)]
  ## growth rates of spp for communities that worked
  grs <- grs[!is.na(dyn)]#[sapply(dyn[!is.na(dyn)], function(x) sum(is.na(x)) == 0)]
  ## dynamics of communities that worked
  dyn <- dyn[!is.na(dyn)]#[sapply(dyn[!is.na(dyn)], function(x) sum(is.na(x)) == 0)]
  ## species with positive biomass
  spp1 <- lapply(dyn, function(x) (x[times,] != 0))
  ## equilibrium matrices
  eqmat <- lapply(1:ncomm, function(x) mtmats[[x]][spp1[[x]], spp1[[x]]])
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



################################################################################################
################################################################################################
################################################################################################
#tatoosh <- as.matrix(read.csv("C:/Users/jjborrelli/Desktop/GitHub/rKeystone/tatoosh.csv", header = F))
tatoosh <- as.matrix(read.csv("~/Desktop/GitHub/rKeystone/tatoosh.csv", header = F))
S = 500
multityp <- lapply(1:200, function(x){
  p1 <- runif(1,0,1)
  p2 <- runif(1, p1, 1)
  c1 <- runif(1, .05, .3)
  mats <- get.adjacency(erdos.renyi.game(S, c1, "gnp", directed = F), sparse = F)
  #tat <- tatoosh*sample(c(1,-1), length(tatoosh), replace = T, prob = c(p1,1-p1))
  tat <- mats*sample(c(-1,1,0), length(mats), replace = T, prob = c(p1,p2-p1,1-(p2)))
  return((tat))
})

multitat <- lapply(1:200, function(x){
  p1 <- runif(1,0,1)
  tat <- tatoosh*sample(c(1,-1), length(tatoosh), replace = T, prob = c(p1,1-p1))
  #tat <- mats*sample(c(-1,1,0), length(mats), replace = T, prob = c(p1,p2-p1,1-(p2)))
  return((tat))
})

multihub <- lapply(1:200, function(x){
  p1 <- runif(1,0,1)
  p2 <- runif(1, p1, 1)
  pow <- rbeta(1, 4, 1)
  mats <- get.adjacency(barabasi.game(S, pow, m = 50, directed = F), sparse = F)
  #tat <- tatoosh*sample(c(1,-1), length(tatoosh), replace = T, prob = c(p1,1-p1))
  tat <- mats*sample(c(-1,1,0), length(mats), replace = T, prob = c(p1,p2-p1,1-(p2)))
  return((tat))
})

## Simulations   #### 

t1 <- Sys.time()
ge.mult <- get_eq(multityp, times = 1000, INTs = INTs, Ki = "val")
ge.mult2 <- get_eq(multityp, times = 1000, INTs = INTs, Ki = "rand")
ge.mult3 <- get_eq(multityp, times = 1000, INTs = INTs, Ki = "set")
t2 <- Sys.time()
t2-t1

ge.tat <- get_eq(multitat, times = 1000, INTs = INTs, Ki = "val")
ge.tat2 <- get_eq(multitat, times = 1000, INTs = INTs, Ki = "rand")
ge.tat3 <- get_eq(multitat, times = 1000, INTs = INTs, Ki = "set")
t3 <- Sys.time()

ge.hub <- get_eq(multihub, times = 1000, INTs = INTs, Ki = "val")
ge.hub2 <- get_eq(multihub, times = 1000, INTs = INTs, Ki = "rand")
ge.hub3 <- get_eq(multihub, times = 1000, INTs = INTs, Ki = "set")
t4 <- Sys.time()

## Fit Power Law   ####
fpl <- lapply(ge.mult$eqst, function(x) fit_power_law(x/sum(x)))
fpl2 <- lapply(ge.mult2$eqst, function(x) fit_power_law(x/sum(x)))
fpl3 <- lapply(ge.mult3$eqst, function(x) fit_power_law(x/sum(x)))

fplt <- lapply(ge.tat$eqst, function(x) fit_power_law(x/sum(x)))
fplt2 <- lapply(ge.tat2$eqst, function(x) fit_power_law(x/sum(x)))
fplt3 <- lapply(ge.tat3$eqst, function(x) fit_power_law(x/sum(x)))


fplh <- lapply(ge.hub$eqst, function(x) fit_power_law(x/sum(x)))
fplh2 <- lapply(ge.hub2$eqst, function(x) fit_power_law(x/sum(x)))
fplh3 <- lapply(ge.hub3$eqst, function(x) fit_power_law(x/sum(x)))

fpall <- list(fpl, fpl2, fpl3)
fptall <- list(fplt, fplt2, fplt3)
fphall <- list(fplh, fplh2, fplh3)

a1 <- lapply(fpall, function(x) sapply(x, "[", c(2,5)))
a2 <- lapply(fptall, function(x) sapply(x, "[", c(2,5)))
a3 <- lapply(fphall, function(x) sapply(x, "[", c(2,5)))

allpow <- do.call(rbind, list(do.call(rbind,sapply(1:3, function(x) cbind(as.matrix(t(a1[[x]])), comm = x, typ = "A"))), do.call(rbind,sapply(1:3, function(x) cbind(as.matrix(t(a2[[x]])), comm = x, typ = "B"))), do.call(rbind,sapply(1:3, function(x) cbind(as.matrix(t(a3[[x]])), comm = x, typ = "C")))))

head(allpow)
unlist(allpow[,"KS.stat"])

## Looking at fits    ####
median(sapply(fpl, "[[", 2)) # 6.5
median(sapply(fpl2, "[[", 2)) # 2.4
median(sapply(fpl3, "[[", 2)) #1.5

boxplot(log10(c(sapply(fpl, "[[", 2),sapply(fpl2, "[[", 2),sapply(fpl3, "[[", 2)))~c(rep(c("A", "B", "C"), c(length(fpl), length(fpl2), length(fpl3)))))

plot(sapply(fpl2, "[[", 2)~sapply(ge.mult2$eqm, function(x) sqrt(nrow(x) * sum(x != 0)/nrow(x)^2)))
cor.test(sapply(fpl3, "[[", 2),sapply(ge.mult3$eqm, function(x) sqrt(nrow(x) * sum(x != 0)/nrow(x)^2)))

fplk <- lapply(ge.mult$eqkv, fit_power_law)
fplk2 <- lapply(ge.mult2$eqkv, fit_power_law)
fplk3 <- lapply(ge.mult3$eqkv, fit_power_law)

alpha <- c(sapply(fplk2, "[[", 2) - sapply(fpl2, "[[", 2),sapply(fplk3, "[[", 2) - sapply(fpl3, "[[", 2))
beta <- rep(c("A", "B"), c(length(fplk2), length(fplk3)))
boxplot(alpha~beta)

sum(sapply(fpl, "[[", 5) <= 0.05)/length(fpl)
sum(sapply(fpl2, "[[", 5) <= 0.05)/length(fpl2)
sum(sapply(fpl3, "[[", 5) <= 0.05)/length(fpl3)



################################################################################################
################################################################################################
library(poweRlaw)
# make continuous power law object
m_ge <- conpl$new(ge.mult3$eqst[[3]])
# estimate lower bound
est <- estimate_xmin(m_ge)
# update distribution object
m_ge$setXmin(est)
# plot with best fit line
plot(m_ge)
lines(m_ge, col = 2, lwd = 2)
# fit lognormal
m_ge_ln <- conlnorm$new(ge.mult3$eqst[[3]])
est <- estimate_xmin(m_ge_ln)
m_ge_ln$setXmin(est)
lines(m_ge_ln, col = 3, lwd = 2)
# fit exp
m_ge_ex <- conexp$new(ge.mult3$eqst[[3]])
est <- estimate_xmin(m_ge_ex)
m_ge_ex$setXmin(est)
lines(m_ge_ex, col = 4, lwd = 2)

dist_ll(m_ge)
dist_ll(m_ge_ln)
dist_ll(m_ge_ex)
################################################################################################
################################################################################################
################################################################################################
## Real Data
### HMP 
otu2 <- read.csv("~/Desktop/otu_table_psn_v13.csv", row.names = 1)
metadat <- read.csv("~/Desktop/v13_map_uniquebyPSN.csv")

stoolsamp <- which(metadat$HMPbodysubsite == "Stool")
spptab <- colnames(otu2) %in% paste0("X",metadat[stoolsamp,]$SampleID)
otu3 <- otu2[-which(rowSums(otu2[,spptab]) == 0),spptab]

otu3 <- t(apply(otu3, 2, function(x) x/sum(x)))

fpl.hmp <- lapply(1:nrow(otu3), function(x) fit_power_law(otu3[x, otu3[x,] != 0]))
hist(sapply(fpl.hmp, "[[", 2))
fpl.al <- sapply(fpl.hmp, "[[", 2)
plot(fpl.al~apply(otu3, 1, function(x) sum(x!=0)))

test2 <- list()
for(i in 1:nrow(otu3)){
  test <- matrix(ncol = 4, nrow = nrow(otu3))
  for(j in 1:nrow(otu3)){
    d <- dist(rbind(otu3[i,], otu3[j,]))
    test[j,] <- c(i, j, d, fpl.al[i] - fpl.al[j])
  }
  test2[[i]] <- test
}

plot(do.call(rbind, test2)[,3:4])


### Caporaso

mg1 <- read.csv("~/Desktop/m3sppdat.csv")

sum(apply(mg1, 2, function(x) median(x[x !=0]))*20)
r1 <- rlnorm(ncol(mg1))
plot(apply(mg1, 2, function(x) median(x[x !=0]))*20, (r1/sum(r1)*20))
