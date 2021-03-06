#save.image("C:/Users/jjborrelli/Desktop/GitHub/microbial-dyn/Data/sim.Rdata")
#mg1 <- read.csv("C:/Users/jjborrelli/Desktop/GitHub/microbial-dyn/Data/m3sppdat.csv")
#mg1 <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/m3sppdat.csv")


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


ext2 <- function (times, states, parms){
  with(as.list(states), {
    states[states < 10^-5] <- 0 
    states <- states + (runif(length(states), 10^-5, 10^-2)*sample(c(0,1), length(states), prob = c(.9,.1), replace = T))
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
    t1[t1 == 1] <- abs(rnorm(sum(t1 == 1), 0, 1))
    #abs(rnorm(sum(t1 == 1), mean(INTs), sd(INTs))) #runif(sum(t1 == 1), 0, 1) 
    t1[t1 == -1] <- -abs(rnorm(sum(t1 == -1), 0, 1))
    #-abs(rnorm(sum(t1 == -1), mean(INTs), sd(INTs))) # runif(sum(t1 == -1), -1, 0) 
    
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
  
  return(list(spp = spp1, eqm = eqmat, eqgr = eqgrs, eqst = eqst, eqkv = eqkv, wrk = wcomm))
  #return(list(spp = spp1, eqgr = eqgrs, eqst = eqst, eqkv = eqkv, wrk = wcomm))
}



################################################################################################
################################################################################################
################################################################################################
#tatoosh <- as.matrix(read.csv("C:/Users/jjborrelli/Desktop/GitHub/rKeystone/tatoosh.csv", header = F))
tatoosh <- as.matrix(read.csv("~/Desktop/GitHub/rKeystone/tatoosh.csv", header = F))
S = 500
multityp <- lapply(1:1000, function(x){
  p1 <- runif(1,0,1)
  p2 <- runif(1, p1, 1)
  c1 <- runif(1, .1, .3)
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

mg1 <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/m3sppdat.csv")

sum(apply(mg1, 2, function(x) median(x[x !=0]))*20)
r1 <- rlnorm(ncol(mg1))
plot(apply(mg1, 2, function(x) median(x[x !=0]))*20, (r1/sum(r1)*20))



##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
ext2 <- function (times, states, parms){
  with(as.list(states), {
    states[states < 10^-5] <- 0 
    #states <- states + (runif(length(states), 10^-5, 10^-2)*sample(c(0,1), length(states), prob = c(.9,.1), replace = T))
    #states <- states + (rlnorm(length(states), -5, 1)*sample(c(0,1), length(states), prob = c(.8,.2), replace = T))
    states[states > 20] <- 20
    return(c(states))
  })
}


multityp <- lapply(1:10, function(x){
  S <- 500# sample(seq(500, 1000, 100), 1)
  p1 <- runif(1,0,1)
  p2 <- runif(1, p1, 1)
  c1 <- runif(1, .1, .3)
  mats <- get.adjacency(erdos.renyi.game(S, c1, "gnp", directed = F), sparse = F)
  tat <- mats*sample(c(-1,1,0), length(mats), replace = T, prob = c(p1,p2-p1,1-(p2)))
  return((tat))
})

multihub <- lapply(1:5, function(x){
  S = 100
  p1 <- runif(1,0,1)
  p2 <- runif(1, p1, 1)
  pow <- rbeta(1, 4, 1)
  mats <- get.adjacency(barabasi.game(S, pow, m = 50, directed = F), sparse = F)
  #tat <- tatoosh*sample(c(1,-1), length(tatoosh), replace = T, prob = c(p1,1-p1))
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
multihub.fill <- fill_mats(multihub, sdevn = 2, sdevp = .5)

dfin <- list()
dfin2 <- list()
for(i in 1:10){
  diag(multityp.fill[[i]]) <- -abs(rnorm(500, -10, 4))
  par1 <- list(alpha = runif(nrow(multityp.fill[[i]])), m = multityp.fill[[i]], K = 20)
  dyn <- ode(runif(nrow(multityp.fill[[i]]),.01,.04), times = 1:1000, func = lvmodK, parms = par1, events = list(func = ext2, time =  1:1000))
  matplot(dyn[,-1], typ = "l", main = i)
  if(nrow(dyn) == 1000){dfin[[i]] <- dyn[1000,-1]}else{dfin[[i]] <- NA}
  if(nrow(dyn) == 1000){dfin2[[i]] <- apply(dyn[,-1], 2, mean)}else{dfin2[[i]] <- NA}
}


dfin <- list()
dfin2 <- list()
for(i in 1:10){
  par1 <- list(alpha = runif(nrow(multihub.fill[[i]])), m = multihub.fill[[i]], K = 20)
  dyn <- ode(runif(nrow(multihub.fill[[i]]),.01,.1), times = 1:1000, func = lvmodK, parms = par1, events = list(func = ext2, time =  1:1000))
  matplot(dyn[,-1], typ = "l", main = i)
  if(nrow(dyn) == 1000){dfin[[i]] <- dyn[1000,-1]}else{dfin[[i]] <- NA}
  if(nrow(dyn) == 1000){dfin2[[i]] <- apply(dyn[,-1], 2, mean)}else{dfin2[[i]] <- NA}
}


x <- 1
abund <- dfin[[x]][dfin[[x]] != 0]
#plot(sort(abund, decreasing = T), ylim = c(0, .25), xlim = c(0,200), typ = "l")
#plot(sort(get_abundvec(abund), decreasing = T), xlim = c(0,200))
plot(vegan::radfit(get_abundvec(abund)))
for(x in 2:5){
  if(is.na(dfin[[x]])){next}
  abund <- dfin[[x]][dfin[[x]] != 0]
  #points(sort(abund, decreasing = T), typ = "l")
  #points(sort(get_abundvec(abund), decreasing = T))
  #print(min(abund))
  plot(vegan::radfit(get_abundvec(abund)), main = i)
}



hist(dyn[1000,-1])
plot(sort(dyn[1000,-1], decreasing = T))
plot(as.vector(table(get_abundvec(dyn[1000,-1][dyn[1000,-1] != 0]/sum(dyn[1000,-1][dyn[1000,-1] != 0]), 10000))))




#################################################################################
#################################################################################

m.mats <- lapply(1:length(dfin), function(x){if(any(is.na(dfin[[x]]))){NA}else{multityp.fill[[x]][dfin[[x]] > 0,dfin[[x]] > 0]}})

lity1 <- lapply(multityp.fill, itystr)
lity2 <- lapply(m.mats[!is.na(m.mats)], itystr)

dat1 <- do.call(rbind, lity2)
ggplot(dat1, aes(x = num, y = m1)) + geom_point() + facet_wrap(~typ)

i = 3
plot(log10(sort(dfin2[[i]][dfin2[[i]] != 0], decreasing = T))~log10(seq(1,length(dfin2[[i]][dfin2[[i]] != 0]),1)))
abline(lm(log10(sort(dfin2[[i]][dfin2[[i]] != 0], decreasing = T))~log10(seq(1,length(dfin2[[i]][dfin2[[i]] != 0]),1))))
summary(lm(log10(sort(dfin2[[i]][dfin2[[i]] != 0], decreasing = T))~log10(seq(1,length(dfin2[[i]][dfin2[[i]] != 0]),1))))

is.na(dfin)

lity <- lapply(multityp.fill, itystr)
sapply(lity, function(x) x[,3])

#################################################################################
#################################################################################

lf1 <- list.files("~/Desktop/Data")
lf2 <- grep("ge", lf1)
lf3 <- grep("mat", lf1)
#for(i in 1:length(lf1[lf2])){
ist1 <- list()
for(i in 1:20){
  ge1 <- readRDS(paste("~/Desktop/Data/", "ge", i, ".rds", sep = ""))
  mat1 <- readRDS(paste("~/Desktop/Data/",  "mat", i, ".rds", sep = ""))
  if(any(is.na(ge1))){next}
  ist1[[i]] <- itystr(mat1[ge1$spp, ge1$spp])
}
