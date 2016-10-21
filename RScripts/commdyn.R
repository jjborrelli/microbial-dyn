### SAVED WORK
# save.image("~/Desktop/simul-example.Rdata")
# load("~/Desktop/simul-example.Rdata")



###
### LIBRARIES
###


library(igraph)
library(NetIndices)
library(deSolve)
library(ggplot2)


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
    states[states < 10^-30] <- 0 
    if(sum(states >= 100) >= 1){states<-rep(0, length(states))} 
    return(c(states))
  })
}

# computes correlation between interaction strength pairs
icor <- function(x){
  i1 <- x[upper.tri(x)]
  i2 <- t(x)[upper.tri(x)] 
  
  ct <- cor.test(i1, i2)
  return(c(ct$statistic, ct$p.value))
}

# compute frequency of different interaction types for whole matrix
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

# compute frequency of different interaction types each spp participates in
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

# compute mean strength of each type interaction a spp participates in
istr.sp <- function(x){
  mm1 <- matrix(nrow = nrow(x), ncol = 3)
  for(i in 1:nrow(x)){
    i1 <- x[i, -i]
    i2 <- x[-i, i]
    
    comp <- which(i1 < 0 & i2 < 0)
    mut <- which(i1 > 0 & i2 > 0)
    pred <- which(i1 > 0 & i2 < 0 | i1 < 0 & i2 > 0)
    amens <- which(i1 < 0 & i2  == 0 | i1 == 0 & i2 < 0)
    comm <- which(i1 > 0 & i2  == 0 | i1 == 0 & i2 > 0)
    
    mm1[i,] <- c(comp = mean(c(i1[comp],i2[comp])), mut = mean(c(i1[mut],i2[mut])), pred = mean(c(i1[pred],i2[pred])))
    mm1[i,][is.nan(mm1[i,])] <- 0
  }
  return(mm1)
}

# function describes how species removal affects the local stability (magnitude of Max(Re(Lambda))) of the equilibrium comm
eigenkey <- function(mat){
  ev.init <- max(Re(eigen(mat)$values))                             # initial eigenvalue for the community
  ev.key <- c()
  for(i in 1:nrow(mat)){
    newmat <- mat[-i,-i]                                            # species removed from comm
    ev.key[i] <- max(Re(eigen(newmat)$values))                      # eigenvalue of perturbed community
  }
  return(ev.key)                                                    # return new eigenvalue (should this be the difference?)
}

# function to get network properties (for now just the network)
getgraph <- function(mat){
  mat[mat != 0] <- 1
  diag(mat) <- 0
  g <- graph.adjacency(mat)
  return(g)
}

# function to simulate impact of independent removals of all species in equilibrium comm
keystone <- function(x, dyn, eqcomm, mats, growth){
  rem <- list()                                                     # list for ode results
  pers <- c()                                                       # vector for persistence (could be done outside the loop)
  
  dyna <- dyn[[x]]                                                  # which initial community dynamics are we using
  spp1 <- eqcomm[[x]]                                               # what species from initial community are present at equilibrium
  initial1 <- mats[spp1, spp1]                                      # interaction matrix for equilibrium community
  for(i in 1:nrow(initial1)){
    
    parms <- list(alpha = growth[spp1[-i]], m = initial1[-i,-i])    # parameters of community with species i removed
    states <- dyna[1000,-1][dyna[1000,-1] > 0]                      # initial (eq) abundances from prior dynamics
    states <- states[-i]                                            # remove species i from initial (eq) abundances
    
    # simulation of new perturbed community
    rem[[i]] <- ode(states, 1:1000, parms = parms, func = lvmod, events = list(func = ext1, time =  1:1000))
    matplot(rem[[i]][,-1], typ = "l")                               # plot dynamics
    
    if(nrow(rem[[i]]) == 1000){                                     # if statement to determine if the run worked or crapped out
      pers[i] <- sum(rem[[i]][1000,-1] > 0)/(nrow(initial1) -1)     # what fraction of species have positive abundance
    }else{
      pers[i] <- NA                                                 # if run crashed gives NA
    }
    
    
    print(i)
  }
  
  # initial abundances for the removal sim
  init.biom <- dyna[1000,-1][dyna[1000,-1] > 0]
  # change in mean abundance following each removal
  delta.biom <- sapply(rem, function(x) if(nrow(x) == 1000){mean(x[1000,-1] - mean(init.biom))}else{NA})
  # change in equilibrium abundance following removal for each species
  delta.eq <- sapply(1:length(rem), function(x) if(nrow(rem[[x]]) == 1000){rem[[x]][1000,-1] - init.biom[-x]}else{rep(NA, length(init.biom[-1]))})
  
  # get coefficient of variation for each species following each removal (last 200 timesteps)
  vary <- lapply(rem, function(x) if(nrow(x) == 1000){apply(x[800:1000,-1], 2, function(y) sd(y)/mean(y))}else{NA})
  # mean CV following each removal
  mean.vary <- sapply(vary, function(x){x[is.nan(x)] <- 0; mean(x)})
  # get coefficient of variation for each species following each removal (first 50 timesteps)
  init.vary <- lapply(rem, function(x) if(nrow(x) == 1000){apply(x[1:50,-1], 2, function(y) sd(y)/mean(y))}else{NA})
  # mean CV immediately following each removal
  m.init.vary <- sapply(init.vary, mean)
  
  # does each species have positive abundance at equilibrium following each removal
  is.eq <- t(sapply(rem, function(x) if(nrow(x) == 1000){(x[1000,-1] > 0)*1}else{NA}))
  
  # get data matrix for abundance change, variability, and persistence
  dat <- cbind(delta.biom, mean.vary, m.init.vary, pers)  
  
  return(list(dat, t(delta.eq), is.eq))
}



###
### SIMULATION
###


t.start <- Sys.time()                                               # begin timing of complete simulation run

S <- 200                                                            # total number of species in the pool 
growth <- runif(S, .01, 1)                                          # growth rate for each spp in the pool
K <- quantile(1:100, rbeta(S, 1, 2))                                # carrying capacities (not used)
mats <- get.adjacency(erdos.renyi.game(S, .1, "gnp"), sparse = F)   # create interaction matrix of all spp in pool
nINT <- sum(mats)                                                   # total number of interactions 

INTstr <- runif(nINT, -1, 1)                                        # sample interaction strengths 

mats[mats != 0] <- INTstr                                           # fill in interaction strengths 
diag(mats) <- -1                                                    # self limitation set to -1

## 
## Begin simulation of subsampled community dynamics

# initialize lists of species and dynamics
spp <- list()
r2 <- list()
# simulation
for(i in 1:1000){
  spp[[i]] <- c(sample(71:200, 20), sample(1:70, 30))               # sample species from the pool (two samples ensure some overlap)
  isp <- spp[[i]]                                                   # local species community
  parms <- list(alpha = growth[isp], m = mats[isp,isp], k = K[isp]) # named parameter list (growth rates and int mat; K not used in sim)
  
  # numerical integration of ODE, simulates dynamics of local community
  r2[[i]] <- ode(runif(length(isp), .1, 1), 1:1000, parms = parms, func = lvmod, events = list(func = ext1, time =  1:1000))
  matplot(r2[[i]][,-1], typ = "l")                                  # plot community dynamics
  print(i)                                                          # which iteration are we on again? 
}

# which runs did not get messed up
use <- sapply(r2, nrow) == 1000 & sapply(r2, function(x) sum(tail(x, 1)[-1] > 0) > 0)
sapply(r2[use], function(x) sum(x[1000,-1] > 0))                    # how many spp coexisting in each local comm
# which species are present in equilibrial communities
eqcomm <- sapply(1:sum(use), function(x) spp[use][[x]][which(r2[use][[x]][1000,-1] > 0)])

# how equilibrial are the communities?
cv.eq <- sapply(r2[use], function(x) apply(x[990:1000,-1], 2, sd)/colMeans(x[990:1000, -1]))
cv.eq[is.nan(cv.eq)] <- 0
hist(colMeans(cv.eq))

# matrix of species found in each local equilibrium community
# can be used to determine compositional similarity of communities
eqmat <- matrix(0, nrow = sum(use), ncol = S)                       # initialize eqmat
for(i in 1:sum(use)){
  eqmat[i,eqcomm[[i]]] <- 1                                         # if the species is present in local comm i it gets a one, 0 otherwise
}

# equilibrium interaction matrices for all iterations that didn't mess up
matuse <- lapply(1:sum(use), function(i) mats[eqcomm[[i]], eqcomm[[i]]])

# compute frequency of interaction types in each equilibrium matrix
ity <- t(sapply(matuse, itypes))
plot(ity[,1]/ity[,2])                                               # plot ratio of competition:mutualism

# get frequency of interaction types for each species in each equilibrial comm
itySP <- lapply(matuse, itypes.sp)
# get mean strength of each interaction type species participate in
istrSP <- lapply(matuse, istr.sp)

# quick test to see if interaction participation influences equilibrium abundance
summary(lm(unlist(lapply(1:sum(use), function(x) r2[use][[x]][1000,-1][r2[use][[x]][1000,-1] > 0]))~do.call(rbind, itySP)))


eigkey <- lapply(matuse, eigenkey)                                  # how each spp removal affects eigenval of comm
summary(lm(unlist(eigkey)~do.call(rbind, itySP)))                   # relationship btwn eig and interaction types
summary(lm(unlist(eigkey)~do.call(rbind, istrSP)))                  # relationship btwn eig and interaction strengths


allg <- lapply(matuse, getgraph)                                    # get the network for each local eq comm
plot(unlist(eigkey)~unlist(sapply(allg, degree)))                   # look at relationship between degree and eig

t.simend <- Sys.time()                                              # note time initial comm sim ends



###
### CHECK FOR KEYSTONE SPECIES
###

t.key <- Sys.time()                                                 # note time keystone simm starts


dyn <- r2[use]                                                      # get new object of dynamics list for runs that worked

# quick test of function
#system.time(
#ks2 <- keystone(2, dyn = r2[use], eqcomm, mats)
#)


# Simulated removal for each species in each equilibrium community
ks1 <- list()
ks2 <- list()
ks3 <- list()
for(i in 1:sum(use)){
  KS <- keystone(i, dyn = r2[use], eqcomm, mats, growth)            # keystone species simulation
  ks1[[i]] <- KS[[1]]                                               # biomass variability and persistence
  ks2[[i]] <- KS[[2]]                                               # change in spp biomass with removal
  ks3[[i]] <- KS[[3]]                                               # who went extinct
  
  cat(paste("\n ------------------|   ", i, "   |------------------ \n"))
}


t.end <- Sys.time()                                                 # what time does the simulation end
t.end - t.start                                                     # total time spent simulating from start to finish



###
### ANALYSIS
###


itySP2 <- do.call(rbind, itySP)                                     # get species level participation in interaction types
istrSP2 <- do.call(rbind, istrSP)                                   # get species level interaction type strengths
allks <- do.call(rbind, lapply(ks1, function(x) x[,1:4]))           # all biomass, variation, and persistence

## Correlations among different stability measures
stabi <- data.frame(allks[complete.cases(allks),], eigen = unlist(eigkey)[complete.cases(allks)])
stabi2 <- data.frame(allks[complete.cases(allks),][allks[complete.cases(allks),4] > 0,], eigen = unlist(eigkey)[complete.cases(allks)][allks[complete.cases(allks),4] > 0])
ggpairs(stabi2)


plot(allks[complete.cases(allks),1][allks[complete.cases(allks),4] != 0] ~ istrSP2[complete.cases(allks),2][allks[complete.cases(allks),4] != 0])
abline(fit1, col = "blue")

p.key <- sapply(ks1, function(x) which.min(x[,4][which(x[,4] > 0)]))
eqcomm
allkeys <- sapply(1:sum(use), function(x) eqcomm[[x]][p.key[x]])
ggplot(data.frame(x = allkeys), aes(x = x)) + geom_bar() 


destab.sp <- lapply(1:sum(use), function(x) eqcomm[[x]][which(eigkey[[x]] > 0)])
destab <- lapply(1:sum(use), function(x) istrSP[[x]][which(eigkey[[x]] > 0),])

cdist <- dist(eqmat)
dim(as.matrix(cdist))
allkeys
table(allkeys)
cdist[which(as.character(allkeys) %in% names(which.max(table(allkeys))))]
mean(cdist)
hist(cdist)
