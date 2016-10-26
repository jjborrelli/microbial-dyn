### Notes to self
# need to fix eigenvalue function to use jacobian, not interaction matrix
# re ran simulation using parametric bootstrapped data from Stein et al. 


### SAVED WORK
# last saved 10-25-16
# save.image("~/Desktop/simul-example.Rdata") 
# load("~/Desktop/simul-example.Rdata")



###
### LIBRARIES
###


library(igraph)
library(NetIndices)
library(deSolve)
library(ggplot2)
library(MuMIn)
library(rootSolve)


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
  mm1 <- matrix(nrow = nrow(x), ncol = 5)
  for(i in 1:nrow(x)){
    i1 <- x[i, -i]
    i2 <- x[-i, i]
    
    comp <- sum(i1 < 0 & i2 < 0)
    mut <- sum(i1 > 0 & i2 > 0)
    pred <- sum(i1 > 0 & i2 < 0 | i1 < 0 & i2 > 0)
    amens <- sum(i1 < 0 & i2  == 0 | i1 == 0 & i2 < 0)
    comm <- sum(i1 > 0 & i2  == 0 | i1 == 0 & i2 > 0)
    
    mm1[i,] <- c(comp = comp, mut = mut, pred = pred, amens = amens, comm = comm)
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
eigenkey <- function(mat, growth, isp, dyna){
  eq.biom <- dyna[1000,-1][dyna[1000,-1] > 0]
  j1 <- jacobian.full(eq.biom, lvmod, parms = list(alpha = growth[isp], m = mat[isp,isp]))
  
  ev.init <- max(Re(eigen(j1)$values))                              # initial eigenvalue for the community
  ev.key <- c()
  for(i in 1:length(isp)){
    ispR <- isp[-i]                                                 # species removed from comm
    
    j2 <- jacobian.full(eq.biom[-i], lvmod, parms = list(alpha = growth[ispR], m = mat[ispR,ispR]))
    
    ev.key[i] <- max(Re(eigen(j2)$values))                          # eigenvalue of perturbed community
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
  
  # change in equilibrium abundance following removal for each species
  delta.eq <- sapply(1:length(rem), function(x) if(nrow(rem[[x]]) == 1000){rem[[x]][1000,-1] - init.biom[-x]}else{rep(NA, length(init.biom[-1]))})
  # change in mean abundance following each removal
  #delta.biom <- sapply(rem, function(x) if(nrow(x) == 1000){mean(x[1000,-1] - mean(init.biom))}else{NA})
  delta.biom <- colMeans(delta.eq)
  
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
mats <- get.adjacency(erdos.renyi.game(S, .2, "gnp", directed = T), sparse = F)   # create interaction matrix of all spp in pool
# note: altered to directed net, .15 = C
nINT <- sum(mats)                                                   # total number of interactions 

# using stein data
SteinInt <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/ecomod-ints.csv", row.names = 1)
INTs <- c(SteinInt[upper.tri(SteinInt)],SteinInt[lower.tri(SteinInt)])
INTstr <- rnorm(nINT, mean(INTs), sd(INTs))
#INTstr <- runif(nINT, -1, 1)                                        # sample interaction strengths 

mats[mats != 0] <- INTstr                                           # fill in interaction strengths 
#diag(mats) <- -1                                                    # self limitation set to -1
diag(mats) <- -rbeta(nrow(mats), 1.1, 5)*5
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
sum(use)
hist(sapply(r2[use], function(x) sum(x[1000,-1] > 0))   )           # how many spp coexisting in each local comm
# which species are present in equilibrial communities
eqcomm <- sapply(1:sum(use), function(x) spp[use][[x]][which(r2[use][[x]][1000,-1] > 0)])

# how equilibrial are the communities?
dyn <- lapply(r2[use], function(x){x[x < 0] <- 0; x})               # get new object of dynamics list for runs that worked
cv.eq <- sapply(dyn, function(x) apply(x[990:1000,-1], 2, sd)/(colMeans(x[990:1000, -1])))
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

# histogram of resulting matrix connectances
hist(sapply(matuse, function(x) sum(x != 0)/(nrow(x)*(nrow(x)-1))))

# compute frequency of interaction types in each equilibrium matrix
ity <- t(sapply(matuse, itypes))
plot(ity[,1]/ity[,2])                                               # plot ratio of competition:mutualism

# get frequency of interaction types for each species in each equilibrial comm
itySP <- lapply(matuse, itypes.sp)
# get mean strength of each interaction type species participate in
istrSP <- lapply(matuse, istr.sp)

# quick test to see if interaction participation influences equilibrium abundance
summary(lm(unlist(lapply(1:sum(use), function(x) r2[use][[x]][1000,-1][r2[use][[x]][1000,-1] > 0]))~do.call(rbind, itySP)))


# how each spp removal affects eigenval of comm
eigkey <- lapply(1:sum(use), function(x) eigenkey(mat = mats, growth = growth, isp = eqcomm[[x]], dyna = dyn[[x]]))        
summary(lm(unlist(eigkey)~do.call(rbind, itySP)))                   # relationship btwn eig and interaction types
summary(lm(unlist(eigkey)~do.call(rbind, istrSP)))                  # relationship btwn eig and interaction strengths


allg <- lapply(matuse, getgraph)                                    # get the network for each local eq comm
plot(unlist(eigkey)~unlist(sapply(allg, degree)))                   # look at relationship between degree and eig

betw <- lapply(allg, betweenness)                                   # get betweenness of each node
clocent <- lapply(allg, closeness)                                  # get closeness centrality
# get neighborhood of each spp going out 2 links
g.neighbors2 <- lapply(1:length(allg), function(x){sapply(graph.neighborhood(allg[[x]], 2), function(y) length(V(y)))})
ecent <- lapply(allg, function(x) eigen_centrality(x)$vector)       # get eigenvector centrality
hscore <- lapply(allg, function(x) hub_score(x)$vector)             # get hub score
p.rank <- lapply(allg, function(x) page_rank(x)$vector)             # get page rank algo



t.simend <- Sys.time()                                              # note time initial comm sim ends

###
### CHECK FOR KEYSTONE SPECIES
###

t.key <- Sys.time()                                                 # note time keystone simm starts



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

itySP3 <- t(apply(itySP2, 1, function(x) x/sum(x)))
itySP3[is.nan(itySP3)] <- 0

# put all data together in single matrix
allks <- cbind(allks, eig = unlist(eigkey), sp.id = unlist(eqcomm), n.comp = itySP2[,1], n.mut = itySP2[,2], n.pred = itySP2[,3], 
               n.amen = itySP2[,4], n.com = itySP2[,5], 
               s.comp = istrSP2[,1], s.mut = istrSP2[,2], s.pred = istrSP2[,3], bet = unlist(betw), close = unlist(clocent),
               neigh = unlist(g.neighbors2),  ec = unlist(ecent), hub = unlist(hscore), pr = unlist(p.rank))
ccak <- complete.cases(allks)                                       # only use complete cases

dim(allks[ccak,])

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




pcdat <- allks[ccak,1:5]                                            # pull out stability data for PCA
pcdat.norm <- apply(pcdat, 2, function(x) (x-mean(x))/sd(x))        # normalize data
pcA <- princomp(pcdat.norm)                                         # PCA on normalized stability data
summary(pcA)                                                        # summary, how much variation explained per axis
loadings(pcA)                                                       # what is on each axis
plot(pcA$scores[,1:2])                                              # PCA scores for first two axes of variation

pcdat2 <- allks[ccak, 6:18]
pcdat.norm2 <- apply(pcdat2, 2, function(x) (x-mean(x))/sd(x))
pcB <- princomp(pcdat.norm2)
summary(pcB)
loadings(pcB)



## Modeling
#### with MuMIn package


mydat <- as.data.frame(allks[ccak,])

fit1 <- glm(delta.biom~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+hub+pr, family = "gaussian", data = mydat, na.action = "na.fail")
fit2 <- glm(mean.vary~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+hub+pr, family = "gaussian", data = mydat, na.action = "na.fail")
fit3 <- glm(m.init.vary~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+hub+pr, 
            family = "gaussian", data = mydat, na.action = "na.fail")
fit4 <- glm(cbind(ceiling(pers*100), (100-ceiling(pers*100)))~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+hub+pr,
            family = "binomial", data = mydat, na.action = "na.fail")
fit5 <- glm(eig~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+hub+pr,
            family = "gaussian", data = mydat, na.action = "na.fail")
fit5.1 <- glm(eig~sp.id+n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+hub+pr,
            family = "gaussian", data = mydat, na.action = "na.fail")

fit6 <- glm(cbind(ceiling(pers*100), (100-ceiling(pers*100)))~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+hub+pr+sp.id,
            family = "binomial", data = mydat, na.action = "na.fail")

d1.fit <- dredge(fit1)
d2.fit <- dredge(fit2)
d3.fit <- dredge(fit3)
d4.fit <- dredge(fit4)
d5.fit <- dredge(fit5)
d6.fit <- dredge(fit6)

head(d1.fit)
head(d2.fit)
head(d3.fit)
head(d4.fit)
head(d5.fit)
head(d6.fit)

dmat1 <- matrix(c(colMeans(d1.fit[d1.fit$delta < 2,], na.rm = T),colMeans(d2.fit[d2.fit$delta < 2,], na.rm = T),colMeans(d3.fit[d3.fit$delta < 2,], na.rm = T),colMeans(d4.fit[d4.fit$delta < 2,], na.rm = T),colMeans(d5.fit[d5.fit$delta < 2,], na.rm = T)), nrow = 5, byrow = T)
colnames(dmat1) <- names(colMeans(d4.fit[d4.fit$delta < 2,]))
rownames(dmat1) <- c("biomass", "meanvary", "initvary", "persist", "eigen")
dmat1
dmat2 <- matrix(c((d1.fit[1,]),(d2.fit[1,]),(d3.fit[1,]),(d4.fit[1,]),(d5.fit[1,])), nrow = 5, byrow = T)
colnames(dmat2) <- names(colMeans(d4.fit[d4.fit$delta < 2,]))
rownames(dmat2) <- c("biomass", "meanvary", "initvary", "persist", "eigen")
dmat2

#####################################
#####################################
fit1 <- glm(delta.biom~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred, family = "gaussian", data = mydat, na.action = "na.fail")
fit2 <- glm(mean.vary~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred, family = "gaussian", data = mydat, na.action = "na.fail")
fit3 <- glm(m.init.vary~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred, 
            family = "gaussian", data = mydat, na.action = "na.fail")
fit4 <- glm(cbind(ceiling(pers*100), (100-ceiling(pers*100)))~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred,
            family = "binomial", data = mydat, na.action = "na.fail")
fit5 <- glm(eig~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred, family = "gaussian", data = mydat, na.action = "na.fail")


d1.fit <- dredge(fit1)
d2.fit <- dredge(fit2)
d3.fit <- dredge(fit3)
d4.fit <- dredge(fit4)
d5.fit <- dredge(fit5)


head(d1.fit)
head(d2.fit)
head(d3.fit)
head(d4.fit)
head(d5.fit)

dmat1 <- matrix(c(colMeans(d1.fit[d1.fit$delta < 2,], na.rm = T),colMeans(d2.fit[d2.fit$delta < 2,], na.rm = T),colMeans(d3.fit[d3.fit$delta < 2,], na.rm = T),colMeans(d4.fit[d4.fit$delta < 2,], na.rm = T),colMeans(d5.fit[d5.fit$delta < 2,], na.rm = T)), nrow = 5, byrow = T)
colnames(dmat1) <- names(colMeans(d4.fit[d4.fit$delta < 2,]))
rownames(dmat1) <- c("biomass", "meanvary", "initvary", "persist", "eigen")
dmat1
dmat2 <- matrix(c((d1.fit[1,]),(d2.fit[1,]),(d3.fit[1,]),(d4.fit[1,]),(d5.fit[1,])), nrow = 5, byrow = T)
colnames(dmat2) <- names(colMeans(d4.fit[d4.fit$delta < 2,]))
rownames(dmat2) <- c("biomass", "meanvary", "initvary", "persist", "eigen")
dmat2



########################
# 
# maybe do something like ranking each species in each comm by impact measure (e.g., SpA is 1 in biomass change, 2 in eigenvalue) and comparing ranks across impact measures and whether spp have similar ranks in similar communities

ks1

impactrank <- lapply(ks1, function(x) apply(x, 2, function(y) order(abs(y), decreasing = T)))
mdist <- dist(eqmat[which(sapply(eqcomm, function(x) 27 %in% x))])

plot(allks[allks[,"sp.id"] == 27,5]~mdist[which(sapply(eqcomm, function(x) 27 %in% x)),which(sapply(eqcomm, function(x) 27 %in% x))])


resmat <- matrix(ncol = 3, nrow = 70)
for(i in 1:70){
  mdist <- dist(eqmat[which(sapply(eqcomm, function(x) i %in% x))])
  d1 <- dist(allks[allks[,"sp.id"] == i,])
  mcor <- vegan::mantel(d1, mdist)
  resmat[i,] <- c(stat = mcor$statistic, pval = mcor$signif, numcom = sum(sapply(eqcomm, function(x) i %in% x)))
}

resmat
