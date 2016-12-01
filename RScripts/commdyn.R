### Notes to self
# need to fix eigenvalue function to use jacobian, not interaction matrix
# re ran simulation using parametric bootstrapped data from Stein et al. 


### SAVED WORK
# last saved 10-25-16
# alt save 10-26-16
# ms save 11-8-16 == example3
# save.image("~/Desktop/simul-example4.Rdata") 
# load("~/Desktop/simul-example3.Rdata")



###
### LIBRARIES
###


library(igraph)
library(NetIndices)
library(deSolve)
library(ggplot2)
library(MuMIn)
library(rootSolve)
library(DAAG)

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

spints <- function(x){
  mm1 <- matrix(nrow = nrow(x), ncol = 10)
  for(i in 1:nrow(x)){
    i1 <- x[i, -i] # effect on self
    i2 <- x[-i, i] # effect on others
    
    pos1 <- sum(i1 > 0)
    neg1 <- sum(i1 < 0)
    non1 <- sum(i1 == 0)
    
    spos1 <- mean(i1[i1 > 0])
    sneg1 <- mean(i1[i1 < 0])
    
    
    pos2 <- sum(i2 > 0)
    neg2 <- sum(i2 < 0)
    non2 <- sum(i2 == 0)
    
    spos2 <- mean(i2[i2 > 0])
    sneg2 <- mean(i2[i2 < 0])
    
    
    mm1[i,] <- c(pos1, neg1, non1, spos1, sneg1, pos2, neg2, non2, spos2, sneg2)
  }
  mm1[is.nan(mm1)] <- 0
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
      pers[i] <- sum(rem[[i]][1000,-1] > 0)                         # what fraction of species have positive abundance
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
  # is.eq <- t(sapply(rem, function(x) if(nrow(x) == 1000){(x[1000,-1] > 0)*1}else{NA}))
  is.eq <- matrix(NA, length(spp1), length(spp1))
  for(i in 1:length(spp1)){if(nrow(rem[[i]] == 1000)){is.eq[i, -i] <- rem[[i]][1000, -1] > 0}else{is.eq[i,-i] <- rep(NA, length(spp1)[-1])}}
  
  
  # get data matrix for abundance change, variability, and persistence
  dat <- cbind(delta.biom, mean.vary, m.init.vary, pers)  
  
  return(list(dat, t(delta.eq), is.eq))
}

newfunc <- function(x){
  x
}

###
### SIMULATION
###


t.start <- Sys.time()                                               # begin timing of complete simulation run

S <- 200                                                            # total number of species in the pool 
growth <- runif(S, .01, 1)                                          # growth rate for each spp in the pool
K <- quantile(1:100, rbeta(S, 1, 2))                                # carrying capacities (not used)
mats <- get.adjacency(erdos.renyi.game(S, .2, "gnp", directed = F), sparse = F)   # create interaction matrix of all spp in pool
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
median(sapply(r2[use], function(x) sum(x[1000,-1] > 0))   )           # how many spp coexisting in each local comm
# which species are present in equilibrial communities
eqcomm <- sapply(1:sum(use), function(x) spp[use][[x]][which(r2[use][[x]][1000,-1] > 0)])

# how equilibrial are the communities?
dyn <- lapply(r2[use], function(x){x[x < 0] <- 0; x})               # get new object of dynamics list for runs that worked
cv.eq <- sapply(dyn, function(x) apply(x[990:1000,-1][,x[1000,-1] > 0], 2, sd)/(colMeans(x[990:1000, -1][,x[1000,-1] > 0])))
#cv.eq[is.nan(cv.eq)] <- 0
#range(colMeans(cv.eq))
hist(unlist(cv.eq))

# matrix of species found in each local equilibrium community
# can be used to determine compositional similarity of communities
eqmat <- matrix(0, nrow = sum(use), ncol = S)                       # initialize eqmat
for(i in 1:sum(use)){
  eqmat[i,eqcomm[[i]]] <- 1                                         # if the species is present in local comm i it gets a one, 0 otherwise
}

# initial interaction matrices
inmatuse <- lapply(1:sum(use), function(x) mats[spp[use][[x]],spp[use][[x]]])
# equilibrium interaction matrices for all iterations that didn't mess up
matuse <- lapply(1:sum(use), function(i) mats[eqcomm[[i]], eqcomm[[i]]])

# histogram of resulting matrix connectances
plot(sapply(matuse, function(x) sum(x != 0)/(nrow(x)*(nrow(x)-1))),sapply(inmatuse, function(x) sum(x != 0)/(nrow(x)*(nrow(x)-1))))

# compute frequency of interaction types in each equilibrium matrix
ity <- t(sapply(matuse, itypes))
plot(ity[,1]/ity[,2])                                               # plot ratio of competition:mutualism
iity <- t(sapply(inmatuse, itypes))

ity2 <- t(apply(ity, 1, function(x) x/sum(x)))
iity2 <- t(apply(iity, 1, function(x) x/sum(x)))

plot(iity2[,1], ity2[,1])

spi <- lapply(matuse, spints)
spi2 <- do.call(rbind, spi)
es <- spi2[,5] + spi2[,4]
eo <- spi2[,10] + spi2[,9]

ityinsp <- lapply(inmatuse, itypes.sp)
ityinsp2 <- do.call(rbind, ityinsp)
summary(glm(unlist(lapply(r2[use], function(x) x[1000, -1] > 0))~ityinsp2, family = "quasibinomial"))

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
fit4 <- glm(cbind(pers,rep(sapply(eqcomm, length), sapply(eqcomm, length))[ccak])~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+
              bet+close+neigh+ec+hub+pr,
            family = "binomial", data = mydat, na.action = "na.fail")
fit5 <- glm(eig~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+hub+pr,
            family = "gaussian", data = mydat, na.action = "na.fail")
#fit5.1 <- glm(eig~sp.id+n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+hub+pr,
#            family = "gaussian", data = mydat, na.action = "na.fail")

#fit6 <- glm(cbind(ceiling(pers*100), (100-ceiling(pers*100)))~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+hub+pr+sp.id,
#            family = "binomial", data = mydat, na.action = "na.fail")

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
fit4 <- glm(cbind(pers,rep(sapply(eqcomm, length), sapply(eqcomm, length))[ccak])~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred,
            family = "binomial", data = mydat, na.action = "na.fail")
fit5 <- glm(eig~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred, family = "gaussian", data = mydat, na.action = "na.fail")


d1.fit <- dredge(fit1)
d2.fit <- dredge(fit2)
d3.fit <- dredge(fit3)
d4.fit <- dredge(fit4)
d5.fit <- dredge(fit5)


model.avg(d1.fit, subset = delta < 5)
model.avg(d2.fit, subset = delta < 5)
model.avg(d3.fit, subset = delta < 5)
model.avg(d4.fit, subset = delta < 5)
model.avg(d5.fit, subset = delta < 5)

dmat1 <- matrix(c(colMeans(d1.fit[d1.fit$delta < 2,], na.rm = T),colMeans(d2.fit[d2.fit$delta < 2,], na.rm = T),colMeans(d3.fit[d3.fit$delta < 2,], na.rm = T),colMeans(d4.fit[d4.fit$delta < 2,], na.rm = T),colMeans(d5.fit[d5.fit$delta < 2,], na.rm = T)), nrow = 5, byrow = T)
colnames(dmat1) <- names(colMeans(d4.fit[d4.fit$delta < 2,]))
rownames(dmat1) <- c("biomass", "meanvary", "initvary", "persist", "eigen")
dmat1
dmat2 <- matrix(c((d1.fit[1,]),(d2.fit[1,]),(d3.fit[1,]),(d4.fit[1,]),(d5.fit[1,])), nrow = 5, byrow = T)
colnames(dmat2) <- names(colMeans(d4.fit[d4.fit$delta < 2,]))
rownames(dmat2) <- c("biomass", "meanvary", "initvary", "persist", "eigen")
dmat2


eq.abund <- unlist(lapply(dyn, function(x) x[1000,-1][x[1000,-1] > 0]))[ccak]
eq.abund2 <- (lapply(dyn, function(x) x[1000,-1][x[1000,-1] > 0]))

plot(lm(mydat$delta.biom~eq.abund)$residuals~eq.abund)


spi2 <- do.call(rbind, spi)
fit1 <- glm(mydat$delta.biom~spi2[ccak,1]+spi2[ccak,2]+spi2[ccak,3]+spi2[ccak,4]+spi2[ccak,5]+spi2[ccak,6]+spi2[ccak,7]+spi2[ccak,8]+spi2[ccak,9]+spi2[ccak,10], family = "gaussian", data = mydat, na.action = "na.fail")

dat2 <- cbind(rowMeans(spi2[,c(4,5)]),rowMeans(spi2[,c(9,10)]))[ccak,]
fit1 <- glm(mydat$delta.biom~dat2[,1]+dat2[,2], family = "gaussian", na.action = "na.fail")
summary(fit1)

########################
# 
ev.init <- c()
for(i in 1:sum(use)){
  dyna <- r2[use][[i]]
  isp <- eqcomm[[i]]
  eq.biom <- dyna[1000,-1][dyna[1000,-1] > 0]
  j1 <- jacobian.full(eq.biom, lvmod, parms = list(alpha = growth[isp], m = mats[isp,isp]))
  
  ev.init[i] <- max(Re(eigen(j1)$values))
}

evinit <- rep(ev.init, sapply(eqcomm, length))[ccak]
destab <- (mydat$eig >= evinit)

eab <- unlist(lapply(eq.abund2, function(x) x/sum(x)))[ccak]
icv <- rep(sapply(cv.eq, mean), sapply(eqcomm, length))[ccak]
neq <- rep(sapply(eqcomm, length), sapply(eqcomm, length))[ccak]
evi <- rep(ev.init, sapply(eqcomm, length))[ccak]

CI.abund <- (mydat$delta.biom/rep(sapply(eq.abund2, mean), sapply(eqcomm, length))[ccak]) * (1/eab)
CI.ivary <- ((icv - mydat$m.init.vary)/icv) * (1/eab) 
CI.pers <- ((neq - mydat$pers)/neq) * (1/eab)
CI.eig <- ((evi - mydat$eig)/evi) * (1/eab)


fit1 <- glm(CI.abund~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred, family = "gaussian", data = mydat, na.action = "na.fail")

fit3 <- glm(CI.ivary~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred, 
            family = "gaussian", data = mydat, na.action = "na.fail")
fit4 <- glm(CI.pers~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred,
            family = "gaussian", data = mydat, na.action = "na.fail")
fit5 <- glm(CI.eig~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred, family = "gaussian", data = mydat, na.action = "na.fail")


d1.fit <- dredge(fit1)
d3.fit <- dredge(fit3)
d4.fit <- dredge(fit4)
d5.fit <- dredge(fit5)


head(d1.fit)
head(d3.fit)
head(d4.fit)
head(d5.fit)

########################
# 
# maybe do something like ranking each species in each comm by impact measure (e.g., SpA is 1 in biomass change, 2 in eigenvalue) and comparing ranks across impact measures and whether spp have similar ranks in similar communities

ks1

impactrank <- lapply(ks1, function(x) apply(x, 2, function(y) order(abs(y), decreasing = T)))
mdist <- dist(eqmat[which(sapply(eqcomm, function(x) 27 %in% x))])

plot(allks[allks[,"sp.id"] == 27,5]~mdist[which(sapply(eqcomm, function(x) 27 %in% x)),which(sapply(eqcomm, function(x) 27 %in% x))])


resmat <- matrix(ncol = 3, nrow = 70)
for(i in 1:70){
  mdist <- dist(eqmat[which(sapply(eqcomm, function(x) i %in% x)),])
  d1 <- dist(allks[allks[,"sp.id"] == i,])
  mcor <- vegan::mantel(d1, mdist)
  resmat[i,] <- c(stat = mcor$statistic, pval = mcor$signif, numcom = sum(sapply(eqcomm, function(x) i %in% x)))
}

resmat


########################
# 
# Co-occurs? 
curve_ball<-function(m){
  RC=dim(m)
  R=RC[1]
  C=RC[2]
  hp=list()
  for (row in 1:dim(m)[1]) {hp[[row]]=(which(m[row,]==1))}
  l_hp=length(hp)
  for (rep in 1:5*l_hp){
    AB=sample(1:l_hp,2)
    a=hp[[AB[1]]]
    b=hp[[AB[2]]]
    ab=intersect(a,b)
    l_ab=length(ab)
    l_a=length(a)
    l_b=length(b)
    if ((l_ab %in% c(l_a,l_b))==F){
      tot=setdiff(c(a,b),ab)
      l_tot=length(tot)
      tot=sample(tot, l_tot, replace = FALSE, prob = NULL)
      L=l_a-l_ab
      hp[[AB[1]]] = c(ab,tot[1:L])
      hp[[AB[2]]] = c(ab,tot[(L+1):l_tot])}
    
  }
  rm=matrix(0,R,C)
  for (row in 1:R){rm[row,hp[[row]]]=1}
  rm
}
curving <- function(adjmat, n){
  newmat <- adjmat
  
  d.all <- matrix(nrow = length(dist(adjmat)), ncol = n)
  
  for(i in 1:n){
    newmat <- curve_ball(newmat)
    d.all[,i] <- dist(newmat)
  }
  return(d.all)
}

d1 <- dist((eqmat))

allD <- curving((eqmat), 500)

sig1 <- sapply(1:nrow(allD), function(x) sum(allD[x,] < d1[x])/500)
sum(sig1 < 0.05)/(371*370)
which(sig1 <= 0.05)

m1 <- matrix(0,nrow = 371, ncol = 371)
m1[lower.tri(m1)][which(sig1 <= 0.05)] <- 1



################################
################################

eqsp <- rep(sapply(eqcomm, length), sapply(eqcomm, length))


eqs1 <- lapply(1:sum(use), function(x) cbind(sp = spp[[x]], eq = r2[use][[x]][1000,-1], x))
eqs2 <- do.call(rbind, eqs1)

cors <- matrix(nrow = 200, ncol = 3)
dis1 <- list()
dis2 <- list()
for(i in 1:200){
  dat <- eqs2[eqs2[,"sp"] == i,]
  c1 <- combn(1:nrow(dat), 2)
  d1 <- c()
  d2 <- c()
  for(j in 1:ncol(c1)){
    temp <- dat[c1[,j],]
    d1[j] <- vegan::vegdist(temp[1:2,2], method = "euclidean")
    d2[j] <- sum(colSums(eqmat[temp[1:2,3],]) == 2)/sum(colSums(eqmat[temp[1:2,3],]) > 0)
  }
  dis1[[i]] <- d1
  dis2[[i]] <- d2
  
  ctest <- cor.test(d1, d2)
  cors[i,] <- c(ctest$estimate, ctest$p.value, nrow(dat)) 
}

cors

plot(cors[,1], cors[,3], pch = 20, col = factor(cors[,2] < 0.05))



################################################
## LDA ? 
library(MASS)
library(dplyr)

G1 <- (abs(CI.pers) > 10 & abs(CI.abund) > 10 & abs(CI.eig) > 10 & abs(CI.ivary) > 10)*1
G1 <- (abs(CI.abund) > 20)*1
G1 <- ((CI.pers) > 10 & (CI.abund) < 10 & (CI.eig) < 10 & (CI.ivary) < 10)*1
G1[((CI.pers) > 1 & (CI.abund) < -1 & (CI.eig) < -1 & (CI.ivary) < -1)] <- 2


quant1 <- .5
G1 <- (abs(CI.pers) > quantile(abs(CI.pers), probs = quant1) & abs(CI.abund) > quantile(abs(CI.abund), probs = quant1) & abs(CI.eig) > quantile(abs(CI.eig), probs = quant1) & abs(CI.ivary) > quantile(abs(CI.ivary), probs = quant1))*1
sum(G1)

length(which(G1 == 0))
length(which(G1 == 1))

newd <- select(mydat, n.comp:ec, pr)
newd <- (sapply(newd, function(x) (x-mean(x))/sd(x)))

newd<- as.data.frame(cbind(newd, G = G1[-which(mydat$pers == 0)]))
neval <- 250+(287-150)

errlin <- c()
for(i in 1:100){
  #s1 <- sample(which(newd2$G == 0), 4000)
  #s2 <- sample(which(newd2$G == 1), 150)
  testcase <- sample(1:nrow(newd2), 2000)
  
  fitLDA <- lda(G~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+pr, data = newd2)
  
  pLDA <- predict(fitLDA, newdata = newd[-testcase,])
  tabl <- table(newd$G[-testcase], pLDA$class)
  errlin[i]=(neval-sum(diag(tabl)))/neval
  
  print(i)
}
mean(errlin)




fitLDA <- lda(G~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred, data = (newd))
fitLDA

fitLDAvalues <- predict(fitLDA)
ldahist(fitLDAvalues$x[,1], g = G1[testcase])

co.id <- rep(1:sum(use), sapply(eqcomm, length))
newd$co.id <- co.id


quant1 <- .9
G1 <- (abs(CI.pers) > quantile(abs(CI.pers), probs = quant1) & abs(CI.abund) > quantile(abs(CI.abund), probs = quant1) & (CI.eig) > quantile((CI.eig), probs = quant1) & abs(CI.ivary) > quantile(abs(CI.ivary), probs = quant1))*1
sum(G1)
G2 <- (abs(CI.pers) > quantile(abs(CI.pers), probs = quant1))*1
sum(G2)
G3 <- (abs(CI.abund) > quantile(abs(CI.abund), probs = quant1))*1
sum(G3)
G4 <- ((CI.eig) > quantile((CI.eig), probs = quant1))*1
sum(G4)
G5 <- (abs(CI.ivary) > quantile(abs(CI.ivary), probs = quant1))*1
sum(G5)


length(which(G1 == 0))
length(which(G1 == 1))

newd <-  select(mydat, n.comp:ec, pr)# select(mydat[-which(mydat$pers == 0),], n.comp:ec, pr)
newd <- (sapply(newd, function(x) (x-mean(x))/sd(x)))

newd<- as.data.frame(cbind(newd, G = G1))# as.data.frame(cbind(newd, G = G1[-which(mydat$pers == 0)]))

newd2 <- mydat
newd2$G <- G1
newd2$G2 <- G2
newd2$G3 <- G3
newd2$G4 <- G4
newd2$G5 <- G5

test <- sample(1:4537, 4000)

fitCI <- glm(G~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+pr, data = newd2, family = "binomial", na.action = "na.fail")
fitCI2 <- glm(G2~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+pr, data = newd2, family = "binomial", na.action = "na.fail")
fitCI3 <- glm(G3~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+pr, data = newd2, family = "binomial", na.action = "na.fail")
fitCI4 <- glm(G4~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+pr, data = newd2, family = "binomial", na.action = "na.fail")
fitCI5 <- glm(G5~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+pr, data = newd2, family = "binomial", na.action = "na.fail")

dCI <- dredge(fitCI)
dCI2 <- dredge(fitCI2)
dCI3 <- dredge(fitCI3)
dCI4 <- dredge(fitCI4)
dCI5 <- dredge(fitCI5)

model.avg(dCI, subset = delta < 2)$coefficients[1,]
model.avg(dCI2, subset = delta < 2)$coefficients
model.avg(dCI3, subset = delta < 2)$coefficients
model.avg(dCI4, subset = delta < 2)$coefficients
model.avg(dCI5, subset = delta < 2)$coefficients

ma.boot <- lapply(1:5, function(x){m <- matrix(0, nrow = 200, ncol = 12); colnames(m) <- c("(Intercept)","n.comp","n.mut","n.pred","s.comp","s.mut","s.pred","bet","close","neigh","ec","pr");return(m)})
for(i in 1:200){
  resamp <- sample(1:nrow(newd2), nrow(newd2), replace = T)
  newd3 <- newd2[resamp, ]
  
  fitCI <- glm(G~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+pr, data = newd3, 
               family = "binomial", na.action = "na.fail")
  fitCI2 <- glm(G2~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+pr, data = newd3, 
                family = "binomial", na.action = "na.fail")
  fitCI3 <- glm(G3~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+pr, data = newd3, 
                family = "binomial", na.action = "na.fail")
  fitCI4 <- glm(G4~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+pr, data = newd3, 
                family = "binomial", na.action = "na.fail")
  fitCI5 <- glm(G5~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+pr, data = newd3, 
                family = "binomial", na.action = "na.fail")
  
  dCI <- dredge(fitCI)
  dCI2 <- dredge(fitCI2)
  dCI3 <- dredge(fitCI3)
  dCI4 <- dredge(fitCI4)
  dCI5 <- dredge(fitCI5)
  
  #ma1 <- model.avg(dCI, subset = delta < 2)$coefficients[1,]
  ma.boot[[1]][i,names(dCI[1,1:12])] <- unlist(dCI[1,1:12])     #[i,names(ma1)] <- ma1
  #ma2 <- model.avg(dCI2, subset = delta < 2)$coefficients[1,]
  ma.boot[[2]][i,names(dCI2[1,1:12])] <- unlist(dCI2[1,1:12])     #[i,names(ma2)] <- ma2
  #ma3 <-  model.avg(dCI3, subset = delta < 2)$coefficients[1,]
  ma.boot[[3]][i,names(dCI3[1,1:12])] <- unlist(dCI3[1,1:12])     #[i,names(ma3)] <- ma3
  #ma4 <- model.avg(dCI4, subset = delta < 2)$coefficients[1,]
  ma.boot[[4]][i,names(dCI4[1,1:12])] <- unlist(dCI4[1,1:12])     #[i,names(ma4)] <- ma4
  #ma5 <- model.avg(dCI5, subset = delta < 2)$coefficients[1,]
  ma.boot[[5]][i,names(dCI5[1,1:12])] <- unlist(dCI5[1,1:12])     #[i,names(ma5)] <- ma5
  
  print(i)
}


####################################
####################################
ev.init <- c()
for(i in 1:sum(use)){
  dyna <- r2[use][[i]]
  isp <- eqcomm[[i]]
  eq.biom <- dyna[1000,-1][dyna[1000,-1] > 0]
  j1 <- jacobian.full(eq.biom, lvmod, parms = list(alpha = growth[isp], m = mats[isp,isp]))
  
  ev.init[i] <- max(Re(eigen(j1)$values))
}

evinit <- rep(ev.init, sapply(eqcomm, length))[ccak]
destab <- (mydat$eig >= evinit)

eab <- unlist(lapply(eq.abund2, function(x) x/sum(x)))[ccak]
icv <- rep(sapply(cv.eq, mean), sapply(eqcomm, length))[ccak]
neq <- rep(sapply(eqcomm, length), sapply(eqcomm, length))[ccak]
evi <- rep(ev.init, sapply(eqcomm, length))[ccak]

CI.abund <- ((mydat$delta.biom/rep(sapply(eq.abund2, mean), sapply(eqcomm, length))[ccak]) * (1/eab))[(mydat$delta.biom/rep(sapply(eq.abund2, mean), sapply(eqcomm, length))[ccak]) < .1]
CI.ivary <- (((icv - mydat$m.init.vary)/icv) * (1/eab))[((icv - mydat$m.init.vary)/icv) < .1]
CI.pers <- (((neq - mydat$pers)/neq) * (1/eab))[((neq - mydat$pers)/neq) < .1]
CI.eig <- (((evi - mydat$eig)/evi) * (1/eab))[((evi - mydat$eig)/evi) < .1]

plot(((neq - mydat$pers)/neq)[((evi - mydat$eig)/evi) < .1]~eab[((evi - mydat$eig)/evi) < .1])


quant1 <- .9
G1 <- (abs(CI.pers) > quantile(abs(CI.pers), probs = quant1) & abs(CI.abund) > quantile(abs(CI.abund), probs = quant1) & (CI.eig) > quantile((CI.eig), probs = quant1) & abs(CI.ivary) > quantile(abs(CI.ivary), probs = quant1))*1
sum(G1)
G2 <- (abs(CI.pers) > quantile(abs(CI.pers), probs = quant1))*1
sum(G2)
G3 <- (abs(CI.abund) > quantile(abs(CI.abund), probs = quant1))*1
sum(G3)
G4 <- ((CI.eig) > quantile((CI.eig), probs = quant1))*1
sum(G4)
G5 <- (abs(CI.ivary) > quantile(abs(CI.ivary), probs = quant1))*1
sum(G5)


G1 <- (abs(CI.pers) > 100 & abs(CI.abund) > 100 & (CI.eig) > 100 & abs(CI.ivary) > 100)
sum(G1)
G2 <- (abs(CI.pers) > 100)*1
sum(G2)
G3 <- (abs(CI.abund) > 100)*1
sum(G3)
G4 <- ((CI.eig) > 100)*1
sum(G4)
G5 <- (abs(CI.ivary) > 100)*1
sum(G5)

t1 <- unlist(lapply(matuse, colSums))[ccak]
t2 <- unlist(lapply(matuse, rowSums))[ccak]
summary(glm(G~t1+t2, data = newd2, family = "binomial", na.action = "na.fail"))


newd2 <- mydat[eab > 10e-5,]
newd2$G <- destab[eab > 10e-5]#G1
newd2$G2 <- CI.pers#G2
newd2$G3 <- CI.abund#G3
newd2$G4 <- CI.eig#G4
newd2$G5 <- CI.ivary#G5

fitCI <- glm(G~n.comp+abs(s.comp)+n.mut+s.mut+n.pred+s.pred, data = newd2, family = "binomial", na.action = "na.fail")
fitCI2 <- glm(G2~n.comp+abs(s.comp)+n.mut+s.mut+n.pred+s.pred, data = newd2[((neq - mydat$pers)/neq) < .1,], family = "gaussian", na.action = "na.fail")
fitCI3 <- glm(G3~n.comp+abs(s.comp)+n.mut+s.mut+n.pred+s.pred, data = newd2[(mydat$delta.biom/rep(sapply(eq.abund2, mean), sapply(eqcomm, length))[ccak]) < .1,], family = "gaussian", na.action = "na.fail")
fitCI4 <- glm(G4~n.comp+abs(s.comp)+n.mut+s.mut+n.pred+s.pred, data = newd2[((evi - mydat$eig)/evi) < .1,], family = "gaussian", na.action = "na.fail")
fitCI5 <- glm(G5~n.comp+abs(s.comp)+n.mut+s.mut+n.pred+s.pred, data = newd2[((icv - mydat$m.init.vary)/icv) < .1,], family = "gaussian", na.action = "na.fail")

dCI <- dredge(fitCI)
dCI2 <- dredge(fitCI2)
dCI3 <- dredge(fitCI3)
dCI4 <- dredge(fitCI4)
dCI5 <- dredge(fitCI5)

ma1 <- model.avg(dCI, subset = delta < 2)
summary(ma1)
#bf1 <-  glm(G~n.comp, data = newd2, family = "binomial", na.action = "na.fail")
#cv1 <- cv.binary(bf1)

# pers
ma2 <- model.avg(dCI2, subset = delta < 2)
summary(ma2)
#bf2 <- glm(G2~n.comp+n.mut+s.comp+s.mut+s.pred+ec, data = newd2, family = "binomial", na.action = "na.fail")
#cv2 <- cv.binary(bf2)

# abund
ma3 <- model.avg(dCI3, subset = delta < 2)
summary(ma3)
#bf3 <- glm(G3~+n.mut+s.comp+s.mut+s.pred+close+pr, data = newd2, family = "binomial", na.action = "na.fail")
#cv3 <- cv.binary(bf3)

# eig
ma4 <- model.avg(dCI4, subset = delta < 2)
summary(ma4)
#bf4 <- glm(G4~n.mut+n.pred+s.comp+s.pred+ec, data = newd2, family = "binomial", na.action = "na.fail")
#cv4 <- cv.binary(bf4)

# init vary
ma5 <- model.avg(dCI5, subset = delta < 2)
summary(ma5)
#bf5 <- glm(G5~n.comp+n.mut+s.comp+s.mut+bet+close+ec+pr, data = newd2, family = "binomial", na.action = "na.fail")
#cv5 <- cv.binary(bf5)


length(unlist(eqcomm)[ccak][newd2$G2 == 1])
df0 <- data.frame(confint(ma1, level = .95), rownames(confint(ma1)), ma1$coefficients[1,], "Stab")
df1 <- data.frame(confint(ma2, level = .95), rownames(confint(ma2)), ma2$coefficients[1,], "Persistence")
df2 <- data.frame(confint(ma3, level = .95), rownames(confint(ma3)), ma3$coefficients[1,], "Change in Abundance")
df3 <- data.frame(confint(ma4, level = .95), rownames(confint(ma4)), ma4$coefficients[1,], "Local Stability")
df4 <- data.frame(confint(ma5, level = .95), rownames(confint(ma5)), ma5$coefficients[1,], "Initial Variation")

colnames(df1) <- c("lower", "upper", "met", "coef", "mod")
colnames(df2) <- c("lower", "upper", "met", "coef", "mod")
colnames(df3) <- c("lower", "upper", "met", "coef", "mod")
colnames(df4) <- c("lower", "upper", "met", "coef", "mod")

df1$sig <- df1$lower < 0 & df1$upper < 0 | df1$lower > 0 & df1$upper > 0
df2$sig <- df2$lower < 0 & df2$upper < 0 | df2$lower > 0 & df2$upper > 0
df3$sig <- df3$lower < 0 & df3$upper < 0 | df3$lower > 0 & df3$upper > 0
df4$sig <- df4$lower < 0 & df4$upper < 0 | df4$lower > 0 & df4$upper > 0

df1$impt <- ma2$importance[rownames(df1)]
df2$impt <- ma3$importance[rownames(df2)]
df3$impt <- ma4$importance[rownames(df3)]
df4$impt <- ma5$importance[rownames(df4)]

dfall <- rbind(df1, df2, df3, df4)

ggplot(dfall) + geom_segment(aes(x = lower, y = met, xend = upper, yend = met, col = sig)) + geom_vline(aes(xintercept = 0)) + 
  geom_point(aes(x = coef, y = met, col = sig)) + facet_wrap(~mod, scales = "free_x") + 
  scale_color_manual(name = "Significant", values = c("grey", "blue")) + xlab("Value") + ylab("Variable") + theme_bw()
ggsave(filename = "~/Desktop/modelpar.jpeg", width = 7, height = 5)

ggplot(dfall) + geom_segment(aes(x = lower, y = met, xend = upper, yend = met, col = sig)) + geom_vline(aes(xintercept = 0)) + 
  geom_point(aes(x = coef, y = met, col = sig)) + facet_wrap(~mod) + 
  scale_color_manual(name = "Significant", values = c("grey", "blue")) + scale_x_continuous(limits = c(-5, 5)) +
  xlab("Value") + ylab("Variable") + theme_bw()
ggsave(filename = "~/Desktop/modelparZOOM.jpeg", width = 7, height = 5)


ggplot(dfall, aes(x = met, y = impt, fill = sig)) + geom_bar(stat = "identity") + scale_fill_manual(name = "Significant", values = c("grey", "blue")) + facet_wrap(~mod) + ylab("Importance") + xlab("Variable") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = "~/Desktop/parimpt.jpeg", width = 7, height = 5)
####################################
####################################

####################################
####################################
##  Universal Interactions



eqsp <- rep(sapply(eqcomm, length), sapply(eqcomm, length))


eqs1 <- lapply(1:sum(use), function(x) cbind(sp = spp[use][[x]], eq = r2[use][[x]][1000,-1], x))
eqs2 <- do.call(rbind, eqs1)

c1 <- combn(1:181, 2)
nshare <- c()
nshare2 <- c()
d1 <- c()
d2 <- c()
d3 <- c()
d4 <- c()
d5 <- c()
for(i in 1:ncol(c1)){
  nshare[i] <- sum(colSums(eqmat[c1[,i],]) == 2)/sum(colSums(eqmat[c1[,i],]) != 0)
  shsp <- which(colSums(eqmat[c1[,i],]) == 2)
  l1 <- eqs1[c1[,i]]
  e1 <- l1[[1]][l1[[1]][,"sp"] %in% shsp, "eq"]
  e2 <- l1[[2]][l1[[2]][,"sp"] %in% shsp, "eq"]

  d1[i] <- dist(rbind(e1, e2))
  d2[i] <- dist(rbind(e1/sum(e1), e2/sum(e2)))
  
  l2 <- lapply(l1, function(x) x[,2] <- x[,2]/sum(x[,2]))
  e1 <- l2[[1]][l1[[1]][,"sp"] %in% shsp]
  e2 <- l2[[2]][l1[[2]][,"sp"] %in% shsp]
  
  nshare2[i] <- sum((e1 + e2)/2)
  
  d3[i] <- dist(rbind(e1, e2))
  d4[i] <- dist(rbind(e1/sum(e1), e2/sum(e2)))
  
  eA <- e1/sum(e1)
  eB <- e2/sum(e2)
  m <- (eA + eB)/2
  
  d5[i] <- sqrt((sum(eA * log10(eA/m)) + sum(eB * log10(eB/m)))/2)
}
range(nshare)
hist(nshare)
par(mfrow = c(2,2))
plot(d1~nshare, main = "Abs Abund")
plot(d2~nshare, main = "Scaled Abs Abund")
plot(d3~nshare, main = "Rel Abund")
plot(d4~nshare, main = "Scaled Rel Abund")

plot(d1~nshare2, main = "Abs Abund")
plot(d2~nshare2, main = "Scaled Abs Abund")
plot(d3~nshare2, main = "Rel Abund")
plot(d4~nshare2, main = "Scaled Rel Abund")

ggplot(data.frame(nshare2, d5), aes(x = nshare2, y = d5)) + geom_point(alpha = 0.1) + geom_smooth()
