mats <- readRDS("~/Documents/Data/mats.rds")
psd7 <- readRDS("~/Documents/Data/vdat2.rds")

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
library(R.utils)
####################################################################################################################
####################################################################################################################

wrk <- rep(NA, 5000)
eig <- c()
eig2 <- c()
spp <- list()
eqst <- list()
grs <- list()
eqkv <- c()
eqm <- list()
ico <- c()
Con <- c()
s0 <- Sys.time()
for(I in 761:5000){
  
  a.i <- psd7$eqa[[I]][order(as.numeric(names(psd7$eqa[[I]])))]
  m.i <- mats[[I]]
  k.i <- psd7$kv[[I]]
  r.i <- psd7$rs[[I]]
  
  par1 <- list(alpha = r.i, m = m.i, K = k.i)
  out <- evalWithTimeout(ode(y = a.i, times = 1:2000, func = lvmodK, parms = par1, events = list(func = ext1, time =  1:2000)), timeout = 120, onTimeout = "silent")
  if(is.null(out)){wrk[I] <- FALSE;next}
  matplot(out[,-1], typ = "l", main = I)
  
  jf <- jacobian.full(a.i, lvmodK, parms = par1)
  eig[I] <- max(Re(eigen(jf)$values))
  
  if(nrow(out) < 2000){wrk[I] <- FALSE;next}
  
  
  spp[[I]] <- out[2000,-1] > 10^-10
  if(sum(abs(sign(m.i[spp[[I]],spp[[I]]]))) == nrow(m.i[spp[[I]],spp[[I]]])){wrk[I] <- FALSE;next}
  
  ity <- itystr(m.i[spp[[I]],spp[[I]]])
  
  ico[I] <- is.connected(graph.adjacency(abs(sign(m.i[spp[[I]],spp[[I]]]))))
  Con[I] <- sum(abs(sign(m.i[spp[[I]],spp[[I]]])))/(sum(spp[[I]])*sum(spp[[I]]))
  
  eqst[[I]] <- out[2000,-1][spp[[I]]]
  grs[[I]] <- r.i[spp[[I]]]
  eqkv[I] <- k.i
  eqm[[I]] <- ity
  
  par2 <- list(alpha = grs[[I]], m = m.i[spp[[I]],spp[[I]]], K = k.i)
  jf2 <- jacobian.full(eqst[[I]], lvmodK, parms = par2)
  eig2[I] <- max(Re(eigen(jf2)$values))
  print(c(I, eig[[I]], eig2[[I]]))
  wrk[I] <- TRUE
}
s1 <- Sys.time()
s1-s0


# wrk[is.na(wrk)] <- TRUE
#25264
