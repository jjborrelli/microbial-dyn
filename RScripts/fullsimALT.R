############################################################################
#### Dynamics Functions

lvmodK <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha * (1 - state/(parms$K/sum(state > 10^-10))) + state * parms$m %*% state 
   
    list(dB)
  })
}

lvmod <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha + state * parms$m %*% state + 10^-5
    list(dB)
  })
}

ext1 <- function (times, states, parms){
  with(as.list(states), {
    states[states < 10^-10] <- 0 
    
    return(c(states))
  })
}


############################################################################


library(parallel)
library(doSNOW)


filepath1 <- "~/Abundance/Data/"
iter <- 1:2
cl <- makeCluster(detectCores() - 1)


clusterExport(cl, c("filepath1"))
registerDoSNOW(cl)

# Parameters
S <- 600        # Number of Species
beta1 <- 1      # first beta param - interactions
beta2 <- 4      # second beta param - interactions
d1p <- 2        # first beta param - diagonal
d2p <- 1        # second beta param - diagonal
Rmax <- .1      # maximum growth rate for unif dist
Kcomm <- 20     # community carrying capacity
tmax <- 2000    # max number of timesteps


iter = 50
reslst <- list()
parlst <- list()
inab <- matrix(nrow = 600, ncol = 50)
tfz <- list()
for(i in 41:iter){
  # Load necessary functions and libraries
  #source("~/Desktop/GitHub/microbial-dyn/RScripts/abund-fxns.R")

  cond <- FALSE
  while(!cond){
    # Create community intereaction network
    
    ## Probability of negative, positive, and no interaction
    p1 <- runif(1,0,1)
    p2 <- runif(1, p1, 1)
    
    ## Connectance
    c1 <- runif(1, .1, .8)
    
    ## Adjacency matrix
    mats <- get.adjacency(barabasi.game(S, 1, 100, directed = F), sparse = F)
    mats2 <- mats*sample(c(-1,1,0), length(mats), replace = T, prob = c(p1,p2-p1,1-(p2)))
    
    ## Interaction matrix
    mats2 <- fill_mat(mats2, dis = "beta", p1 = beta1, p2 = beta2)
    diag(mats2) <- -rbeta(nrow(mats), d1p, d2p)
    
    # Dynamics
    
    ## Initial Abundance
    inab[,i] <- runif(nrow(mats),.001,.02)
    
    ## Growth rates
    gr <- runif(nrow(mats), 0, Rmax)
    
    parlst[[i]] <- list(alpha = gr, m = mats2, K = 20)
    dyn <- ode(inab[,i], times = 1:tmax, func = lvmod, parms = parlst[[i]], events = list(func = ext1, time =  1:tmax))
    
    if(nrow(dyn) == tmax){
      spp <- dyn[tmax, -1] > 0
      mcon <- graph.adjacency(abs(sign(parlst[[i]]$m[spp, spp])))
      cond <- ifelse(is.connected(mcon), TRUE, FALSE)
    }else{
      cond <- FALSE
    }
  }
  
  matplot(dyn[,-1], typ = "l", main = i)
  
  spp <- dyn[tmax, -1] > 0
  dfin <- dyn[tmax, -1]
  #grfin <- parlst$alpha
  #D <- diag(parlst$m)
  
  #df <- data.frame(spp, ab = dfin, ab.i = inab, gr = grfin, d = D, K = Kcomm, K2 = -gr/D)
  
  #sprF <- sp_role(parlst$m[spp,spp])
  #spiF <- int_sp(parlst$m[spp,spp])
  #sprI <- sp_role(parlst$m)
  #spiI <- sp_role(parlst$m)
  tga <- apply(dyn[,-1], 1, function(x) get_abundvec(x[x > 0], 2000))
  tfz[[i]] <- do.call(rbind, lapply(tga, fzmod))
  
  conn <- sum(parlst$m[spp, spp] != 0)/(sum(spp)*sum(spp))
  
  #fin.df <- data.frame(df[spp,], sprF, spiF)
  #init.df <- data.frame(df, sprI, spiI)
  
  reslst2[[i]] <- list(dfin, conn)
  #reslst[[i]] <- list(init.df, fin.df, conn)
  #return(reslst)
  print(i)
}

tga <- apply(dyn[,-1], 1, function(x) get_abundvec(x, 2000))
tfz <- do.call(rbind, lapply(tga, fzmod))

plot(log10(allf$s)~log10(allf$N), ylim = c(-.5,1))
#ga <- lapply(lapply(reslst, "[[", 2), get_abundvec)
#do.call(rbind,lapply(ga, fzmod))
rb1 <- do.call(rbind, lapply(lapply(lapply(reslst2, "[[",1), get_abundvec, N = 2000), fzmod))
points(rb1[,1:2], col = "green4", pch = 20)

points(do.call(rbind ,sapply(tfz, function(x) x[500:2000,1:2])), pch = 20, col = "green4")
points(do.call(rbind ,sapply(tfz, function(x) x[2000,1:2])), pch = 20, col = "blue")

reslst3 <- list()
tfz2 <- list()
for(i in 1:length(parlst)){
  dyn <- ode(inab[,i], times = 1:tmax, func = lvmodK, parms = parlst[[i]], events = list(func = ext1, time =  1:tmax))
  if(nrow(dyn) != 2000){
    reslst3[[i]] <- list(NA, NA)
    tfz2[[i]] <- data.frame(N = NA, s = NA, nll = NA, r2 = NA)
    next
  }else{
    spp <- dyn[tmax, -1] > 0
    dfin <- dyn[tmax, -1]
    conn <- sum(parlst$m[spp, spp] != 0)/(sum(spp)*sum(spp))
    reslst3[[i]] <- list(dfin, conn)
  }
  tga <- apply(dyn[,-1], 1, function(x) get_abundvec(x[x > 0], 2000))
  tfz2[[i]] <- do.call(rbind, lapply(tga, fzmod))
  
  print(i)
}

rb2 <- do.call(rbind, lapply(lapply(lapply(reslst3, "[[",1), get_abundvec, N = 2000), fzmod))
points(rb2[,1:2], col = "darkgrey", pch = 20)


points(do.call(rbind, lapply(tfz2, function(x) x[500:2000,1:2])), pch = 20, col = "purple4")
points(t(sapply(tfz2, function(x) x[2000,1:2])), pch = 20, col = "gray")


Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
}  

plot()
Plot_ConvexHull(allf$N, allf$s, "black")
pt1 <- log10(do.call(rbind ,sapply(tfz, function(x) x[500:2000,1:2])))
Plot_ConvexHull(pt1[,1], pt1[,2], "blue")
pt2 <- log10(do.call(rbind, lapply(tfz2, function(x) x[500:2000,1:2])))
Plot_ConvexHull(pt2[,1], pt2[,2], "gray")
