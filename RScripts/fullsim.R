library(parallel)
library(doSNOW)


filepath1 <- "~/Abundance/Data/"
iter <- 1:2
cl <- makeCluster(detectCores() - 1)


clusterExport(cl, c("filepath1"))
registerDoSNOW(cl)

RES <- foreach(x = iter) %do% {
  # Load necessary functions and libraries
  source("~/Desktop/GitHub/microbial-dyn/RScripts/abund-fxns.R")
  
  # Parameters
  
  S <- 600        # Number of Species
  beta1 <- 1      # first beta param - interactions
  beta2 <- 4      # second beta param - interactions
  d1p <- 2        # first beta param - diagonal
  d2p <- 1        # second beta param - diagonal
  Rmax <- .1      # maximum growth rate for unif dist
  Kcomm <- 20     # community carrying capacity
  tmax <- 2000    # max number of timesteps
  
  cond <- FALSE
  while(!cond){
    # Create community intereaction network
    
    ## Probability of negative, positive, and no interaction
    p1 <- runif(1,0,1)
    p2 <- runif(1, p1, 1)
    
    ## Connectance
    c1 <- runif(1, .1, .8)
    
    ## Adjacency matrix
    mats <- get.adjacency(erdos.renyi.game(S, c1, "gnp", directed = F), sparse = F)
    mats <- mats*sample(c(-1,1,0), length(mats), replace = T, prob = c(p1,p2-p1,1-(p2)))
    
    ## Interaction matrix
    mats <- fill_mat(mats, dis = "beta", p1 = beta1, p2 = beta2)
    diag(mats) <- -rbeta(nrow(mats), d1p, d2p)
    
    # Dynamics
    
    ## Initial Abundance
    inab <- runif(nrow(mats),.001,.02)
    
    ## Growth rates
    gr <- runif(nrow(mats), 0, Rmax)
    
    parlst <- list(alpha = gr, m = mats, K = 20)
    dyn <- ode(inab, times = 1:tmax, func = lvmodK, parms = parlst, events = list(func = ext1, time =  1:tmax))
    
    if(nrow(dyn) == tmax){
      spp <- dyn[tmax, -1] > 0
      mcon <- graph.adjacency(abs(sign(parlst$m[spp, spp])))
      cond <- ifelse(is.connected(mcon), TRUE, FALSE)
    }else{
      cond <- FALSE
    }
  }
  
  #matplot(dyn[,-1], typ = "l", main = 1)
  
  spp <- dyn[tmax, -1] > 0
  dfin <- dyn[tmax, -1]
  grfin <- parlst$alpha
  D <- diag(parlst$m)
  
  df <- data.frame(spp, ab = dfin, ab.i = inab, gr = grfin, d = D, K = Kcomm, K2 = -gr/D)
  
  sprF <- sp_role(parlst$m[spp,spp])
  spiF <- int_sp(parlst$m[spp,spp])
  sprI <- sp_role(parlst$m)
  spiI <- sp_role(parlst$m)
  
  
  fin.df <- data.frame(df[spp,], sprF, spiF)
  init.df <- data.frame(df, sprI, spiI)
  
  reslst <- list(init.df, fin.df)
  return(reslst)
}

