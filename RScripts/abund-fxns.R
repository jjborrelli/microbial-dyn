
############################################################################
############################################################################
#### Libraries
library(deSolve)
library(igraph)

############################################################################
############################################################################
#### Dynamics Functions

lvmodK <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha * (1 - state/(parms$K/sum(state > 10^-10))) + state * parms$m %*% state 
    #dB <- state * parms$alpha * (1 - state/(parms$K)) + state * parms$m %*% state 
    
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
############################################################################
#### Setup Functions

fill_mat <- function(mat, dis, p1 = .5, p2 = 1){
  t1 <- mat
  diag(t1) <- 0  
  
  if(dis == "beta"){
    t1[t1 != 0] <- rbeta(sum(t1 != 0), p1, p2) * sign(mat[mat != 0])
  }else if(dis == "unif"){
    t1[t1 == 1] <- runif(sum(t1 == 1), 0, p1)
    t1[t1 == -1] <- runif(sum(t1 == -1), -p2, 0) 
  }else if(dis == "norm"){
    t1[t1 == 1] <- abs(rnorm(sum(t1 == 1), 0, p1))
    t1[t1 == -1] <- -abs(rnorm(sum(t1 == -1), 0, p2))
  }else{
    t1[t1 == 1] <- runif(sum(t1 == 1), 0, 1)
    t1[t1 == -1] <- runif(sum(t1 == -1), -1, 0)
  }
  
  return(t1)
}


############################################################################
############################################################################
#### Data Functions

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

# info on species level interactions
int_sp <- function(mat){
  d1 <- diag(mat)
  diag(mat) <- 0
  mm1 <- matrix(nrow = nrow(mat), ncol = 18)
  colnames(mm1) <- c("self", "nComp", "CompIn", "CompOut", "nMut", "MutIn", "MutOut", "nPred", "PredIn", "PredOut", "nAmens", "AmensIn", "AmensOut", "nComm", "CommIn", "CommOut", "allIn", "allOut")
  for(i in 1:nrow(mat)){
    i1 <- mat[,i]
    i2 <- mat[i,]
    
    comp <- which(i1 < 0 & i2 < 0)
    mut <- which(i1 > 0 & i2 > 0)
    pred <- which(i1 > 0 & i2 < 0 | i1 < 0 & i2 > 0)
    amens <- which(i1 < 0 & i2  == 0 | i1 == 0 & i2 < 0)
    comm <- which(i1 > 0 & i2  == 0 | i1 == 0 & i2 > 0)
    
    mm1[i, "nComp"] <- length(comp)
    mm1[i, "nMut"] <- length(mut)
    mm1[i, "nPred"] <- length(pred)
    mm1[i, "nAmens"] <- length(amens)
    mm1[i, "nComm"] <- length(comm)
    
    mm1[i, "CompIn"] <- ifelse(is.na(mean(i1[comp])), 0, mean(i1[comp]))
    mm1[i, "MutIn"] <- ifelse(is.na(mean(i1[mut])), 0, mean(i1[mut]))
    mm1[i, "PredIn"] <- ifelse(is.na(mean(i1[pred])), 0, mean(i1[pred]))
    mm1[i, "AmensIn"] <- ifelse(is.na(mean(i1[amens])), 0, mean(i1[amens]))
    mm1[i, "CommIn"] <- ifelse(is.na(mean(i1[comm])), 0, mean(i1[comm]))
    
    mm1[i, "CompOut"] <- ifelse(is.na(mean(i2[comp])), 0, mean(i2[comp]))
    mm1[i, "MutOut"] <- ifelse(is.na(mean(i2[mut])), 0, mean(i2[mut]))
    mm1[i, "PredOut"] <- ifelse(is.na(mean(i2[pred])), 0, mean(i2[pred]))
    mm1[i, "AmensOut"] <- ifelse(is.na(mean(i2[amens])), 0, mean(i2[amens]))
    mm1[i, "CommOut"] <- ifelse(is.na(mean(i2[comm])), 0, mean(i2[comm]))
    
    mm1[i, "self"] <- d1[i]
    mm1[i, "allIn"] <- ifelse(sum(i1 != 0) > 0, mean(i1[i1 != 0]), 0)
    mm1[i, "allOut"] <- ifelse(sum(i2 != 0) > 0, mean(i2[i2 != 0]), 0)
  }
  
  
  return(mm1)
}

# network level data

sp_role <- function(mat){
  # convert to unweighted and weighted graphs 
  g <- graph.adjacency(abs(sign(mat)))
  g2 <- graph.adjacency(abs(mat), weighted = T)
  # unweighted and weighted betweenness
  b.uw <- igraph::betweenness(g)
  b.w <- igraph::betweenness(g2)
  # degree
  d.in <- igraph::degree(g, mode = "in")
  d.out <- igraph::degree(g, mode = "out")
  d.tot <- igraph::degree(g, mode = "total")
  # clustering
  cc.uw <- igraph::transitivity(g2, "local")
  cc.w <- igraph::transitivity(g2, "weighted")
  # average path length
  apl.uw.mean <- colMeans(igraph::distances(g))
  #apl.uw.median <- apply(distances(g), 2, median)
  apl.w.mean <- colMeans(igraph::distances(g2))
  #apl.w.median <- apply(distances(g2), 2, median)
  wc.uw <- walktrap.community(g)
  wc.mod.uw <- wc.uw$modularity
  wc.mem.uw <- wc.uw$membership
  wc.w <- walktrap.community(g2)
  wc.mod.w <- wc.w$modularity
  wc.mem.w <- wc.w$membership
  
  #res <- matrix(c(b.uw, b.w, d.in, d.out, d.tot, cc.uw, cc.w, apl.uw.mean, apl.uw.median, apl.w.mean, apl.w.median),
  #              nrow = nrow(mat), ncol = 11)
  res <- matrix(c(b.uw, b.w, d.in, d.out, d.tot, cc.uw, cc.w, apl.uw.mean, apl.w.mean,
                  wc.mod.uw, wc.mem.uw, wc.mod.w, wc.mem.w),
                nrow = nrow(mat), ncol = 13)
  colnames(res) <- c("bet.uw", "bet.w", "d.in", "d.out", "d.tot", "cc.uw", "cc.w", 
                     "apl.uw.mu", "apl.w.mu", "mod.uw", "mem.uw", "mod.w", "mem.w")
  
  #res <- cbind(res, mod)
  return(res)
}

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

webprops <- function(mat){
  # convert to unweighted and weighted graphs 
  g <- graph.adjacency(abs(sign(mat)))
  g2 <- graph.adjacency(abs(mat), weighted = T)
  
  cc.w <- transitivity(g, "global")
  
  conn <- edge_density(g, loops = T)
  
  diam.uw <- diameter(g)
  diam.w <- diameter(g2)
  
  mod.uw <- modularity(g, walktrap.community(g)$membership)
  mod.w <- modularity(g2, walktrap.community(g2)$membership)
  
  apl <- average.path.length(g)
  nclust <- no.clusters(g2, "strong")
  
  m.int <- mean(mat[mat != 0])
  m.n.int <- mean(mat[mat < 0])
  m.p.int <- mean(mat[mat > 0])
  
  #itypes
  ity <- itypes(mat)
  cnames <- names(ity)
  ity <- data.frame(matrix(ity, ncol = 5))
  colnames(ity) <- cnames
  return(data.frame(cc.w, conn, diam.uw, diam.w, apl, mod.uw, mod.w, nclust, m.int, m.n.int, m.p.int, ity))
}

localstability <- function(tp, par, dyn){
  mrev <- c()
  for(i in 1:length(tp)){
    t1 <- tp[i]
    jf <- rootSolve::jacobian.full(dyn[t1,dyn[t1,] > 0], func = lvmod, parms = par)
    mrev[i] <- max(Re(eigen(jf)$values))
  }
  
  return(mrev)
}



fitzipf_r <- function(x, N, trunc, start.value, upper = 20, ...){
  if (any(x <= 0)) stop ("All x must be positive")
  if(class(x)!="rad") rad.tab <- rad(x)
  else rad.tab <- x
  # y <- rep(rad.tab$rank, rad.tab$abund)
  dots <- list(...)
  if (!missing(trunc)){
    if (min(rad.tab$rank)<=trunc) stop("truncation point should be lower than the lowest rank") #
  }
  if(missing(N)){
    N <- max(rad.tab$rank) #
  }
  if(missing(start.value)){
    p <- rad.tab$abund/sum(rad.tab$abund)
    lzipf <- function(s, N) -s*log(1:N) - log(sum(1/(1:N)^s))
    opt.f <- function(s) sum((log(p) - lzipf(s, length(p)))^2)
    opt <- optimize(opt.f, c(0.5, length(p)))
    sss <- opt$minimum
  }
  else{
    sss <- start.value
  }
  if(missing(trunc)){
    LL <- function(N, s) -sum(rad.tab$abund*dzipf(rad.tab$rank, N=N, s=s, log = TRUE)) #
  }
  else{
    LL <- function(N, s) -sum(rad.tab$abund*dtrunc("zipf", x = rad.tab$rank, coef = list(N = N, s = s), trunc = trunc, log = TRUE)) #
  }
  result <- do.call("mle2", c(list(LL, start = list(s = sss), data = list(x = rad.tab$rank), fixed=list(N=N), method = "Brent", lower = 0, upper = upper), dots))
  if(abs(as.numeric(result@coef) - upper) < 0.001)
    warning("mle equal to upper bound provided. \n Try increase value for the 'upper' argument")
  new("fitrad", result, rad="zipf", distr = "zipf of relative abundance", trunc = ifelse(missing(trunc), NaN, trunc), rad.tab=rad.tab)
}


# Function to get interaction numbers and strengths for a given community
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


# My r2 function
get_r2 <- function(o, p){
  1 - sum((o-p)^2)/sum((o - mean(o))^2)
}


# Yingnan's modified r2 function
r2modified <- function(x,y,log=FALSE){
  if(log){
    rm = 1-sum((log10(x)-log10(y))^2)/sum((log10(x)-mean(log10(x)))^2)
  }
  else{
    rm = 1-sum((x-y)^2)/sum((x-mean(x))^2)
  }
  return(rm)  
}


get_abundvec <- function(abund, N = 10000){
  r.ab <- abund/sum(abund)
  samp2 <- sample(1:length(abund), N, replace = T, prob = r.ab)
  t1 <- table(samp2)
  return(as.numeric(t1))
}


fzmod <- function(x){
  fz1 <- fitzipf_r(x)
  fc1 <- fz1@fullcoef
  nll <- fz1@minuslogl(N = fz1@fullcoef[1], s = fz1@fullcoef[2])
  r2 <- r2modified(sort(x, decreasing = T), radpred(fz1)$abund)
  
  return(data.frame(t(fc1), nll, r2))
}

