get_abundvec <- function(abund, N = 10000){
  r.ab <- abund/sum(abund)
  samp2 <- sample(1:length(abund), N, replace = T, prob = r.ab)
  t1 <- table(samp2)
  return(as.numeric(t1))
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

fzmod <- function(x){
  fz1 <- fitzipf_r(x)
  fc1 <- fz1@fullcoef
  nll <- fz1@minuslogl(N = fz1@fullcoef[1], s = fz1@fullcoef[2])
  r2 <- r2modified(sort(x, decreasing = T), radpred(fz1)$abund)
  
  return(data.frame(t(fc1), nll, r2))
}

r2modified <- function(x,y,log=FALSE){
  if(log){
    rm = 1-sum((log10(x)-log10(y))^2)/sum((log10(x)-mean(log10(x)))^2)
  }
  else{
    rm = 1-sum((x-y)^2)/sum((x-mean(x))^2)
  }
  return(rm)  
}

library(deSolve)

lvmodK <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha * (1 - state/(parms$K/sum(state > 10^-10))) + state * parms$m %*% state 
    #dB <- state * parms$alpha * (1 - state/(parms$K)) + state * parms$m %*% state 
    
    list(dB)
  })
}

lvmodKi <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha * (1 - state/(parms$K)) + state * parms$m %*% state 
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

######################################################################################
######################################################################################
######################################################################################

multityp <- lapply(1:20, function(x){
  S <- 600
  p1 <- runif(1,0,1)
  p2 <- runif(1, p1, 1)
  c1 <- runif(1, .1, .3)
  mats <- get.adjacency(erdos.renyi.game(S, c1, "gnp", directed = F), sparse = F)
  tat <- mats*sample(c(-1,1,0), length(mats), replace = T, prob = c(p1,p2-p1,1-(p2)))
  return((tat))
})

multityp.fill <- lapply(multityp, function(x) fill_mat(x, dis = "beta", p1 = 1, p2 = 2))

par1 <- list()
dyn <- list()
for(i in 1:20){

  diag(multityp.fill[[i]]) <- (runif(nrow(multityp.fill[[i]]), -1, 0))

  par1[[i]] <- list(alpha = runif(nrow(multityp.fill[[i]]), 0,.1), m = multityp.fill[[i]], K = 20)
  dyn[[i]] <-(ode(runif(nrow(multityp.fill[[i]]),.001,.02), times = 1:2000, func = lvmodK, parms = par1[[i]], events = list(func = ext1, time =  1:2000)))
  
  matplot(dyn[[i]][,-1], typ = "l", main = i)
}


eqab <- lapply(dyn, function(x) if(nrow(x) == 2000){x[2000,-1][x[2000,-1]>0]}else{NA})
sapply(eqab, function(x) if(any(is.na(x))){NA}else{fzmod((x))})[[6]]


p1 <- par1[[6]]
d1 <- dyn[[6]]


par3 <- list()
dyn2 <- list()
#tops <- c(545, 205, 116, 587, 191)
for(i in 1:100){
  par2 <- p1
  #par2$m[par2$m != 0] <- abs(rnorm(sum(p1$m != 0), p1$m[p1$m != 0], abs(p1$m[p1$m != 0]*.1)))*sign(p1$m[p1$m != 0])
  #diag(par2$m) <- diag(p1$m)
  par2$K <- rep(p1$K/sum(d1[2000,-1] > 10^-10), 600)
  #par2$alpha <- rev(sort(par2$alpha))
  par2$K[7] <- par2$K[7]*(2) 
  #r1 <- rlnorm(600, -1, 1)
  #par2$K <- par2$K[order(d1[2000,-1], decreasing = T)]*sort(r1, decreasing = T)
  
  par3[[i]] <- par2
  dyn2[[i]] <-(ode(d1[1,-1], times = 1:2000, func = lvmodKi, parms = par3[[i]], events = list(func = ext1, time =  1:2000)))
  
  matplot(dyn2[[i]][,-1], typ = "l", main = i)
}


eqab <- lapply(dyn2, function(x) if(nrow(x) == 2000){x[2000,-1][x[2000,-1]>0]}else{NA})
sapply(eqab, function(x) if(any(is.na(x))){NA}else{fzmod((x))})[,1]
mean(unlist(sapply(eqab, function(x) if(any(is.na(x))){NA}else{fzmod((x))})[1,]))

par4 <- list()
dyn3 <- list()
#tops <- c(545, 205, 116, 587, 191)
for(i in 1:10){
  par2 <- p1
  par2$m[par2$m != 0] <- abs(rnorm(sum(p1$m != 0), p1$m[p1$m != 0], abs(p1$m[p1$m != 0]*.1)))*sign(p1$m[p1$m != 0])
  diag(par2$m) <- diag(p1$m)
  #par2$K <- rep(p1$K/sum(d1[2000,-1] > 10^-10), 600)
  #par2$alpha <- rev(sort(par2$alpha))
  #par2$K[545] <- par2$K[545]*(i/2) 
  #r1 <- rlnorm(600, -1, 1)
  #par2$K <- par2$K[order(d1[2000,-1], decreasing = T)]*sort(r1, decreasing = T)
  
  par4[[i]] <- par2
  dyn3[[i]] <-(ode(d1[1,-1], times = 1:2000, func = lvmodK, parms = par4[[i]], events = list(func = ext1, time =  1:2000)))
  
  matplot(dyn3[[i]][,-1], typ = "l", main = i)
}



eqab <- lapply(dyn3, function(x) if(nrow(x) == 2000){x[2000,-1][x[2000,-1]>0]}else{NA})
sapply(eqab, function(x) if(any(is.na(x))){NA}else{fzmod((x))})
mean(unlist(sapply(eqab, function(x) if(any(is.na(x))){NA}else{fzmod((x))})[1,]))
