#### Effect of truncation on s value
####################################################################################################################################
####################################################################################################################################
#### Functions
resamp <- function(X){
  g1 <- apply(otu3, 2, get_abundvec, N = X)
  g2 <- lapply(psd2$eqa, get_abundvec, X)
  g3 <- lapply(psd4$eqa, function(x) get_abundvec(x[x>0], N = X))
  
  return(list(g1, g2, g3))
}

get_trunc_s <- function(rs, trunc = seq(10, 400, 10)){
  strt <- Sys.time()
  hmp1 <- rs[[1]]
  gavsim <- rs[[2]]
  gavsim4 <- rs[[3]]
  hmp.sran <- list()
  sim.sran1 <- list()
  sim.sran2 <- list()
  for(i in 1:length(trunc)){
    hmp.sran[[i]] <- sapply(hmp1[sapply(hmp1, length) >= trunc[i]], function(x) fitzipf_r(x[1:trunc[i]])@coef)
    sim.sran1[[i]] <- sapply(gavsim[sapply(gavsim, length) >= trunc[i]], function(x) fitzipf_r(x[1:trunc[i]])@coef)
    sim.sran2[[i]] <- sapply(gavsim4[sapply(gavsim4, length) >= trunc[i]], function(x) fitzipf_r(x[1:trunc[i]])@coef)
  }
  fin <- Sys.time()
  cat(fin-strt)
  
  return(list(hmp.sran, sim.sran1, sim.sran2))
}

plot_trunc <- function(tlists, X = 2000){
  hmp.sran <- tlists[[1]]
  sim.sran1 <- tlists[[2]]
  sim.sran2 <- tlists[[3]]
  
  plot(sapply(hmp.sran, median)~trunc, ylim = c(0,4), typ = "l", lwd = 2, main = X, xlab = "Truncation", ylab = "Parameter Value (s)")
  points(sapply(hmp.sran, min)~trunc, typ = "l", lwd = 2)
  points(sapply(hmp.sran, max)~trunc, typ = "l", lwd = 2)
  
  points(sapply(sim.sran1, median)~trunc, typ = "l", lty = 2, col = "blue", lwd = 2)
  points(sapply(sim.sran1, min)~trunc, typ = "l", lty = 2, col = "blue", lwd = 2)
  points(sapply(sim.sran1, max)~trunc, typ = "l", lty = 2, col = "blue", lwd = 2)
  
  points(sapply(sim.sran2, median)~trunc, typ = "l", lty = 3, col = "darkgreen", lwd = 2)
  points(sapply(sim.sran2, min)~trunc, typ = "l", lty = 3, col = "darkgreen", lwd = 2)
  points(sapply(sim.sran2, max)~trunc, typ = "l", lty = 3, col = "darkgreen", lwd = 2)
  
}


####################################################################################################################################
####################################################################################################################################
#### Simulation

t1 <- Sys.time()
test <- resamp(2000)
srans <- get_trunc_s(test)
plot_trunc(srans)
t2 <- Sys.time()
t2 - t1
# 22 min

t1 <- Sys.time()
test2 <- resamp(5000)
srans2 <- get_trunc_s(test2)
plot_trunc(srans2)
t2 <- Sys.time()
t2 - t1
# 23.5 min

t1 <- Sys.time()
test3 <- resamp(10000)
srans3 <- get_trunc_s(test3)
plot_trunc(srans3)
t2 <- Sys.time()
t2 - t1
# 25 min

t1 <- Sys.time()
test4 <- resamp(15000)
srans4 <- get_trunc_s(test4)
plot_trunc(srans3)
t2 <- Sys.time()
t2 - t1
# 25 min

####################################################################################################################################
####################################################################################################################################
#### Plot

par(mfrow = c(1,4))
plot_trunc(srans, X = 2000)
plot_trunc(srans2, X = 5000)
plot_trunc(srans3, X = 10000)
plot_trunc(srans4, X = 15000)