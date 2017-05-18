# using methods from Berry and Widder 2014

lv_bw <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha * (1 + (parms$m %*% state)/parms$K)
    
    list(dB)
  })
}


multityp <- lapply(1:1000, function(x){
  S <- 500
  p1 <- runif(1,0,1)
  p2 <- runif(1, p1, 1)
  c1 <- runif(1, .1, .2)
  mats <- get.adjacency(erdos.renyi.game(S, c1, "gnp", directed = T), sparse = F)
  tat <- mats*sample(c(-1,1,0), length(mats), replace = T, prob = c(p1,p2-p1,1-(p2)))
  return((tat))
})

multityp.fill <- fill_mats(multityp, sdevn = -1, sdevp = 1)

s1 <- Sys.time()
dfin <- list()
dfin2 <- list()
dfin3 <- list()
dfin4 <- list()
cv <- c()
k2 <- list()
for(i in 1:1000){
  
  diag(multityp.fill[[i]]) <- -1#(runif(nrow(multityp.fill[[i]]), -2, -1))
  
  k1 <- rlnorm(nrow(multityp.fill[[i]]), 0, 1)
  k2[[i]] <- k1/max(k1)*100
  par1 <- list(alpha = runif(nrow(multityp.fill[[i]]), 0, 1), m = multityp.fill[[i]], K = k2[[i]])
  dyn <-(ode(runif(nrow(multityp.fill[[i]]),1,10), times = 1:1000, func = lv_bw, parms = par1, events = list(func = ext1, time =  1:1000)))
  
  matplot(dyn[,-1], typ = "l", main = i)
  if(nrow(dyn) == 1000){
    dfin[[i]] <- dyn[1000,-1]
    dfin2[[i]] <- dyn[10,-1]
    dfin3[[i]] <- dyn[50,-1]
    dfin4[[i]] <- dyn[100,-1]
    cv[[i]] <- sd(apply(dyn[900:1000,-1], 1, function(x) mean(x[x>0])))/mean(apply(dyn[900:1000,-1], 1, function(x) mean(x[x>0])))
  }else{
    dfin[[i]] <- NA
    dfin2[[i]] <- NA
    dfin3[[i]] <- NA
    dfin4[[i]] <- NA
    cv[i] <- NA
  }
  
}
s2 <- Sys.time()
s2-s1

t2k <- rad_info(dfin, multityp.fill)
t10 <- rad_info(dfin2, multityp.fill)
t50 <- rad_info(dfin3, multityp.fill)
t100 <- rad_info(dfin4, multityp.fill)
pairs(cbind(t2k$sV, t10$sV, t50$sV, t100$sV))

median(t2k$sV)
median(t100$sV)
median(t50$sV)
median(t10$sV)