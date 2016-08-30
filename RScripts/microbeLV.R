library(deSolve)



mouse <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/murineMICROBEint.csv", row.names = 1)
mouseSE <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/murineMICROBEintSE.csv", row.names = 1)
alphas <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/murineMICROBEa.csv")

mse <- as.matrix(mouseSE)
mse[is.na(mse)] <- 0

m <- as.matrix(mouse)

lvmod <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha * (1 - state/1) + state * parms$m %*% state
    
    list(dB)
  })
}

ext1 <- function (times, states, parms){
  with(as.list(states), {
    states[states < 10^-30] <- 0 
    return(c(states))
  })
}


parms <- list(alpha = (alphas$Alpha), m = m)

system.time(
res1 <- ode(runif(17), 1:1000, parms = parms, func = lvmod, events = list(func = ext1, time =  1:1000))
)
res1
matplot(res1[,-1], typ = "l", lwd = 2)


B <- matrix(runif(17 * 1000, .5, 1), nrow = 17, ncol = 1000)

strt <- Sys.time()
out <- matrix(nrow = 1000, ncol = 17)
tte <- matrix(nrow = 1000, ncol = 17)
for(i in 1:1000){
  res <- ode(B[,i], 1:1000, parms = parms, func = lvmod, events = list(func = ext1, time = 1:1000))
  tte[i,] <- apply(res[,-1], 2, function(x) 1000 - sum(x == 0))
  out[i,] <- res[1000, -1]
  print(i)
}
end <- Sys.time()
end - strt

persist <- apply(out, 1, function(x){sum(x > 0)})
hist(persist)
p2 <- persist/17
fit1 <- glm(p2~t(B), family = "quasibinomial")
summary(fit1)
plot(fit1$coefficients[-1], apply(out, 2, median))


m2 <- m
m2[which(abs(m) <= 1)] <- 0
parms2 <- list(alpha = (alphas$Alpha)[-1], m = m[-1,-1])

out2 <- matrix(nrow = 100, ncol = 16)
for(i in 1:100){
  out2[i,] <- ode(B[-1,i], 1:1000, parms = parms2, func = lvmod, events = list(func = ext1, time = 1:1000))[1000, -1]
  print(i)
}

apply(out2, 1, function(x){sum(x > 0)})
plot(apply(out2, 1, function(x){sum(x > 0)})~persist[1:100])

dropoutSP <- function(x, iter = 100, t = 1000, B, par, mod = lvmod, ev = ext1){
  par$alpha <- par$alpha[-x]
  par$m <- par$m[-x, -x]
  
  out2 <- matrix(nrow = iter, ncol = nrow(par$m))
  for(j in 1:iter){
    out2[j,] <- ode(B[-x,j], 1:t, parms = par, func = mod, events = list(func = ev, time = 1:t))[t, -1]
  }
  return(out2)
}


dropoutSP2 <- function(x, iter = 100, t = 1000, B, par, mod = lvmod, ev = ext1){
  B[x,] <- 0
  out2 <- list()#matrix(nrow = iter, ncol = nrow(par$m))
  for(j in 1:iter){
    out2[[j]] <- ode(B[,j], 1:t, parms = par, func = mod, events = list(func = ev, time = 1:t))[-1,]
  }
  return(out2)
}


dropoutINT <-  function(x, iter = 100, t = 1000, B, par, mod = lvmod, ev = ext1){
  par$m[m !=0][x] <- 0
  
  out2 <- matrix(nrow = iter, ncol = nrow(par$m))
  for(j in 1:iter){
    out2[j,] <- ode(B[,j], 1:t, parms = par, func = mod, events = list(func = ev, time = 1:t))[t, -1]
    cat(j, " ")
  }
  return(out2)
}


B <- matrix(runif(17 * 1000, .5, 1), nrow = 17, ncol = 1000)
parms <- list(alpha = (alphas$Alpha), m = m)


p1mu <- c()
p1med <- c()
temp <- list()
p1 <- list()
for(i in 1:17){
  temp[[i]] <- dropoutSP(x = i, iter = 500, t = 1000, B = B, par = parms)
  p1[[i]] <- apply(temp[[i]], 1, function(x) sum(x > 0))
  
  p1mu[i] <- mean(p1[[i]])
  p1med[i] <- median(p1[[i]])
  
  print(i)
}


presences <- lapply(lapply(temp, function(x) unlist(apply(x, 1, function(j) which(j != 0)))), table)
1-p1mu/17

nz <- t(sapply(temp2, function(x){apply(x, 2, function(j) sum(j != 0))}))
diag(nz) <- NA
boxplot(nz/200)

strt <- Sys.time()
#p2mu <- c()
#p2med <- c()
temp2 <- lapply(1:17, function(x) lapply(1:500, function(y) matrix(NA, nrow = 1000, ncol = 17)))
#p2 <- list()
for(i in 1:17){
  temp2[[i]] <- dropoutSP2(x = i, iter = 500, t = 1000, B = B, par = parms)
  #p2[[i]] <- apply(temp2[[i]], 1, function(x) sum(x > 0))
  
  #p2mu[i] <- mean(p2[[i]])
  #p2med[i] <- median(p2[[i]])
  
  print(i)
}
end <- Sys.time()
end-strt

length(temp2[[1]])
#get time to extinction
lapply(temp2, function(x){sapply(x, function(y){apply(y[,-1], 2, function(x) 1000 - sum(x == 0))})})
#get equilibrium comms
lapply(temp2, function(x){sapply(x, function(y){y[1000,-1]})})



library(rootSolve)
eig <- c()
for(i in 1:500){
  jac <- jacobian.full(B[,i], lvmod, parms = parms)
  eig[i]<- max(Re(eigen(jac)$values))
  
}

hist(eig)
parms$m







# Stochastic

t.max <- 1000
B1 <- matrix(runif(17), nrow = t.max, ncol = 17, byrow = T)

for(i in 2:t.max){
  m.sto <- matrix(rnorm(length(m), m, mse), 17, 17)
  parms <- list(alpha = (alphas$Alpha), m = m.sto)
  
  res1 <- ode(B1[i-1,], 1:2, parms = parms, func = lvmod, events = list(func = ext1, time =  1:2))
  B1[i,] <- res1[2,-1]
}

matplot(B1, typ = "l")
B1[t.max,]


stodyn <- function(t.max, ints = m, intSE = mse, a = alphas$Alpha){
  
  B1 <- matrix(runif(17), nrow = t.max, ncol = 17, byrow = T)
  
  for(i in 2:t.max){
    m.sto <- matrix(rnorm(length(ints), ints, intSE), 17, 17)
    parms <- list(alpha = (a), m = m.sto)
    
    res1 <- ode(B1[i-1,], 1:2, parms = parms, func = lvmod, events = list(func = ext1, time =  1:2))
    B1[i,] <- res1[2,-1]
  }
  return(B1)
}

matplot(stodyn(100), typ = "l")

test <- list()
for(i in 1:100){
  test[[i]] <- stodyn(100)
}

test2 <- lapply(1:17, function(x) sapply(test, function(q) q[,x]))
matplot(sapply(test2, function(x) apply(x, 1, median)))


xi <- runif(17)

mre <- c()
for(i in 1:1000){
  xi <- B[,i]
  mre[i] <- max(Re(eigen(jacobian.full(xi, func = lvmod, parms = parms))$values))
}
plot(mre)

out1 <- ode(xi, 1:1000, parms = parms, func = lvmod, events = list(func = ext1, time = 1:1000))
matplot(out1[,-1], typ = "l")
