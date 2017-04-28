library(igraph)
library(NetIndices)
library(deSolve)
library(ggplot2)
library(rootSolve)
library(reshape2)


#################################################################################################
#################################################################################################
#################################################################################################

# Lotka-Volterra model with evenly distributed K
lvmodK <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha * (1 - state/(parms$K/sum(state > 10^-10))) + state * parms$m %*% state 
    #dB <- state * parms$alpha * (1 - state/(parms$K)) + state * parms$m %*% state 
    
    list(dB)
  })
}


# Lotka-Volterra model with set K
lvmodK2 <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha + state * parms$m %*% state 
    #dB <- state * parms$alpha * (1 - state/(parms$K)) + state * parms$m %*% state 
    
    list(dB)
  })
}

lvmodN <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    imat <- parms$m
    for(i in 1:nrow(parms$m)){
      for(j in 1:ncol(parms$m)){
        if(i == j){next} 
        if(state[i] < 0){state[i] <- 0}
        if(state[j] < 0){state[j] <- 0}
        p1 <- (state[i]/sum(state))*(state[j]/sum(state)); imat[i,j] <- imat[i,j] * sample(c(1,0), 1, prob = c(p1, 1-p1))
      }
      }
    dB <- state * parms$alpha * (1 - state/(parms$K/sum(state > 0))) + state * imat %*% state 
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

# Function to detect extinction and add immigration
ext2 <- function(times, states, parms){
  with(as.list(states), {
    states[states < 10^-5] <- 0 
    #states <- states + (runif(length(states), 10^-5, 10^-2)*sample(c(0,1), length(states), prob = c(.9,.1), replace = T))
    states <- abs(rnorm(1:length(states), states, .01))
    return(c(states))
  })
}


# fill interaction matrices with strengths
fill_mats <- function(mats, sdevp = .5, sdevn = 1){
  t2 <- list()
  for(i in 1:length(mats)){
    t1 <- mats[[i]]
    diag(t1) <- 0  #-rbeta(length(diag(t1)), 1.1, 5)*5
    t1[t1 == 1] <- abs(rnorm(sum(t1 == 1), .1, sdevp))
    t1[t1 == -1] <- -abs(rnorm(sum(t1 == -1), .1, sdevn))
    t2[[i]] <- t1 
    
    #t1[t1 == 1] <- runif(sum(t1 == 1), 0, sdevp)
    #t1[t1 == -1] <- runif(sum(t1 == -1), sdevn, 0)
    #t2[[i]] <- t1     
    
        
  }
  return(t2)
}



#################################################################################################
#################################################################################################
#################################################################################################

rad_info <- function(dfin, multityp.fill, thres = 5){
  eqa <- lapply(dfin[!is.na(dfin)], function(x) sort(x[x>0], decreasing = T))
  eqa <- eqa[!sapply(eqa, function(x) length(x) == 0)]
  eqm <- lapply((1:length(dfin[!is.na(dfin)])), function(x) multityp.fill[[x]][dfin[!is.na(dfin)][[x]] > 0, dfin[!is.na(dfin)][[x]] > 0])
  eqm <- eqm[!sapply(eqm, function(x) any(is.na(x)))]
  ico <- sapply(eqm, function(x) is.connected(graph.adjacency(abs(sign(x)))))
  eqm1 <- lapply(eqm, itystr)
  eqm1 <- lapply(1:length(eqm1), function(x) data.frame(eqm1[[x]], N = nrow(eqm[[x]])))
  dstr <- sapply(eqm, function(x) mean(diag(x)))
  eqa2 <- lapply(eqa, get_abundvec, 2000)
  fzN <- t(sapply(eqa2[sapply(eqa2, length) > thres], fzmod))
  
  testdat <- data.frame(prepdat(eqa2[sapply(eqa2, length) > thres], eqm1[sapply(eqa2, length) > thres], fzN[,"s"], fzN[,"r2"], dstr[sapply(eqa2, length) > thres]), abs = sapply(eqa2, sum)[sapply(eqa2, length) > thres], ico = ico[sapply(eqa2, length) > thres])
  
  return(testdat)
}


multihub <- lapply(1:5, function(x){
  S = 500
  p1 <- runif(1,0,1)
  p2 <- runif(1, p1, 1)
  pow <- rbeta(1, 4, 1)
  mats <- get.adjacency(barabasi.game(S, pow, m = 50, directed = F), sparse = F)
  #tat <- tatoosh*sample(c(1,-1), length(tatoosh), replace = T, prob = c(p1,1-p1))
  tat <- mats*sample(c(-1,1,0), length(mats), replace = T, prob = c(p1,p2-p1,1-(p2)))
  return((tat))
})


multityp <- lapply(1:10, function(x){
  S <- 1000
  p1 <- runif(1,0,1)
  p2 <- runif(1, p1, 1)
  c1 <- runif(1, .1, .3)
  mats <- get.adjacency(erdos.renyi.game(S, c1, "gnp", directed = F), sparse = F)
  tat <- mats*sample(c(-1,1,0), length(mats), replace = T, prob = c(p1,p2-p1,1-(p2)))
  return((tat))
})

multityp.fill <- fill_mats(multityp, sdevn = 1, sdevp = 1)
multityp.fill <- multityp.fill[sapply(multityp.fill, function(x) is.connected(graph.adjacency(abs(sign(x)))))]
#multityp.fill <- append(multityp.fill, multityp.fill)
#multihub.fill <- fill_mats(multihub, sdevn = -2, sdevp = 1)
s1 <- Sys.time()
dfin <- list()
dfin2 <- list()
dfin3 <- list()
dfin4 <- list()
alphas <- list()
cv <- c()
par1 <- list()
dyn <- list()
for(i in 1:10){
  #if(i < 501){
  #  diag(multityp.fill[[i]]) <- -2#(runif(nrow(multityp.fill[[i]]), sample(c(-1,-2,-3),1), 0))
  #}else{
    diag(multityp.fill[[i]]) <- (runif(nrow(multityp.fill[[i]]), -1, 0))
  #}
  #diag(multityp.fill[[i]]) <- (runif(nrow(multityp.fill[[i]]), sample(c(-1,-2,-3),1), 0))
  par1[[i]] <- list(alpha = runif(nrow(multityp.fill[[i]]), 0,.1), m = multityp.fill[[i]], K = 20)
  dyn[[i]] <-(ode(runif(nrow(multityp.fill[[i]]),.0001,.002), times = 1:2000, func = lvmodK, parms = par1[[i]], events = list(func = ext1, time =  1:2000)))
  
  matplot(dyn[[i]][,-1], typ = "l", main = i)
}

s2 <- Sys.time()
s2-s1

mat <- multityp.fill[[19]]
diag(mat) <- 0
mat2 <- melt(mat)
mat2 <- mat2[mat2$value != 0,]
#d <- dyn[[19]]

eig <- c()
for(i in 1:2000){
  eig[i] <- max(Re(eigen(rootSolve::jacobian.full(d2[i,-1], func = lvmodK, parms = par1[[19]]))$values))
}
plot(eig, typ = "l")

strt <- dyn[[13]][2000,-1]
strt[7] <- 0
d2 <- (ode(strt, times = 1:500, func = lvmodK, parms = par1[[13]], events = list(func = ext1, time =  1:500)))
#matplot(d2[,-1], typ = "l")

### PLOT EXTINCTION
par(mar = c(5,5,3,3))
matplot(rbind(d[,-1], d2[,-1]), typ = "l", lwd = 3, xlab = "", ylab = "")
title(ylab = list("Abundance", cex = 2), xlab = list("Time", cex = 2))
points(x = 500, y = dyn[[13]][2000,-1][7], pch = 4, cex = 2, lwd = 2, col = "red")
text(x = 600, y = dyn[[13]][2000,-1][7], label = "<- Extinction", cex = 1.5)


strt <- dyn[[13]][2000,-1]
strt <- c(strt, .2)
par2 <- par1[[13]]
par2$alpha <- c(par2$alpha, .1)
par2$m <- rbind(cbind(par2$m, runif(10, -.1, .1)*sample(c(0,1), 10, replace = T)), c(runif(10, -.1,1)*sample(c(0,1), 10, replace = T), -1))
d3 <- (ode(strt, times = 1:500, func = lvmodK, parms = par2, events = list(func = ext1, time =  1:500)))
#matplot(d3[,-1], typ = "l")

### PLOT INVASION
par(mar = c(5,5,3,3))
matplot(rbind(cbind(d[,-1],NA), d3[,-1]), typ = "l", lwd = 3, xlab = "", ylab = "")
title(ylab = list("Abundance", cex = 2), xlab = list("Time", cex = 2))
points(x = 500, y = .2, pch = 10, cex = 2, lwd = 2)
text(x = 400, y = .2, label = "Invasion ->", cex = 1.5)

### PLOT 1000sp
par(mar = c(5,5,3,3))
matplot(dyn[[9]][,-1], typ = "l", lwd = 2, xlab = "", ylab = "")
title(ylab = list("Abundance", cex = 2), xlab = list("Time", cex = 2))


d.rem <- rbind(d[,-1], d2[,-1])
d.inv <- rbind(cbind(d[,-1],NA), d3[,-1])

mat <- par2$m
diag(mat) <- 0
mat2 <- reshape2::melt(mat)
mat2 <- mat2[mat2$value != 0,]
g <- graph.data.frame(mat2)
E(g)$weights <- abs(mat2$value)/max(abs(mat2$value))
E(g)$sign <- ifelse(sign(mat2$value) < 0, "red", "blue") 
V(g)$ab <- d[10,-1]*10

#gif.rem
g2 <- induced_subgraph(g, c(1:3, 5:11))
g3 <- induced_subgraph(g, c(1:3, 5,6,8:11))
kklay <- layout.kamada.kawai(g)
for(i in 1:1000){
  V(g2)$ab <- d.rem[i,as.numeric(names(V(g2)))]*100
  plot(g2, edge.color = E(g2)$sign, edge.curved = F, edge.width = E(g2)$weights*3, edge.arrow.size = 1.5, layout = kklay[1:10,], vertex.size = V(g2)$ab, vertex.color = "black", vertex.label.color = "blue")
  if(i >=500){points(x = kklay[6,1], y = kklay[6,2], pch = 10, cex = 2, lwd = 2)}
}

library(animation)
setwd("~/Desktop/gifs")
saveGIF(
  {
    for(i in seq(5,1000,5)){
     if(i < 501){
        V(g2)$ab <- d.rem[i,as.numeric(names(V(g2)))]*100
        plot(g2, edge.color = E(g2)$sign, edge.curved = F, edge.width = E(g2)$weights*3, edge.arrow.size = 1.5, layout = kklay[as.numeric(names(V(g2))),], vertex.size = V(g2)$ab, vertex.color = "black", vertex.label = NA)
      }else{
        V(g3)$ab <- d.inv[i,as.numeric(names(V(g3)))]*100
        
        plot(g3, edge.color = E(g3)$sign, edge.curved = F, edge.width = E(g3)$weights*3, edge.arrow.size = 1.5, 
             layout = kklay[as.numeric(names(V(g3))),], vertex.size = V(g3)$ab, vertex.color = "black", vertex.label = NA)
      }
    }
  }, 
  movie.name = "netdyn1a.gif", interval = 0.001, ani.width = 600, ani.height = 600, outdir = "~/Desktop/"
)

saveGIF(
  {
    for(i in seq(5,1000,5)){
      if(i < 501){
        V(g2)$ab <- d.rem[i,as.numeric(names(V(g2)))]*100
        plot(g2, edge.color = E(g2)$sign, edge.curved = F, edge.width = E(g2)$weights*3, edge.arrow.size = 1.5, layout = kklay[as.numeric(names(V(g2))),], vertex.size = V(g2)$ab, vertex.color = "black", vertex.label = NA)
      }else{
        V(g)$ab <- d.inv[i,]*100
        
        plot(g, edge.color = E(g)$sign, edge.curved = F, edge.width = E(g)$weights*3, edge.arrow.size = 1.5, 
             layout = kklay[as.numeric(names(V(g))),], vertex.size = V(g)$ab, vertex.color = "black", vertex.label = NA)
      }
      
    }
  }, 
  movie.name = "netdyn2a.gif", interval = 0.001, ani.width = 600, ani.height = 600, outdir = "~/Desktop/"
)



saveGIF(
  {
    for(i in seq(5,1000,5)){
      if(i == 1){
        par(mar = c(5,5,3,3))
        matplot(d.rem[1:2,], typ = "l", lwd = 3, xlab = "", ylab = "", xlim = c(0,1000), ylim = c(0,.2))
        title(ylab = list("Abundance", cex = 2), xlab = list("Time", cex = 2))
        points(x = 500, y = dyn[[13]][2000,-1][7], pch = 4, cex = 2, lwd = 2, col = "red")
        text(x = 600, y = dyn[[13]][2000,-1][7], label = "<- Extinction", cex = 1.5)
        abline(v = i)
      }else{
        par(mar = c(5,5,3,3))
        matplot(d.rem[1:i,], typ = "l", lwd = 3, xlab = "", ylab = "", xlim = c(0,1000), ylim = c(0,.2))
        title(ylab = list("Abundance", cex = 2), xlab = list("Time", cex = 2))
        points(x = 500, y = dyn[[13]][2000,-1][7], pch = 4, cex = 2, lwd = 2, col = "red")
        text(x = 600, y = dyn[[13]][2000,-1][7], label = "<- Extinction", cex = 1.5)
        abline(v = i, xpd = F)
      }
    }
  }, 
  movie.name = "extdyn1.gif", interval = 0.001, ani.width = 950, ani.height = 800, outdir = "~/Desktop/"
)


saveGIF(
  {
    for(i in seq(5,1000,5)){
      if(i == 1){
        par(mar = c(5,5,3,3))
        matplot(d.inv[1:2,], typ = "l", lwd = 3, xlab = "", ylab = "", xlim = c(0,1000), ylim = c(0,.3))
        title(ylab = list("Abundance", cex = 2), xlab = list("Time", cex = 2))
        points(x = 500, y = .2, pch = 10, cex = 2, lwd = 2)
        text(x = 400, y = .2, label = "Invasion ->", cex = 1.5)
        abline(v = i)
      }else{
        par(mar = c(5,5,3,3))
        matplot(d.inv[1:i,], typ = "l", lwd = 3, xlab = "", ylab = "", xlim = c(0,1000), ylim = c(0,.3))
        title(ylab = list("Abundance", cex = 2), xlab = list("Time", cex = 2))
        points(x = 500, y = .2, pch = 10, cex = 2, lwd = 2)
        text(x = 400, y = .2, label = "Invasion ->", cex = 1.5)
        abline(v = i, xpd = F)
      }
    }
  }, 
  movie.name = "invdyn1.gif", interval = 0.001, ani.width = 950, ani.height = 800, outdir = "~/Desktop/"
)

saveGIF(
  {
    for(i in seq(5,2000,5)){
      if(i == 1){
        par(mar = c(5,5,3,3))
        matplot(dyn[[9]][1:2,-1], typ = "l", lwd = 2, xlab = "", ylab = "", xlim = c(0,2000))
        title(ylab = list("Abundance", cex = 2), xlab = list("Time", cex = 2))
        abline(v = i, xpd = F)
      }else{
        par(mar = c(5,5,3,3))
        matplot(dyn[[9]][1:i,-1], typ = "l", lwd = 2, xlab = "", ylab = "", xlim = c(0,2000))
        title(ylab = list("Abundance", cex = 2), xlab = list("Time", cex = 2))
        abline(v = i, xpd = F)
      }
    }
  }, 
  movie.name = "thousand.gif", interval = 0.001, ani.width = 950, ani.height = 800, outdir = "~/Desktop/"
)



saveHTML(
  {
    for(i in 1:2000){
      V(g)$ab <- (d[i,-1]*10)
      par(mfrow = c(2,1), mar = c(3,4,0,0))
      plot(g, edge.color = E(g)$sign, edge.curved = F, edge.width = E(g)$weights*3, edge.arrow.size = 1.5, 
           layout = kklay, vertex.size = V(g)$ab, vertex.color = "black", vertex.label = NA, frame = F)
      if(i == 1){matplot(d[1:2,-1], typ = "l", xlim = c(0,2000)); abline(v = i, xpd = F)}else{matplot(d[1:i,-1], typ = "l", xlim = c(0,2000));abline(v = i, xpd = F)}
      
    }
  }, 
  movie.name = "dyn.gif", interval = 0.02, nmax = 1000, ani.width = 1200, ani.height = 1000, outdir = "~/Desktop/"
)


?plot.igraph
library(ggraph)
ggraph(g) + geom_edge_arc(aes(edge_colour = factor(sign)),arrow = arrow(length = unit(5, 'mm'))) + geom_node_point(aes(size = ab), show.legend = F)  

t2k <- rad_info(dfin, multityp.fill)
t10 <- rad_info(dfin2, multityp.fill)
t50 <- rad_info(dfin3, multityp.fill)
t100 <- rad_info(dfin4, multityp.fill)
pairs(cbind(t2k$sV, t10$sV, t50$sV, t100$sV))

median(t2k$sV)
median(t100$sV)
median(t50$sV)
median(t10$sV)

dfin[[24]]


mt1 <- multityp.fill[!is.na(dfin)][[2]][dfin[!is.na(dfin)][[2]] > 0,dfin[!is.na(dfin)][[2]] > 0]
ab1 <- dfin[!is.na(dfin)][[2]][dfin[!is.na(dfin)][[2]]>0]

library(rootSolve)
mre2 <- c()
for(i in 1:10000){
  jf1 <- jacobian.full(ab1, lvmodK, parms = list(alpha = runif(length(ab1),0,.1), m = mt1, K = 100))
  mre2[i] <- max(Re(eigen(jf1)$values))
}

#################################################################################################
#################################################################################################
#################################################################################################
x <- 2
y <- 4

m1 <- multityp.fill[[x]][dfin[[x]] > 0,dfin[[x]] > 0]
m2 <- multityp.fill[[y]][dfin[[y]] > 0,dfin[[y]] > 0]

m3 <- matrix(0, nrow = nrow(m1), ncol = ncol(m2))
m3[sample(1:length(m3), 10)] <- runif(10, -.1, 1)
m4 <- matrix(0, nrow = nrow(m2), ncol = ncol(m1))
m4[sample(1:length(m4), 10)] <- runif(10, -.1, 1)

m.all <- (cbind(rbind(m1, m4),rbind(m3,m2)))

init.ab <- c(dfin[[x]][dfin[[x]] > 0], dfin[[y]][dfin[[y]] > 0])
alph <- c(alphas[[x]][dfin[[x]] > 0], alphas[[y]][dfin[[y]] > 0])

par1 <- list(alpha = alph, m = m.all, K = 100)
dyn <-(ode(init.ab, times = 1:2000, func = lvmodK, parms = par1, events = list(func = ext1, time =  1:2000)))

matplot(dyn[,-1], typ = "l", main = i)

fzmod(dyn[2000,-1][dyn[2000,-1] > 0])
