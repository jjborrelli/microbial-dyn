#in case of error - saved to Docs/Abund/immigration.Rdata
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

lvmodKI <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha * (1 - state/(parms$K/sum(state > 10^-10))) + state * parms$m %*% state + 10^-5
    
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


iter = 150
reslst <- list()
parlst <- list()
inab <- cbind(inab, matrix(nrow = 600, ncol = 100))
tfz <- list()
for(i in 51:iter){
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
  
  #matplot(dyn[,-1], typ = "l", main = i)
  
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

plot(log10(allf$s)~log10(allf$N), ylim = c(-1,1), xlim = c(.5,3))
#ga <- lapply(lapply(reslst, "[[", 2), get_abundvec)
#do.call(rbind,lapply(ga, fzmod))
rb1 <- do.call(rbind, lapply(lapply(lapply(reslst2, "[[",1), get_abundvec, N = 2000), fzmod))
points(rb1[,1:2], col = "green4", pch = 20)

points(do.call(rbind ,sapply(tfz, function(x) x[500:2000,1:2])), pch = 20, col = "green4")
points(log10(do.call(rbind, lapply(tfz, function(x) x[2000,1:2]))), pch = 20, col = "blue")

reslst3 <- list()
tfz2 <- list()
for(i in 63:length(parlst)){
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
points(log10(do.call(rbind, lapply(tfz2, function(x) x[2000,1:2]))), pch = 20, col = "gray")


Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor, lwd = 2)
}  

points(log10(alld$N), log10(alld$fz), col = "red1")
plot()
Plot_ConvexHull(allf$N, allf$s, "black")
pt1 <- (do.call(rbind ,sapply(tfz, function(x) x[500:2000,1:2])))
Plot_ConvexHull(pt1[,1], pt1[,2], "blue")
pt2 <- log10(do.call(rbind, lapply(tfz2, function(x) x[500:2000,1:2])))
Plot_ConvexHull(pt2[,1], pt2[,2], "darkslategray")
Plot_ConvexHull(log10(alld$N), log10(alld$fz), "red")

fit1 <- rq(log10(s) ~ log10(N), tau = .01, data = allf)
fit2 <- rq(log10(s) ~ log10(N), tau = .001, data = allf)
abline(fit2, lty = 2)





test.im <- function(n, inab, parlst, dFunc = "lvmod", tmax = 2000){
  if(n > length(parlst)){stop("Too many iterations")}
  fitfz <- list()
  for(i in 1:n){
    if(dFunc == "lvmod"){
      dyn <- ode(inab[,i], times = 1:tmax, func = lvmod, parms = parlst[[i]], events = list(func = ext1, time =  1:tmax))
    }else if(dFunc == "lvmodK"){
      dyn <- ode(inab[,i], times = 1:tmax, func = lvmodK, parms = parlst[[i]], events = list(func = ext1, time =  1:tmax))
    }else if(dFunc == "lvmodKI"){
      dyn <- ode(inab[,i], times = 1:tmax, func = lvmodKI, parms = parlst[[i]], events = list(func = ext1, time =  1:tmax))
    }else{
      stop("Set dFunc to either lvmod or lvmodK or lvmodKI")
    }
    if(nrow(dyn) != 2000){
      #reslst3[[i]] <- list(NA, NA)
      tfz2[[i]] <- data.frame(N = NA, s = NA, nll = NA, r2 = NA)
      next
    }else{
      spp <- dyn[tmax, -1] > 0
      dfin <- dyn[tmax, -1]
      #conn <- sum(parlst$m[spp, spp] != 0)/(sum(spp)*sum(spp))
      #reslst3[[i]] <- list(dfin, conn)
    }
    
    tga <- apply(dyn[,-1], 1, function(x) get_abundvec(x[x > 0], 2000))
    fitfz[[i]] <- do.call(rbind, lapply(tga, fzmod))
    print(i)
  }
  
  return(fitfz)
}




# Set immigration to 10^-3
lvmod <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha + state * parms$m %*% state + 10^-3
    list(dB)
  })
}
highim <- test.im(n = 150, inab, parlst, "lvmod")
points(log10(do.call(rbind,lapply(highim, function(x) x[2000,1:2]))), pch = 0, col = "green2")
saveRDS(do.call(rbind,lapply(highim, function(x) x[2000,1:2])), "~/Documents/AbundData/IM_10_3.rds")


# Set immigration to 10^-4
lvmod <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha + state * parms$m %*% state + 10^-4
    list(dB)
  })
}
highim2 <- test.im(n = 150, inab, parlst, "lvmod")
points(log10(do.call(rbind,lapply(highim2, function(x) x[2000,1:2]))), pch = 20, col = "coral2")
saveRDS(do.call(rbind,lapply(highim2, function(x) x[2000,1:2])), "~/Documents/AbundData/IM_10_4.rds")


# Set immigration to 10^-5
#lvmod <- function(times, state, parms){
#  with(as.list(c(state, parms)), {
#    dB <- state * parms$alpha + state * parms$m %*% state + 10^-5
#    list(dB)
#  })
#}
#medim <- test.im(n = 20, inab, parlst, "lvmod")
#points(log10(do.call(rbind,lapply(medim, function(x) x[2000,1:2]))), pch = 20, col = "gray3")
#saveRDS(do.call(rbind,lapply(medim, function(x) x[2000,1:2])), "~/Documents/AbundData/IM_10_5.rds")


# Set immigration to 10^-6
lvmod <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha + state * parms$m %*% state + 10^-6
    list(dB)
  })
}
lowim <- test.im(n = 150, inab, parlst, "lvmod")
points(log10(do.call(rbind,lapply(lowim, function(x) x[2000,1:2]))), pch = 20, col = "gray3")
saveRDS(do.call(rbind,lapply(lowim, function(x) x[2000,1:2])), "~/Documents/AbundData/IM_10_6.rds")


# Set immigration to 10^-7
lvmod <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha + state * parms$m %*% state + 10^-7
    list(dB)
  })
}
lowim2 <- test.im(n = 150, inab, parlst, "lvmod")
points(log10(do.call(rbind,lapply(lowim2, function(x) x[2000,1:2]))), pch = 20, col = "purple")
saveRDS(do.call(rbind,lapply(lowim2, function(x) x[2000,1:2])), "~/Documents/AbundData/IM_10_7.rds")


# Set immigration to 10^-8
lvmod <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha + state * parms$m %*% state + 10^-8
    list(dB)
  })
}
lowim3 <- test.im(n = 150, inab, parlst, "lvmod")
points(log10(do.call(rbind,lapply(lowim3, function(x) x[2000,1:2]))), pch = 20, col = "green4")
saveRDS(do.call(rbind,lapply(lowim3, function(x) x[2000,1:2])), "~/Documents/AbundData/IM_10_8.rds")


tfz2 <- test.im(n = 150, inab, parlst, "lvmodK")

K.im5 <- test.im(n = 150, inab, parlst, "lvmodKI")
points(log10(do.call(rbind,lapply(K.im5, function(x) x[2000,1:2]))), pch = 20, col = "gray1")
saveRDS(do.call(rbind,lapply(K.im5, function(x) x[2000,1:2])), "~/Documents/AbundData/K_IM_10_5.rds")
ptk5 <- log10(do.call(rbind,lapply(K.im5, function(x) x[2000,1:2])))
Plot_ConvexHull(ptk5[,1], ptk5[,2], "gray1")


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
plot(log10(allf$s)~log10(allf$N), ylim = c(-1,1), xlim = c(1.5,3))


# HIGHIM
#points(log10(do.call(rbind, lapply(highim, function(x) x[2000,1:2]))), pch = 20, col = cbPalette[1])
pt1 <- log10(do.call(rbind, lapply(highim, function(x) x[2000,1:2])))
pt1 <- pt1[complete.cases(pt1),]
Plot_ConvexHull(pt1[,1], pt1[,2], cbPalette[1])
text(apply(pt1, 2, median)[1],apply(pt1, 2, median)[2], "1", col = cbPalette[1])

#HIGHIM2
#points(log10(do.call(rbind,lapply(highim2, function(x) x[2000,1:2]))), pch = 0, col = cbPalette[2])
pt2 <- log10(do.call(rbind,lapply(highim2, function(x) x[2000,1:2])))
pt2 <- pt2[complete.cases(pt2),]
Plot_ConvexHull(pt2[,1], pt2[,2], cbPalette[2])
text(apply(pt2, 2, median)[1],apply(pt2, 2, median)[2], "2", col = cbPalette[2])

#TFZ/MEDIM
#points(log10(do.call(rbind,lapply(tfz, function(x) x[2000,1:2]))), pch = 20, col = cbPalette[3])
pt3 <- log10(do.call(rbind,lapply(tfz, function(x) x[2000,1:2])))
pt3 <- pt3[complete.cases(pt3),]
Plot_ConvexHull(pt3[,1], pt3[,2], cbPalette[3])
text(apply(pt3, 2, median)[1],apply(pt3, 2, median)[2], "3", col = cbPalette[3])

#LOWIM
#points(log10(do.call(rbind,lapply(lowim, function(x) x[2000,1:2]))), pch = 20, col = cbPalette[4])
pt4 <- log10(do.call(rbind,lapply(lowim, function(x) x[2000,1:2])))
pt4 <- pt4[complete.cases(pt4),]
Plot_ConvexHull(pt4[,1], pt4[,2], cbPalette[4])
text(apply(pt4, 2, median)[1],apply(pt4, 2, median)[2], "4", col = cbPalette[4])

#LOWIM2
#points(log10(do.call(rbind, lapply(lowim2, function(x) x[2000,1:2]))), pch = 20, col = cbPalette[5])
pt5 <- log10(do.call(rbind, lapply(lowim2, function(x) x[2000,1:2])))
pt5 <- pt5[complete.cases(pt5),]
Plot_ConvexHull(pt5[,1], pt5[,2], cbPalette[5])
text(apply(pt5, 2, median)[1],apply(pt5, 2, median)[2], "5", col = cbPalette[5])

#LOWIM3
#points(log10(do.call(rbind,lapply(lowim3, function(x) x[2000,1:2]))), pch = 20, col = cbPalette[6])
pt6 <- log10(do.call(rbind,lapply(lowim3, function(x) x[2000,1:2])))
pt6 <- pt6[complete.cases(pt6),]
Plot_ConvexHull(pt6[,1], pt6[,2], cbPalette[6])
text(apply(pt6, 2, median)[1],apply(pt6, 2, median)[2], "6", col = cbPalette[6])

#NOIM K
#points(log10(do.call(rbind,lapply(tfz2, function(x) x[2000,1:2]))), pch = 20, col = cbPalette[6])
pt7 <- log10(do.call(rbind,lapply(tfz2, function(x) x[2000,1:2])))
pt7 <- pt6[complete.cases(pt6),]
Plot_ConvexHull(pt7[,1], pt7[,2], cbPalette[7])
text(apply(pt7, 2, median)[1],apply(pt7, 2, median)[2], "7", col = cbPalette[7])

#IM K
#points(log10(do.call(rbind,lapply(K.im5, function(x) x[2000,1:2]))), pch = 20, col = cbPalette[7])
ptk5 <- log10(do.call(rbind,lapply(K.im5, function(x) x[2000,1:2])))
Plot_ConvexHull(ptk5[,1], ptk5[,2], cbPalette[8])
text(apply(ptk5, 2, median)[1],apply(ptk5, 2, median)[2], "8", col = cbPalette[8])