load("~/Downloads/sim.Rdata")

############################################################################################################
############################################################################################################
############################################################################################################
x <- fill_mat(tat)

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
  
  #mean(apply(df1[df1$inty == "amensalism",1:2], 1, function(x) (x[x < 0])))
  #mean(apply(df1[df1$inty == "commensalism",1:2], 1, function(x) (x[x < 0])))
  #mean(apply(df1[df1$inty == "competition",1:2], 1, function(x) (x[x < 0])))
  #mean(apply(df1[df1$inty == "mutualism",1:2], 1, function(x) (x[x < 0])))
  #mean(apply(df1[df1$inty == "predation",1:2], 1, function(x) (x[x < 0])))
  
  num[iL %in% aggregate(df1$inty, list(df1$inty), length)$Group.1] <- aggregate(df1$inty, list(df1$inty), length)$x
  umu[iL %in% aggregate(df1$inty, list(df1$inty), length)$Group.1] <- aggregate(df1$i2, list(df1$inty), mean)$x
  lmu[iL %in% aggregate(df1$inty, list(df1$inty), length)$Group.1] <- aggregate(df1$i1, list(df1$inty), mean)$x
  iL <- c("amensalism", "commensalism", "competition", "mutualism", "predation")
  
  df2 <- data.frame(typ = iL, num = num, m1 = umu, m2 = lmu)
  
  return(df2)
}

############################################################################################################
############################################################################################################
############################################################################################################


plot(hist(sort(ge.mult2$eqst[[1]], decreasing = T), plot = F)$counts)

igraph::fit_power_law(hist(sort(ge.mult2$eqst[[1]], decreasing = T), plot = F)$counts)

table(ge.mult2$eqst[[1]])

fpl3[[1]]

par(mfrow = c(5,5))
for(i in 1:25){
  plot(sort(ge.mult2$eqst[[i]], decreasing = T))
  text(100,.1, label = sapply(fpl2,"[[",5))
}

dev.off()


plot(as.vector(table(ceiling(ge.mult2$eqst[[1]]/min(ge.mult2$eqst[[1]]))))~as.numeric(names(table(ceiling(ge.mult2$eqst[[1]]/min(ge.mult2$eqst[[1]]))))))
hist(sort(ge.mult2$eqst[[1]], decreasing = T))


length(ge.mult2$eqst[[]])
p1 <- ge.mult2$eqst[[1]]/sum(ge.mult2$eqst[[1]])
samp1 <- sample(1:327, 100000, prob = p1, replace = T)
spp <- table(samp1)
table(as.vector(spp))

plot(as.vector(table(as.vector(table(samp1))))~as.numeric(names((table(as.vector(table(samp1)))))))
par(mfcol = c(3,10))
for(x in 1:10){
  abund <- ge.mult2$eqst[[x]]#otu3[,x][otu3[,x] !=0]
  r.ab <- signif(abund/sum(abund),6)
  N = c(10^3, 5*10^3, 10^4, 5*10^4, 10^5)
  t3 <- list()
  for(i in 1:3){
    samp2 <- sample(1:length(abund), N[i], replace = T, prob = r.ab)
    t1 <- table(samp2)
    t2 <- table(as.vector(t1))
    plot(as.vector(t2)~as.numeric(names(t2)), main = N[i], xlab = "Abundance", ylab = "Freq")
  }
}

r.ab <- (otu3[,1][otu3[,1] !=0]/sum(otu3[,1][otu3[,1] !=0]))



x <- 3
abund <- otu3[,x][otu3[,x] !=0]
r.ab <- abund/sum(abund)
N = c(10^3, 5*10^3, 10^4, 5*10^4, 10^5)
t3 <- list()
par(mfrow = c(5,2))
for(i in 1:10){
  #s#amp2 <- sample(1:length(r.ab), N[i], replace = T, prob = r.ab)
  #t1 <- table(samp2)
  #t2 <- table(as.vector(t1))
  t2 <- table(otu3[,i][otu3[,i] !=0])
  plot(as.vector(t2)~as.numeric(names(t2)), main = N[i], xlab = "Abundance", ylab = "Freq", xlim = c(0,300))
  
  t3[[i]] <- as.numeric(t1)
}

dev.off()
plot(r.ab[1:length(r.ab) %in% as.numeric(names(t1))], as.vector(t1)/sum(as.vector(t1)))
  
lapply(t3, fitsad, "pareto")

## 

get_abundvec <- function(abund, N = 10000){
  r.ab <- abund/sum(abund)
  samp2 <- sample(1:length(abund), N, replace = T, prob = r.ab)
  t1 <- table(samp2)
  return(as.numeric(t1))
}

avec <- lapply(ge.hub3$eqst, get_abundvec)
get_bestfit <- function(avec){
  fits1 <- lapply(avec, fitsad, sad = "lnorm")
  fits2 <- lapply(avec, fitsad, sad = "power")
  fits3 <- lapply(avec, fitsad, sad = "powbend")
  fits4 <- lapply(avec, fitsad, sad = "mzsm")
  fits5 <- lapply(avec, fitsad, sad = "poilog")
  fits6 <- lapply(avec, fitsad, sad = "bs")
  fits7 <- lapply(avec, fitsad, sad = "ls")
  fits8 <- lapply(avec, fitsad, sad = "weibull")
  fit1 <- sapply(fits1, AIC)
  fit2 <- sapply(fits2, AIC)
  fit3 <- sapply(fits3, AIC)
  fit4 <- sapply(fits4, AIC)
  fit5 <- sapply(fits5, AIC)
  fit6 <- sapply(fits6, AIC)
  fit7 <- sapply(fits7, AIC)
  fit8 <- sapply(fits8, AIC)
  #sapply(fits5, function(x) c(coef(x), AICvol = AIC(x)))
  t.wins <- table(apply(cbind(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8), 1, 
                        function(x){c("lnorm", "power", "powbend", "mzsm", "poilog", "bs", "ls", "weibull")[which.min(x)]}))
  #tAIC <- cbind(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8)
  #colnames(tAIC) <- c("lnorm", "power", "powbend", "mzsm", "poilog", "bs", "ls", "weibull")
  #return(tAIC)
  return(t.wins)
}

gbf1m <- get_bestfit3(lapply(ge.mult$eqst, get_abundvec))
gbf2m <- get_bestfit3(lapply(ge.mult2$eqst, get_abundvec))
gbf3m <- get_bestfit3(lapply(sapply(ge.mult3$eqst, function(x) x[x > 0] ), get_abundvec))

gbf1t <- get_bestfit3(lapply(ge.tat$eqst, get_abundvec))
gbf2t <- get_bestfit3(lapply(ge.tat2$eqst, get_abundvec))
gbf3t <- get_bestfit3(lapply(ge.tat3$eqst, get_abundvec))


gbfh1 <- get_bestfit3(lapply(ge.hub$eqst, get_abundvec))
gbfh2 <- get_bestfit3(lapply(ge.hub2$eqst, get_abundvec))
gbfh3 <- get_bestfit3(lapply(ge.hub3$eqst, get_abundvec))

boxplot(t(apply(gbf3, 1, function(x) x-min(x))))
allgbfA <- rbind(cbind(melt(t(apply(gbf1, 1, function(x) x-min(x)))), typ = 1),
      cbind(melt(t(apply(gbf2, 1, function(x) x-min(x)))), typ = 2),
      cbind(melt(t(apply(gbf3, 1, function(x) x-min(x)))), typ = 3))
ggplot(allgbfA, aes(x = Var2, y = value)) + geom_boxplot(aes(fill = factor(typ))) + facet_grid(~factor(typ))

sad.mods <- c("bs","gamma","geom","lnorm","ls","mzsm","nbinom","pareto", "poilog","power", "powbend", "volkov","weibull")

##################################################################
##################################################################
##################################################################

otu2 <- read.csv("~/Desktop/otu_table_psn_v13.csv", row.names = 1)
metadat <- read.csv("~/Desktop/v13_map_uniquebyPSN.csv")

stoolsamp <- which(metadat$HMPbodysubsite == "Stool")
spptab <- colnames(otu2) %in% paste0("X",metadat[stoolsamp,]$SampleID)
otu3 <- otu2[-which(rowSums(otu2[,spptab]) == 0),spptab]

rad.fitted <- lapply(1:ncol(otu3), function(x) vegan::radfit(otu3[,x][otu3[,x] != 0]))
rad.fitted2 <- lapply(1:ncol(otu3), function(x) vegan::radfit(sort(otu3[,x][otu3[,x] != 0], decreasing = T)[1:50]))
plot(rad.fitted2[[1]])
table(sapply(rad.fitted, function(x) which.min(sapply(x$models, AIC))))
table(sapply(rad.fitted2, function(x) which.min(sapply(x$models, AIC))))

apply(otu3, 2, function(x) sum(x != 0))
nspp <- seq(10, 170, 10)
res1 <- list()
for(i in 1:length(nspp)){
  rad.fitts <- lapply(1:ncol(otu3)[-c(184,166,165)], function(x) vegan::radfit(sort(otu3[,x][otu3[,x] != 0], decreasing = T)[1:nspp[i]]))
  res1[[i]] <- table(sapply(rad.fitts, function(x) which.min(sapply(x$models, AIC))))
}

x <- 2
evalWithTimeout(plot(radpred(fitpoilog((sort(otu3[,x][otu3[,x] != 0], decreasing = T)))), pch = 20), timeout = 120)
points(sort(otu3[,x][otu3[,x] != 0], decreasing = T), col = "blue", pch = 20)



get_bestfit(lapply(1:ncol(otu3), function(x) otu3[,x][otu3[,x] != 0]))

otuAIC <- get_bestfit(lapply(1:ncol(otu3), function(x) otu3[,x][otu3[,x] != 0]))
boxplot(t(apply(otuAIC, 1, function(x) x-min(x))))
otuAIC[,1:2]

##################################################################
##################################################################
##################################################################

ga1 <- get_abundvec(eq1)
fsp1  <- fitsad(ga1, "power")
fsp2  <- fitsad(ga1, "ls")
pred <- dpower(sort(unique(ga1)), fsp1@coef)
pred2 <- dls(sort(unique(ga1)), 10000, alpha = fsp2@coef)
obs <- as.vector(table(ga1))/sum(as.vector(table(ga1)))
sum((pred-obs)^2)
sum((pred2-obs)^2)
AIC(fsp1)
AIC(fsp2)


nsp <- sapply(ge.mult$eqst, length)

ga1 <- lapply(ge.mult$eqst[-4], get_abundvec, N = 1000)
obs1 <- lapply(ga1, function(x)  cbind(as.vector(table(x))/sum(as.vector(table(x))), as.numeric(names(table(x)))))
fsp1 <- lapply(ga1, fitpoilog)
pred1 <- lapply(1:length(ga1), function(x) dpoilog(sort(unique(ga1[[x]])), fsp1[[x]]@coef[1], fsp1[[x]]@coef[2]))

par(mfrow = c(1,5))
for(i in 1:length(ga1)){
  plot(obs1[[i]][,2:1], typ = "o", main = nsp[-4][i])
  points(pred1[[i]]~obs1[[i]][,2], pch = 20, col = "blue", typ = "o")
}


ga2 <- lapply(ge.mult$eqst, get_abundvec, N = 2000)
obs2 <- lapply(ga2, function(x) cbind(as.vector(table(x))/sum(as.vector(table(x))), as.numeric(names(table(x)))))
fsp2 <- lapply(ga2, fitpoilog)
pred2 <- lapply(1:length(ga2), function(x) dpoilog(sort(unique(ga2[[x]])), fsp2[[x]]@coef[1], fsp2[[x]]@coef[2]))

par(mfrow = c(1,5))
for(i in 1:length(ga2)){
  plot(obs2[[i]][,2:1], typ = "o", main = nsp[i])
  points(pred2[[i]]~obs2[[i]][,2], pch = 20, col = "blue", typ = "o")
}

ga3 <- lapply(ge.mult$eqst, get_abundvec, N = 10000)
obs3 <- lapply(ga3, function(x) cbind(as.vector(table(x))/sum(as.vector(table(x))), as.numeric(names(table(x)))))
fsp3 <- lapply(ga3, fitpoilog)
pred3 <- lapply(1:length(ga3), function(x) dpoilog(sort(unique(ga3[[x]])), fsp3[[x]]@coef[1], fsp3[[x]]@coef[2]))

par(mfrow = c(1,5))
for(i in 1:length(ga3)){
  plot(obs3[[i]][,2:1], typ = "o", main = nsp[i])
  points(pred3[[i]]~obs3[[i]][,2], pch = 20, col = "blue", typ = "o")
}

get_r2 <- function(o, p){
  1 - sum((o-p)^2)/sum((o - mean(o))^2)
}

sapply(1:length(ga2), function(x){get_r2(obs2[[x]][,1], pred2[[x]])})
sapply(1:length(ga2), function(x){get_r2(obs3[[x]][,1], pred3[[x]])})

fsp1 <- lapply(ga1, fitpower)
pred1 <- lapply(1:length(ga1), function(x) dpower(sort(unique(ga1[[x]])), fsp1[[x]]@coef))
fsp2 <- lapply(ga1, fitls)
pred2 <- lapply(1:length(ga1), function(x) dls(sort(unique(ga1[[x]])), 10000, fsp1[[x]]@coef))
obs <- lapply(ga1, function(x)  as.vector(table(x))/sum(as.vector(table(x))))

sse1 <- sapply(1:length(obs), function(x) sum((pred1[[x]] - obs[[x]])^2))
sse2 <- sapply(1:length(obs), function(x) sum((pred2[[x]] - obs[[x]])^2))

plot(sse1, sapply(fsp1, AIC))
plot(sse2, sapply(fsp2, AIC))
##################################################################
##################################################################
##################################################################
library(vegan)
get_abundvec <- function(abund, N = 10000){
  r.ab <- abund/sum(abund)
  samp2 <- sample(1:length(abund), N, replace = T, prob = r.ab)
  t1 <- table(samp2)
  return(as.numeric(t1))
}


#lf1 <- list.files("D:/jjborrelli/parSADhub_data/")
lf1 <- list.files("~/Documents/Data/parSAD_data/")
lf2 <- grep("ge", lf1)
#for(i in 1:length(lf1[lf2])){
rads <- list()
rads2 <- list()
nspp <- c()
f <- c()
for(i in 1:length(lf2)){
  ge1 <- readRDS(paste("~/Documents/Data/parSAD_data/", lf1[lf2][[i]], sep = ""))
  if(any(is.na(ge1))){next}
  if(any(is.na(ge1$eqst))){next}
  eqa1 <- ge1$eqst/sum(ge1$eqst)
  eqa1[eqa1 < 0] <- 0
  gav1 <- get_abundvec(eqa1, N = length(ge1$eqst)*2)
  rads[[i]] <- radfit(gav1)
  gav2 <- get_abundvec(eqa1, N = 2000)
  rads2[[i]] <- radfit(gav2)
  #plot(as.vector(table(get_abundvec(ge1$eqst, 1000))))
  nspp[i] <- length(ge1$eqst)
  f[i] <- lf1[lf2][[i]]
  
  gav3 <- ceiling(eqa1*5000)
   
  rads3[[i]] <- radfit(gav3)
  gav4 <- ceiling(eqa1*1000)
  rads4[[i]] <- radfit(gav4)
  
  print(i)
}



raic <- sapply(rads[!sapply(rads, is.null)], AIC)
table(apply(raic, 2, which.min))

raic2 <- sapply(rads2[!sapply(rads2, is.null)], AIC)
table(apply(raic2, 2, which.min))

raic3 <- sapply(rads3[!sapply(rads3, is.null)], AIC)
table(apply(raic3, 2, which.min))

raic4 <- sapply(rads4[!sapply(rads4, is.null)], AIC)
table(apply(raic4, 2, which.min))

nspp1 <- sapply(1:nrow(eqabs), function(x){if(any(is.na(eqabs[x,]))){return(NA);next};length(eqabs[x,eqabs[x,] != 0]/sum(eqabs[x,eqabs[x,] != 0])*5000)})
npred <- sapply(1:length(rads3), function(x){if(is.na(nspp1[x])){return(NA)};nrow(predict(rads3[[x]]))})


SSE <- lapply(1:length(rads3), function(x) if(is.null(rads3[[x]])){return(NA)}else{apply(predict(rads3[[x]]), 2, function(z) sum((z - eqabs[x,eqabs[x,]!=0])^2))})

Rsq <- lapply(1:length(rads3), function(z){
  if(is.na(nspp1[z])){return(NA)}
  obsAB <- ceiling(eqabs[z,eqabs[z,]!=0]/sum(eqabs[z,eqabs[z,]!=0])*5000)
  SST <- apply(predict(rads3[[z]]), 2, function(x) sum((x - mean(obsAB))^2))
  SSE <- apply(predict(rads3[[z]]), 2, function(x) sum((x - obsAB)^2))
  
  return((SST - SSE)/SST)
})

rsqmat <- do.call(rbind, Rsq[!is.na(Rsq)])

rsqdat <- (do.call(rbind, Rsq[!is.na(Rsq)]))
rsqdat <- melt(rsqdat)
rsqdat <- cbind(rsqdat, N = rep(nspp1[!is.na(nspp1)], 5))
domint <- sapply(eqmat, function(x){if(is.null(x)){return(NA)};x$typ[which.max(x$num)]})
alldom <- rep(domint[!is.na(domint)],5)[rsqdat$N > 100]

ggplot(rsqdat[rsqdat$N > 100,], aes(x = N, y = value)) + geom_point(col = alldom) + geom_smooth(method = "loess") + facet_grid(~Var2)
ggplot(rsqdat[rsqdat$N > 100,], aes(x = N, y = value)) + geom_smooth(aes(N, value, colour = factor(alldom))) + facet_grid(~Var2)
rsa1 <- c()
for(z in 1:length(rads3)){
  if(is.na(nspp1[z])){(NA);next}
  obsAB <- ceiling(eqabs[z,eqabs[z,]!=0]/sum(eqabs[z,eqabs[z,]!=0])*5000)
  SST <- apply(predict(rads3[[z]]), 2, function(x) sum((x - mean(obsAB))^2))
  SSE <- apply(predict(rads3[[z]]), 2, function(x) sum((x - obsAB)^2))
  warnings()
  ((SST - SSE)/SST)
  print(z)
}

##################################################################
##################################################################

lf1 <- list.files("~/Documents/Data/parSAD_data/")
lf2 <- grep("ge", lf1)
lf3 <- grep("mat", lf1)

eqmat <- list()
iconn <- c()
for(i in 1:length(lf2)){
  ge1 <- readRDS(paste("~/Documents/Data/parSAD_data/", lf1[lf2][[i]], sep = ""))
  if(any(is.na(ge1))){next}
  if(any(is.na(ge1$eqst))){next}
  mat1 <- readRDS(paste("~/Documents/Data/parSAD_data/", lf1[lf3][[i]], sep = ""))
  
  iconn[i] <- is.connected(graph.adjacency(abs(sign(mat1[ge1$spp, ge1$spp]))))
  
  eqmat[[i]] <- data.frame(itystr(mat1[ge1$spp, ge1$spp]), web = i, N = sum(ge1$spp))
  print(i)
}

int1 <- t(sapply(eqmat, function(x){if(is.null(x)){return(c(NA,NA,NA,NA,NA))};x$num}))
lapply(head(eqmat,10), function(x){is.null(x)})
maxRsq <- apply(rsqmat, 1, which.max)

int2 <- int1[complete.cases(int1),]
cart1 <- (rpart::rpart(maxRsq~int2+nspp1[!is.na(nspp1)]))
plot(cart1)
text(cart1)

istrs <- lapply(eqmat, function(x){x$m1[2] <- x$m2[2]; x$m2[2] <- 0;return(x)})
int3 <- t(sapply(istrs, function(x){if(any(is.na(unlist(x)))){return(c(NA,NA,NA,NA,NA))};x$m1}))
int4 <- int3[complete.cases(int3),]

lnpar <- t(sapply(rads3, function(x){if(is.null(x)){return(c(NA, NA))};x$models$Lognormal$coefficients}))
pepar <- (sapply(rads3, function(x){if(is.null(x)){return(c(NA))};x$models$Preemption$coefficients}))
zmpar <- t(sapply(rads3, function(x){if(is.null(x)){return(c(NA, NA, NA))};x$models$Mandelbrot$coefficients}))
isinf <- apply(zmpar, 1, function(x) is.infinite(x[1]))
zpar <- t(sapply(rads3, function(x){if(is.null(x)){return(c(NA, NA))};x$models$Zipf$coefficients}))


summary(lm(lnpar[complete.cases(lnpar),]~int2+int4+nspp1[!is.na(nspp)]))
summary(lm(pepar[!is.na(pepar)]~int2+nspp1[!is.na(nspp)]))
summary(lm(zmpar[!isinf,][complete.cases(zmpar[!isinf,]),]~int1[!isinf,][complete.cases(zmpar[!isinf,]),]+nspp1[!isinf][complete.cases(zmpar[!isinf,])]))
summary(lm(zpar[complete.cases(zpar),]~int2+int4+nspp1[!is.na(nspp)]))
head(eqmat, 10)

eqabs <- matrix(0, nrow = length(lf2), ncol = 1000)
for(i in 1:length(lf2)){
  ge1 <- readRDS(paste("~/Documents/Data/parSAD_data/", lf1[lf2][[i]], sep = ""))
  if(any(is.na(ge1))){eqabs[i,] <- NA; next}
  if(any(is.na(ge1$eqst))){eqabs[i,] <- NA;next}
  eqabs[i,1:length(ge1$eqst)] <- sort(ge1$eqst, decreasing = T)
  print(i)
}


mandCoef <- sapply(rads2[which(apply(raic2, 2, which.min) == 5)][!sapply(rads2[which(apply(raic2, 2, which.min) == 5)], is.null)], function(x) x$models$Mandelbrot$coefficients)
mandCoef[,1:200]
plot(rads2[[1]])




i = 1
plot(log10(eqabs[i, eqabs[i,] != 0])~log10(1:length(eqabs[i, eqabs[i,] != 0])))
abline(lm(log10(eqabs[i, eqabs[i,] != 0])~log10(1:length(eqabs[i, eqabs[i,] != 0]))))

slp <- c()
rsq <- c()
for(i in 1:nrow(eqabs)){
  if(any(is.na(eqabs[i,eqabs[i,]!=0]))){slp[i] <- NA; rsq[i] <- NA; next}
  fit <- (lm(log10(eqabs[i, eqabs[i,] != 0])~log10(1:length(eqabs[i, eqabs[i,] != 0]))))
  slp[i] <- fit$coefficients[2]
  rsq[i] <- summary(fit)$r.squared
}

plot(slp~rsq)

##################################################################
##################################################################
##################################################################
x1 <- ge.mult
x2 <- ge.tat
x3 <- ge.hub

mt1 <- lapply(1:length(x1$wrk), function(x){multityp[[x1$wrk[[x]]]][x1$spp[[x]],x1$spp[[x]]]})
mt2 <- lapply(1:length(x2$wrk), function(x){multitat[[x2$wrk[[x]]]][x2$spp[[x]],x2$spp[[x]]]})
mt3 <- lapply(1:length(x3$wrk), function(x){multihub[[x3$wrk[[x]]]][x3$spp[[x]],x3$spp[[x]]]})


itym <- t(sapply(mt1, itypes)) # t(sapply(multityp, itypes)[,ge.mult$wrk])
ityt <- t(sapply(mt2, itypes)) # t(sapply(multitat, itypes)[,ge.tat$wrk])
ityh <- t(sapply(mt3, itypes)) # t(sapply(multihub, itypes)[,ge.hub$wrk])


grpsm <- (apply(gbf1m, 1, function(x) names(which.min(x))))
grpst <- (apply(gbf1t, 1, function(x) names(which.min(x))))
grpsh <- (apply(gbfh1, 1, function(x) names(which.min(x))))

datm <- data.frame(grpsm, itym)
datt <- data.frame(grpst, ityt)
dath <- data.frame(grpsh, ityh)

s1m <- sample(1:nrow(datm), nrow(datm)/2)
s1t <- sample(1:nrow(datt), nrow(datt)/2)
s1h <- sample(1:nrow(dath), nrow(dath)/2)

ldafit1 <- MASS::lda(grpsm~comp+mut+pred+amens+comm, data = datm[s1m,])
plda1 <- predict(ldafit1, newdata = datm[-s1m,])
tab1 <- table(pred = as.character(plda1$class), obs = grpsm[-s1m])
sum(diag(tab1))/sum(tab1)
plot(predict(ldafit1)$x, col = predict(ldafit1)$class, pch = 20)

ldafit2 <- MASS::lda(grpst~comp+mut+pred+amens+comm, data = datt[s1t,])
plda2 <- predict(ldafit2, newdata = datt[-s1t,])
tab2 <- table(pred = as.character(plda2$class), obs = grpst[-s1t])
sum(diag(tab2))/sum(tab2)
plot(predict(ldafit2)$x, col = predict(ldafit2)$class, pch = 20)

ldafit3 <- MASS::lda(grpsh~comp+mut+pred+amens+comm, data = dath[s1h,])
plda3 <- predict(ldafit3, newdata = dath[-s1h,])
tab3 <- table(pred = as.character(plda3$class), obs = grpsh[-s1h])
sum(diag(tab3))/sum(tab3)
plot(predict(ldafit3)$x, col = predict(ldafit3)$class, pch = 20)
