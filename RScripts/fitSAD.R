load("~/Downloads/sim.Rdata")


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

get_bestfit(lapply(ge.mult$eqst, get_abundvec))
get_bestfit(lapply(ge.mult2$eqst, get_abundvec))
get_bestfit(lapply(sapply(ge.mult3$eqst, function(x) x[x > 0] ), get_abundvec))

get_bestfit(lapply(ge.tat$eqst, get_abundvec))
get_bestfit(lapply(ge.tat2$eqst, get_abundvec))
get_bestfit(lapply(ge.tat3$eqst, get_abundvec))


gbf1 <- get_bestfit(lapply(ge.hub$eqst, get_abundvec))
gbf2 <- get_bestfit(lapply(ge.hub2$eqst, get_abundvec))
gbf3 <- get_bestfit(lapply(ge.hub3$eqst, get_abundvec))

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

get_bestfit(lapply(1:ncol(otu3), function(x) otu3[,x][otu3[,x] != 0]))

otuAIC <- get_bestfit(lapply(1:ncol(otu3), function(x) otu3[,x][otu3[,x] != 0]))
boxplot(t(apply(otuAIC, 1, function(x) x-min(x))))
