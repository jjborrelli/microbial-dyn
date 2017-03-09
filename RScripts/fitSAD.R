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