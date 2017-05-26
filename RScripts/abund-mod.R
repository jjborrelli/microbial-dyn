library(data.table)
library(rpart)
library(rpart.plot)
library(randomForest)

dat <- readRDS("~/Documents/AbundData/res.rds")
dat2 <- readRDS("~/Documents/AbundData/res2.rds")

dat.fin <- rbindlist(lapply(dat, "[[", 2))
dat.init <- rbindlist(lapply(dat, "[[", 1))
dat.init$com <- rep(1:1000, each = 600)
dat.init$rab <- unlist(lapply(lapply(dat, "[[", 1), function(x) x$ab/sum(x$ab)))
dat.init$rabi <- unlist(lapply(lapply(dat, "[[", 1), function(x) x$ab.i/sum(x$ab.i)))
dat.fin$com <- rep(1:1000, sapply(lapply(dat, "[[", 2), nrow))
dat.fin$rab <- unlist(lapply(lapply(dat, "[[", 2), function(x) x$ab/sum(x$ab)))
dat.fin$rabi <- unlist(lapply(lapply(dat, "[[", 2), function(x) x$ab.i/sum(x$ab.i)))
dat.fin$N <- rep(sapply(lapply(dat, "[[", 2), nrow),sapply(lapply(dat, "[[", 2), nrow))

dat2.fin <- rbindlist(lapply(dat2, "[[", 2))
dat2.init <- rbindlist(lapply(dat2, "[[", 1))
dat2.init$com <- rep(1:1000, each = 600)
dat2.init$rab <- unlist(lapply(lapply(dat2, "[[", 1), function(x) x$ab/sum(x$ab)))
dat2.init$rabi <- unlist(lapply(lapply(dat2, "[[", 1), function(x) x$ab.i/sum(x$ab.i)))
dat2.fin$com <- rep(1:1000, sapply(lapply(dat2, "[[", 2), nrow))
dat2.fin$rab <- unlist(lapply(lapply(dat2, "[[", 2), function(x) x$ab/sum(x$ab)))
dat2.fin$rabi <- unlist(lapply(lapply(dat2, "[[", 2), function(x) x$ab.i/sum(x$ab.i)))


datrel <- rbindlist(lapply(lapply(dat, "[[", 2), function(x){for(i in 1:ncol(x)){x[,i] <- x[,i]/sum(x[,i]); x[,i][is.nan(x[,i])] <- 0};return(x)}))

plot(dat.init$spp~dat.init$ab.i)

fit1 <- lm((ab)~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, data = datrel, x = F, y = F, model = F, na.action = "na.fail")
summary(fit1)
#dfit1 <- MuMIn::dredge(fit1)

fit2 <- lm(rab~K2+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut+allIn+allOut, data = dat.init, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2)
fpart2 <- rpart(log(ab)~K2+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, data = dat.fin, method = "anova")
prp(fpart2, extra = 1)

s1 <- sample(1:nrow(dat.init), 200)
fit3 <- glm(spp~K2+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, data = dat.init[-s1,], family = "binomial",x = F, y = F, model = F, na.action = "na.fail")
summary(fit3)
boxplot(predict.glm(fit3, dat.init[s1,], type = "response")~dat.init$spp[s1])


maxab <- sapply(lapply(dat,"[[", 2), function(x) max(x$ab))
summary(maxab)
maxab2 <- sapply(lapply(dat,"[[", 2), function(x) max(x$ab[-which.max(x$ab)]))

summary(lm(log10(maxab)~log10(sapply(lapply(dat,"[[", 2), function(x) nrow(x)))))
df1 <- as.data.frame(t(sapply(lapply(dat,"[[", 2), function(x) colMeans(x))))
df1$N <- sapply(lapply(dat,"[[", 2), function(x) nrow(x))
df1$ma <- maxab

fit2 <- lm(log10(ma)~log10(N)+K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, data = df1, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2)
sapply(lapply(dat,"[[", 2), function(x) order(x$ab, decreasing = T))

fit2 <- lm(log10(ma)~log10(N)+nComp+nMut++nPred+nAmens+nComm, data = df1, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2)

fit2 <- lm(log10(ma)~log10(N)+d.tot+cc.w+nComp+CompIn+nMut+nPred+PredOut+nAmens+CommIn, data = df1, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2)

topfive <- (allord %in% 1:5)
topten <- (allord %in% 1:10)

dat.fin$t5 <- topfive
dat.fin$t10 <- topten

gt5 <- glm(t5~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, family = "binomial", data = dat.fin, x = F, y = F, model = F, na.action = "na.fail")
summary(gt5)
DAAG::cv.binary(gt5)
gt10 <- glm(t10~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, family = "binomial", data = dat.fin, x = F, y = F, model = F, na.action = "na.fail")
summary(gt10)
DAAG::cv.binary(gt10)

library(MASS)
ldafit <- lda(t10~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, data = dat.fin)


#PCoAres <- vegan::capscale(dplyr::select(dat.fin, K2, bet.w, d.tot, cc.w, apl.w.mu, nComp:CommOut)~1, distance = "bray")
library(ape)
pc1 <- princomp(dplyr::select(dat.fin, K2, bet.w, d.tot, cc.w, apl.w.mu, nComp:CommOut))
loadings(pc1)
dmat <- dist(dplyr::select(dat.fin, K2, bet.w, d.tot, cc.w, apl.w.mu, nComp:CommOut), method = "euclidean")
