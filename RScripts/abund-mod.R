library(data.table)
library(rpart)
library(rpart.plot)
library(randomForest)
library(sads)
library(lme4)

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
datrel$com <- dat.fin$com

datz <- rbindlist(lapply(lapply(dat, "[[", 2), function(x){for(i in 1:ncol(x)){x[,i] <- (x[,i]-mean(x[,i]))/sd(x[,i]); x[,i][is.nan(x[,i])] <- 0};return(x)}))
datz$com <- dat.fin$com
datz$ab <- dat.fin$ab
datz$rab <- dat.fin$rab
dat2rel <- rbindlist(lapply(lapply(dat2, "[[", 2), function(x){for(i in 1:ncol(x)){x[,i] <- x[,i]/sum(x[,i]); x[,i][is.nan(x[,i])] <- 0};return(x)}))

plot(dat.init$spp~dat.init$ab.i)

fit1 <- lmer(rab~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut+(1 | com), data = datz, na.action = "na.fail")
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


maxab <- sapply(lapply(dat,"[[", 2), function(x) max(x$ab/sum(x$ab)))
summary(maxab)
maxab2 <- sapply(lapply(dat,"[[", 2), function(x) max(x$ab[-which.max(x$ab)]))

summary(lm(log10(maxab)~log10(sapply(lapply(dat,"[[", 2), function(x) nrow(x)))))
df1 <- as.data.frame(t(sapply(lapply(dat,"[[", 2), function(x) colMeans(x))))
df1$N <- sapply(lapply(dat,"[[", 2), function(x) nrow(x))
df1$ma <- maxab

fit2 <- lm(cbind(ma, N)~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, data = df1, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2)
sapply(lapply(dat,"[[", 2), function(x) order(x$ab, decreasing = T))

fit2 <- lm(log10(ma)~log10(N)+nComp+nMut++nPred+nAmens+nComm, data = df1, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2)

fit2 <- lm(log10(ma)~log10(N)+d.tot+cc.w+nComp+CompIn+nMut+nPred+PredOut+nAmens+CommIn, data = df1, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2)


allord <- unlist(lapply(lapply(dat,"[[", 2), function(x) order(x$ab, decreasing = T)))
topthree <- (allord %in% 1:3)
topfive <- (allord %in% 1:5)
topten <- (allord %in% 1:10)

dat.fin$t3 <- topthree
dat.fin$t5 <- topfive
dat.fin$t10 <- topten


gt3 <- glm(dat.fin$t5~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, family = "binomial", data = dat.fin, x = F, y = F, model = F, na.action = "na.fail")
summary(gt3)

gt5 <- glm(dat.fin$t5~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, family = "binomial", data = dat.fin, x = F, y = F, model = F, na.action = "na.fail")
summary(gt5)

gt10 <- glm(dat.fin$t10~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, family = "binomial", data = dat.fin, x = F, y = F, model = F, na.action = "na.fail")
summary(gt10)

library(MASS)
ldafit3 <- lda(t3~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, data = dat.fin, CV = TRUE)
tab <- table(factor(dat.fin$t3), ldafit3$class)
conCV1 <- rbind(tab[1, ]/sum(tab[1, ]), tab[2, ]/sum(tab[2, ]))
dimnames(conCV1) <- list(Actual = c("No", "Yes"), "Predicted (cv)" = c("No","Yes"))
print(round(conCV1, 3))

ldafit5 <- lda(factor(t5)~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, data = dat.fin, CV = TRUE)
tab <- table(factor(dat.fin$t5), ldafit5$class)
conCV1 <- rbind(tab[1, ]/sum(tab[1, ]), tab[2, ]/sum(tab[2, ]))
dimnames(conCV1) <- list(Actual = c("No", "Yes"), "Predicted (cv)" = c("No","Yes"))
print(round(conCV1, 3))

ldafit10 <- lda(factor(t10)~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, data = dat.fin, CV = TRUE)
tab <- table(factor(dat.fin$t10), ldafit10$class)
conCV1 <- rbind(tab[1, ]/sum(tab[1, ]), tab[2, ]/sum(tab[2, ]))
dimnames(conCV1) <- list(Actual = c("No", "Yes"), "Predicted (cv)" = c("No","Yes"))
print(round(conCV1, 3))


cbind(coefficients(ldafit5),coefficients(ldafit3),coefficients(ldafit10))


rf3 <- randomForest(factor(t3)~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, 
                    data = datrel, mtry = 15, ntree = 1000)


###########################################################################
###########################################################################
maxab <- sapply(lapply(dat,"[[", 2), function(x) max(x$ab))
summary(maxab)

## 
fz1 <- lapply(lapply(dat, "[[", 2), function(x) fzmod(get_abundvec(x$ab/sum(x$ab), 2000)))
rbfz <- do.call(rbind, fz1)

df1 <- as.data.frame(t(sapply(lapply(dat,"[[", 2), function(x) colMeans(x))))
df2 <- as.data.frame(apply(df1, 2, function(x){(x-mean(x))/sd(x)}))

df1$N <- sapply(lapply(dat,"[[", 2), function(x) nrow(x))
df1$ma <- maxab
df1$fz <- rbfz[,2]

df2$N <- sapply(lapply(dat,"[[", 2), function(x) nrow(x))
df2$ma <- maxab
df2$fz <- rbfz[,2]


fit2 <- lm(cbind(log10(ma), N)~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, data = df2, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2)

fit3 <- lm(fz~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, data = df2, x = F, y = F, model = F, na.action = "na.fail")
summary(fit3)
AIC(fit3)
fit4 <- (lm(rbfz[,2]~fit2$fitted.values))
summary(fit4)
AIC(fit4)
###########################################################################
###########################################################################


#PCoAres <- vegan::capscale(dplyr::select(dat.fin, K2, bet.w, d.tot, cc.w, apl.w.mu, nComp:CommOut)~1, distance = "bray")
library(ape)
pc1 <- princomp(dplyr::select(dat.fin, K2, bet.w, d.tot, cc.w, apl.w.mu, nComp:CommOut))
loadings(pc1)

library(vegan)
nmds1 <- metaMDS(dplyr::select(dat.fin, K2, bet.w, d.tot, cc.w, apl.w.mu, nComp:CommOut), distance = "bray")
