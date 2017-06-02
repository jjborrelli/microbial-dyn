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

fit2 <- MASS::lda(spp~K2+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut+allIn+allOut, data = rbind(dat.init,dat2.init), CV = TRUE)
tab <- table(factor(c(dat.init$spp,dat2.init$spp)), fit2$class)
conCV1 <- rbind(tab[1, ]/sum(tab[1, ]), tab[2, ]/sum(tab[2, ]))
dimnames(conCV1) <- list(Actual = c("No", "Yes"), "Predicted (cv)" = c("No","Yes"))
print(round(conCV1, 3))


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
maxab2 <- sapply(lapply(dat2,"[[", 2), function(x) max(x$ab/sum(x$ab)))

## 
fz1 <- lapply(lapply(dat, "[[", 2), function(x) fzmod(x$ab/sum(x$ab)))
rbfz <- do.call(rbind, fz1)
fz2 <- lapply(lapply(dat2, "[[", 2), function(x) fzmod(x$ab/sum(x$ab)))
rbfz2 <- do.call(rbind, fz2)

df1 <- as.data.frame(t(sapply(lapply(dat,"[[", 2), function(x) colMeans(x))))
df1.1 <- as.data.frame(t(sapply(lapply(dat2,"[[", 2), function(x) colMeans(x))))
df2 <- as.data.frame(apply(df1, 2, function(x){(x-mean(x))/sd(x)}))
df2.1 <- as.data.frame(apply(df1.1, 2, function(x){(x-mean(x))/sd(x)}))

df1$N <- sapply(lapply(dat,"[[", 2), function(x) nrow(x))
df1$ma <- maxab
df1$fz <- rbfz[,2]
df1.1$N <- sapply(lapply(dat2,"[[", 2), function(x) nrow(x))
df1.1$ma <- maxab2
df1.1$fz <- rbfz2[,2]

df2$N <- sapply(lapply(dat,"[[", 2), function(x) nrow(x))
df2$ma <- maxab
df2$fz <- rbfz[,2]
df2.1$N <- sapply(lapply(dat2,"[[", 2), function(x) nrow(x))
df2.1$ma <- maxab2
df2.1$fz <- rbfz2[,2]

#########################
# MODEL MAXIMUM ABUNDANCE
## ALL PRED
fit2ma <- lm(log10(ma)~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, data = rbind(df2, df2.1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit2ma)
AIC(fit2ma)
## WITHOUT COVARYING INT STRENGTHS (IN vs OUT)
fit2ma <- lm(log10(ma)~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut, data = rbind(df2, df2.1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit2ma)
AIC(fit2ma)
## ONLY NETWORK STRUCTURE
fit2ma <- lm(log10(ma)~K2+bet.w+d.tot+cc.w+apl.w.mu, data = rbind(df2, df2.1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit2ma)
AIC(fit2ma)
## ONLY INTERACTIONS
fit2ma <- lm(log10(ma)~nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, data = rbind(df2, df2.1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit2ma)
AIC(fit2ma)
## ONLY INTERACTIONS W/O COVARYING INT STRENGTHS
fit2ma <- lm(log10(ma)~nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut, data = rbind(df2, df2.1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit2ma)
AIC(fit2ma)
## ALL INT STRENGTHS AND NUM INTS
fit2ma <- lm(log10(ma)~d.tot+allOut, data = rbind(df2, df2.1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit2ma)
AIC(fit2ma)
## ONLY COMP AND MUT
fit2ma <- lm(log10(ma)~nComp+CompOut+nMut+MutOut, data = rbind(df2, df2.1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit2ma)
AIC(fit2ma)

#########################
# MODEL NUMBER OF SPECIES
## ALL PRED
fit2N <- lm(N~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, data = rbind(df2, df2.1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit2N)
AIC(fit2N)
## WITHOUT COVARYING INT STRENGTHS (IN vs OUT)
fit2N <- lm(N~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut, data = rbind(df2, df2.1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit2N)
AIC(fit2N)
## ONLY NETWORK STRUCTURE
fit2N <- lm(N~K2+bet.w+d.tot+cc.w+apl.w.mu, data = rbind(df2, df2.1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit2N)
AIC(fit2N)
## ONLY INTERACTIONS
fit2N <- lm(N~nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, data = rbind(df2, df2.1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit2N)
AIC(fit2N)
## ONLY INTERACTIONS W/O COVARYING INT STRENGTHS
fit2N <- lm(N~nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut, data = rbind(df2, df2.1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit2N)
AIC(fit2N)
## ALL INT STRENGTHS AND NUM INTS
fit2N <- lm(N~d.tot+allIn+allOut, data = rbind(df2, df2.1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit2N)
AIC(fit2N)
## ONLY COMP AND MUT
fit2N <- lm(N~nComp+CompOut+nMut+MutOut, data = rbind(df2, df2.1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit2N)
AIC(fit2N)


###########################
# MODEL RANK ABUNDANCE DIST
## ALL PRED
fit3 <- lm(fz~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, data = rbind(df2, df2.1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit3)
AIC(fit3)
## WITHOUT COVARYING INT STRENGTHS (IN vs OUT)
fit3 <- lm(fz~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut, data = rbind(df2, df2.1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit3)
AIC(fit3)
## ONLY NETWORK STRUCTURE
fit3 <- lm(fz~K2+bet.w+d.tot+cc.w+apl.w.mu, data = rbind(df2, df2.1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit3)
AIC(fit3)
## ONLY INTERACTIONS
fit3 <- lm(fz~nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, data = rbind(df2, df2.1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit3)
AIC(fit3)
## ONLY INTERACTIONS W/O COVARYING INT STRENGTHS
fit3 <- lm(fz~nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut, data = rbind(df2, df2.1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit3)
AIC(fit3)
## ALL INT STRENGTHS AND NUM INTS
fit3 <- lm(fz~d.tot+allIn+allOut, data = rbind(df2, df2.1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit3)
AIC(fit3)
## ONLY COMP AND MUT
fit3 <- lm(fz~nComp+CompOut+nMut+MutOut, data = rbind(df2, df2.1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit3)
AIC(fit3)


fit4 <- (lm(c(rbfz[,2],rbfz2[,2])~fit2ma$fitted.values+fit2N$fitted.values))
summary(fit4)
AIC(fit4)
###########################################################################
###########################################################################


df1 <- as.data.frame(t(sapply(lapply(dat,"[[", 1), function(x) colMeans(x[x$spp == 1,]))))
df1.1 <- as.data.frame(t(sapply(lapply(dat2,"[[", 1), function(x) colMeans(x[x$spp == 1,]))))
df2 <- as.data.frame(apply(df1, 2, function(x){(x-mean(x))/sd(x)}))
df2.1 <- as.data.frame(apply(df1.1, 2, function(x){(x-mean(x))/sd(x)}))

df1.1$N <- sapply(lapply(dat2,"[[", 2), function(x) nrow(x))
df1.1$ma <- maxab2
df1.1$fz <- rbfz2[,2]

df2$N <- sapply(lapply(dat,"[[", 2), function(x) nrow(x))
df2$ma <- maxab
df2$fz <- rbfz[,2]
df2.1$N <- sapply(lapply(dat2,"[[", 2), function(x) nrow(x))
df2.1$ma <- maxab2
df2.1$fz <- rbfz2[,2]
