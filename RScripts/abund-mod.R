library(data.table)
library(rpart)
library(rpart.plot)
library(randomForest)
library(sads)
library(lme4)

dat <- readRDS("~/Documents/AbundData/res.rds")
dat2 <- readRDS("~/Documents/AbundData/res2.rds")
dat3 <- readRDS("~/Documents/AbundData/res3.rds")
dat4 <- readRDS("~/Documents/AbundData/res4.rds")
dat5 <- readRDS("~/Documents/AbundData/res5.rds")
dat6 <- readRDS("~/Documents/AbundData/res6.rds")

d6conn <- do.call(rbind, lapply(dat6, "[[", 3))

dat.fin <- rbindlist(lapply(dat, "[[", 2))
dat.init <- rbindlist(lapply(dat, "[[", 1))
dat.init$com <- rep(1:1500, each = 600)
dat.init$rab <- unlist(lapply(lapply(dat, "[[", 1), function(x) x$ab/sum(x$ab)))
dat.init$rabi <- unlist(lapply(lapply(dat, "[[", 1), function(x) x$ab.i/sum(x$ab.i)))
dat.fin$com <- rep(1:1500, sapply(lapply(dat, "[[", 2), nrow))
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
dat2.fin$N <- rep(sapply(lapply(dat2, "[[", 2), nrow),sapply(lapply(dat2, "[[", 2), nrow))


datrel <- rbindlist(lapply(lapply(dat, "[[", 2), function(x){for(i in 1:ncol(x)){x[,i] <- x[,i]/sum(x[,i]); x[,i][is.nan(x[,i])] <- 0};return(x)}))
datrel$com <- dat.fin$com

datz <- rbindlist(lapply(lapply(dat, "[[", 2), function(x){for(i in 1:ncol(x)){x[,i] <- (x[,i]-mean(x[,i]))/sd(x[,i]); x[,i][is.nan(x[,i])] <- 0};return(x)}))
datz$com <- dat.fin$com
datz$ab <- dat.fin$ab
datz$rab <- dat.fin$rab
dat2rel <- rbindlist(lapply(lapply(dat2, "[[", 2), function(x){for(i in 1:ncol(x)){x[,i] <- x[,i]/sum(x[,i]); x[,i][is.nan(x[,i])] <- 0};return(x)}))

plot(dat.init$spp~dat.init$ab.i)
fit1 <- lmer(log10(ab)~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut+(1 | com), data = dat.fin, na.action = "na.fail")
summary(fit1)


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
dat.comp <- readRDS("~/Documents/AbundData/rescomp1.rds")
macom <- sapply(lapply(dat.comp,"[[", 2), function(x) max(x$ab))
fzcom <- do.call(rbind,lapply(lapply(dat.comp, "[[", 2), function(x) fzmod(get_abundvec(x$ab/sum(x$ab), 2000))))

###########################################################################
###########################################################################

dat <- readRDS("~/Documents/AbundData/res.rds")
dat2 <- readRDS("~/Documents/AbundData/res2.rds")
dat3 <- readRDS("~/Documents/AbundData/res3.rds")
dat4 <- readRDS("~/Documents/AbundData/res4.rds")
dat5 <- readRDS("~/Documents/AbundData/res5.rds")
dat6 <- readRDS("~/Documents/AbundData/res6.rds")


###########################################################################

maxab <- sapply(lapply(dat,"[[", 2), function(x) max(x$ab))
maxab2 <- sapply(lapply(dat2,"[[", 2), function(x) max(x$ab))
maxab3 <- sapply(lapply(dat3,"[[", 2), function(x) max(x$ab))
maxab4 <- sapply(lapply(dat4,"[[", 2), function(x) max(x$ab))
maxab5 <- sapply(lapply(dat5,"[[", 2), function(x) max(x$ab))
maxab6 <- sapply(lapply(dat6,"[[", 2), function(x) max(x$ab))


ga1 <- lapply(lapply(dat, "[[", 2), function(x) (get_abundvec(x$ab/sum(x$ab), 2000)))
ga2 <- lapply(lapply(dat2, "[[", 2), function(x) (get_abundvec(x$ab/sum(x$ab), 2000)))
ga3 <- lapply(lapply(dat3, "[[", 2), function(x) (get_abundvec(x$ab/sum(x$ab), 2000)))
ga4 <- lapply(lapply(dat4, "[[", 2), function(x) (get_abundvec(x$ab/sum(x$ab), 2000)))
ga5 <- lapply(lapply(dat5, "[[", 2), function(x) (get_abundvec(x$ab/sum(x$ab), 2000)))
ga6 <- lapply(lapply(dat6, "[[", 2), function(x) (get_abundvec(x$ab/sum(x$ab), 2000)))

sk1 <- sapply(ga1, function(x) sum((x-mean(x))^3/length(x))/sd(x)^3)
sk2 <- sapply(ga2, function(x) sum((x-mean(x))^3/length(x))/sd(x)^3)
sk3 <- sapply(ga3, function(x) sum((x-mean(x))^3/length(x))/sd(x)^3)
sk4 <- sapply(ga4, function(x) sum((x-mean(x))^3/length(x))/sd(x)^3)
sk5 <- sapply(ga5, function(x) sum((x-mean(x))^3/length(x))/sd(x)^3)
sk6 <- sapply(ga6, function(x) sum((x-mean(x))^3/length(x))/sd(x)^3)


## 
rbfz <- lapply(ga1, function(x) fzmod(x))
#rbfz <- lapply(lapply(dat, "[[", 2), function(x) fzmod((x$ab/sum(x$ab))))
rbfz <- do.call(rbind, rbfz)
rbfz2 <- lapply(ga2, function(x) fzmod(x))
#rbfz2 <- lapply(lapply(dat2, "[[", 2), function(x) fzmod((x$ab/sum(x$ab))))
rbfz2 <- do.call(rbind, rbfz2)
rbfz3 <- lapply(ga3, function(x) fzmod(x))
#rbfz3 <- lapply(lapply(dat3, "[[", 2), function(x) fzmod((x$ab/sum(x$ab))))
rbfz3 <- do.call(rbind, rbfz3)
rbfz4 <- lapply(ga4, function(x) fzmod(x))
#rbfz4 <- lapply(lapply(dat4, "[[", 2), function(x) fzmod((x$ab/sum(x$ab))))
rbfz4 <- do.call(rbind, rbfz4)
rbfz5 <- lapply(ga5, function(x) fzmod(x))
#rbfz5 <- lapply(lapply(dat5, "[[", 2), function(x) fzmod((x$ab/sum(x$ab))))
rbfz5 <- do.call(rbind, rbfz5)
rbfz6 <- lapply(ga6, function(x) fzmod(x))
#rbfz6 <- lapply(lapply(dat6, "[[", 2), function(x) fzmod((x$ab/sum(x$ab))))
rbfz6 <- do.call(rbind, rbfz6)

#points(rbfz[,1:2], pch = 20)
#plot(rbfz2[,1:2], pch = 20, col = "blue")

df1.1sd <- as.data.frame(t(sapply(lapply(dat,"[[", 2), function(x) apply(dplyr::select(x, gr:allOut), 2, sd))))
colnames(df1.1sd) <- paste0(colnames(df1.1sd), ".sd")
df1.2sd <- as.data.frame(t(sapply(lapply(dat2,"[[", 2), function(x) apply(dplyr::select(x, gr:allOut), 2, sd))))
colnames(df1.2sd) <- paste0(colnames(df1.2sd), ".sd")
df1.3sd <- as.data.frame(t(sapply(lapply(dat3,"[[", 2), function(x) apply(dplyr::select(x, gr:allOut), 2, sd))))
colnames(df1.3sd) <- paste0(colnames(df1.3sd), ".sd")
df1.4sd <- as.data.frame(t(sapply(lapply(dat4,"[[", 2), function(x) apply(dplyr::select(x, gr:allOut), 2, sd))))
colnames(df1.4sd) <- paste0(colnames(df1.4sd), ".sd")
df1.5sd <- as.data.frame(t(sapply(lapply(dat5,"[[", 2), function(x) apply(dplyr::select(x, gr:allOut), 2, sd))))
colnames(df1.5sd) <- paste0(colnames(df1.5sd), ".sd")
df1.6sd <- as.data.frame(t(sapply(lapply(dat6,"[[", 2), function(x) apply(dplyr::select(x, gr:allOut), 2, sd))))
colnames(df1.6sd) <- paste0(colnames(df1.6sd), ".sd")


df1.1 <- as.data.frame(t(sapply(lapply(dat,"[[", 2), function(x) colMeans(x))))
df1.2 <- as.data.frame(t(sapply(lapply(dat2,"[[", 2), function(x) colMeans(x))))
df1.3 <- as.data.frame(t(sapply(lapply(dat3,"[[", 2), function(x) colMeans(x))))
df1.4 <- as.data.frame(t(sapply(lapply(dat4,"[[", 2), function(x) colMeans(x))))
df1.5 <- as.data.frame(t(sapply(lapply(dat5,"[[", 2), function(x) colMeans(x))))
df1.6 <- as.data.frame(t(sapply(lapply(dat6,"[[", 2), function(x) colMeans(x))))
df2.1 <- as.data.frame(apply(df1.1, 2, function(x){(x-mean(x))/sd(x)}))
df2.2 <- as.data.frame(apply(df1.2, 2, function(x){(x-mean(x))/sd(x)}))
df2.3 <- as.data.frame(apply(df1.3, 2, function(x){(x-mean(x))/sd(x)}))
df2.4 <- as.data.frame(apply(df1.4, 2, function(x){(x-mean(x))/sd(x)}))
df2.5 <- as.data.frame(apply(df1.5, 2, function(x){(x-mean(x))/sd(x)}))
df2.6 <- as.data.frame(apply(df1.6, 2, function(x){(x-mean(x))/sd(x)}))

df1.1$N <- rbfz[,1]#sapply(lapply(dat,"[[", 2), function(x) nrow(x))
df1.1$ma <- maxab
df1.1$fz <- rbfz[,2]
df1.2$N <- rbfz2[,1]#sapply(lapply(dat2,"[[", 2), function(x) nrow(x))
df1.2$ma <- maxab2
df1.2$fz <- rbfz2[,2]
df1.3$N <- rbfz3[,1]#sapply(lapply(dat3,"[[", 2), function(x) nrow(x))
df1.3$ma <- maxab3
df1.3$fz <- rbfz3[,2]
df1.4$N <- rbfz4[,1]#sapply(lapply(dat4,"[[", 2), function(x) nrow(x))
df1.4$ma <- maxab4
df1.4$fz <- rbfz4[,2]
df1.5$N <- rbfz5[,1]#sapply(lapply(dat5,"[[", 2), function(x) nrow(x))
df1.5$ma <- maxab5
df1.5$fz <- rbfz5[,2]
df1.6$N <- rbfz6[,1]#sapply(lapply(dat6,"[[", 2), function(x) nrow(x))
df1.6$ma <- maxab6
df1.6$fz <- rbfz6[,2]


df2.1$N <- rbfz[,1]#sapply(lapply(dat,"[[", 2), function(x) nrow(x))
df2.1$ma <- maxab
df2.1$fz <- rbfz[,2]
df2.2$N <- rbfz2[,1]#sapply(lapply(dat2,"[[", 2), function(x) nrow(x))
df2.2$ma <- maxab2
df2.2$fz <- rbfz2[,2]
df2.3$N <- rbfz3[,1]#sapply(lapply(dat3,"[[", 2), function(x) nrow(x))
df2.3$ma <- maxab3
df2.3$fz <- rbfz3[,2]
df2.4$N <- rbfz4[,1]#sapply(lapply(dat4,"[[", 2), function(x) nrow(x))
df2.4$ma <- maxab4
df2.4$fz <- rbfz4[,2]
df2.5$N <- rbfz5[,1]#sapply(lapply(dat5,"[[", 2), function(x) nrow(x))
df2.5$ma <- maxab5
df2.5$fz <- rbfz5[,2]
df2.6$N <- rbfz6[,1]#sapply(lapply(dat6,"[[", 2), function(x) nrow(x))
df2.6$ma <- maxab6
df2.6$fz <- rbfz6[,2]

alldat <- rbindlist(list(df2.1, df2.2, df2.3, df2.4, df2.5, df2.6))
alldat$typ <- factor(rep(c("rand", "hub"), c(3000, 4500)))
alldat2 <- rbindlist(list(df1.1, df1.2, df1.3, df1.4, df1.5, df1.6))
alldat2$typ <- factor(rep(c("rand", "hub"), c(3000, 4500)))
alldat2$r2 <- c(rbfz[,"r2"],rbfz2[,"r2"], rbfz3[,"r2"],rbfz4[,"r2"],rbfz5[,"r2"],rbfz6[,"r2"])
alldat$N2 <- ((alldat$N-mean(alldat$N))/sd(alldat$N))
alldat2$sk <- c(sk1, sk2, sk3, sk4, sk5, sk6)
alldat$sk <- c(sk1, sk2, sk3, sk4, sk5, sk6)
#########################
# MODEL MAXIMUM ABUNDANCE
fit2ma <- lm(log10(ma)~N2+K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut+typ, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
fitmadf <- data.frame(full1 = c(fit2ma$coefficients, allOut = NA))
coefnames <- c(names(fit2ma$coefficients), "allOut")
fitmap <- data.frame(full1 = c(summary(fit2ma)$coefficients[,4] <= 0.05, allOut = NA))
mafit <- data.frame(full1 = c(NA, NA), full2 = c(NA, NA), struct = c(NA, NA), ints = c(NA, NA), ints2 = c(NA, NA), degstr = c(NA, NA), commut = c(NA, NA))
## ALL PRED
fit2ma <- lm(log10(ma)~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut+typ, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2ma) # r2 = .75
mafit$full1 <- c(AIC(fit2ma), summary(fit2ma)$r.squared) # -1909
fitmadf$full1[coefnames %in% names(fit2ma$coefficients)] <- fit2ma$coefficients
fitmap$full1[coefnames %in% names(fit2ma$coefficients)] <- summary(fit2ma)$coefficients[,4] <= 0.05
## WITHOUT COVARYING INT STRENGTHS (IN vs OUT)
fit2ma <- lm(log10(ma)~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut+typ, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2ma) # r2 = ..74
mafit$full2 <- c(AIC(fit2ma), summary(fit2ma)$r.squared) # -1876
fitmadf$full2[coefnames %in% names(fit2ma$coefficients)] <- fit2ma$coefficients
fitmap$full2[coefnames %in% names(fit2ma$coefficients)] <- summary(fit2ma)$coefficients[,4] <= 0.05
## ONLY NETWORK STRUCTURE
fit2ma <- lm(log10(ma)~K2+bet.w+d.tot+cc.w+apl.w.mu+typ, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2ma) # r2 = .30 
mafit$struct <- c(AIC(fit2ma), summary(fit2ma)$r.squared)  # 103
fitmadf$struct[coefnames %in% names(fit2ma$coefficients)] <- fit2ma$coefficients
fitmap$struct[coefnames %in% names(fit2ma$coefficients)] <- summary(fit2ma)$coefficients[,4] <= 0.05
## ONLY INTERACTIONS
fit2ma <- lm(log10(ma)~nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut+typ, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2ma) # r2 = .59
mafit$ints <- c(AIC(fit2ma), summary(fit2ma)$r.squared)  # -982
fitmadf$ints[coefnames %in% names(fit2ma$coefficients)] <- fit2ma$coefficients
fitmap$ints[coefnames %in% names(fit2ma$coefficients)] <- summary(fit2ma)$coefficients[,4] <= 0.05
## ONLY INTERACTIONS W/O COVARYING INT STRENGTHS
fit2ma <- lm(log10(ma)~nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut+typ, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2ma) # r2 = .56
mafit$ints2 <- c(AIC(fit2ma), summary(fit2ma)$r.squared) # -818
fitmadf$ints2[coefnames %in% names(fit2ma$coefficients)] <- fit2ma$coefficients
fitmap$ints2[coefnames %in% names(fit2ma$coefficients)] <- summary(fit2ma)$coefficients[,4] <= 0.05
## ALL INT STRENGTHS AND NUM INTS
fit2ma <- lm(log10(ma)~d.tot+allOut+typ, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2ma) # r2 = .37
mafit$degstr <- c(AIC(fit2ma), summary(fit2ma)$r.squared)  # -116
fitmadf$degstr[coefnames %in% names(fit2ma$coefficients)] <- fit2ma$coefficients
fitmap$degstr[coefnames %in% names(fit2ma$coefficients)] <- summary(fit2ma)$coefficients[,4] <= 0.05
## ONLY COMP AND MUT
fit2ma <- lm(log10(ma)~nComp+CompOut+nMut+MutOut+typ, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2ma) # r2 = .34
mafit$commut <- c(AIC(fit2ma), summary(fit2ma)$r.squared)  # -10
fitmadf$commut[coefnames %in% names(fit2ma$coefficients)] <- fit2ma$coefficients
fitmap$commut[coefnames %in% names(fit2ma$coefficients)] <- summary(fit2ma)$coefficients[,4] <= 0.05


#########################
# MODEL NUMBER OF SPECIES
fitNdf <- fitmadf
fitNp <- fitmap
Nfit <- data.frame(full1 = c(NA, NA), full2 = c(NA, NA), struct = c(NA, NA), ints = c(NA, NA), ints2 = c(NA, NA), degstr = c(NA, NA), commut = c(NA, NA))
## ALL PRED
fit2N <- lm(N~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut+typ, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2N) # .96
Nfit$full1 <- c(AIC(fit2N), summary(fit2N)$r.squared)  # 18123
fitNdf$full1[coefnames %in% names(fit2N$coefficients)] <- fit2N$coefficients
fitNp$full1[coefnames %in% names(fit2N$coefficients)] <- summary(fit2N)$coefficients[,4] <= 0.05
## WITHOUT COVARYING INT STRENGTHS (IN vs OUT)
fit2N <- lm(N~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut+typ, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2N) # .96
Nfit$full2 <- c(AIC(fit2N), summary(fit2N)$r.squared)  # 18116
fitNdf$full2[coefnames %in% names(fit2N$coefficients)] <- fit2N$coefficients
fitNp$full2[coefnames %in% names(fit2N$coefficients)] <- summary(fit2N)$coefficients[,4] <= 0.05
## ONLY NETWORK STRUCTURE
fit2N <- lm(N~K2+bet.w+d.tot+cc.w+apl.w.mu+typ, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2N) # .96
Nfit$struct <- c(AIC(fit2N), summary(fit2N)$r.squared)  # 18144
fitNdf$struct[coefnames %in% names(fit2N$coefficients)] <- fit2N$coefficients
fitNp$struct[coefnames %in% names(fit2N$coefficients)] <- summary(fit2N)$coefficients[,4] <= 0.05
## ONLY INTERACTIONS
fit2N <- lm(N~nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut+typ, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2N) # .63
Nfit$ints <- c(AIC(fit2N), summary(fit2N)$r.squared)  # 22393
fitNdf$ints[coefnames %in% names(fit2N$coefficients)] <- fit2N$coefficients
fitNp$ints[coefnames %in% names(fit2N$coefficients)] <- summary(fit2N)$coefficients[,4] <= 0.05
## ONLY INTERACTIONS W/O COVARYING INT STRENGTHS
fit2N <- lm(N~nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut+typ, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2N) # .58
Nfit$ints2 <- c(AIC(fit2N), summary(fit2N)$r.squared) # 22595
fitNdf$ints2[coefnames %in% names(fit2N$coefficients)] <- fit2N$coefficients
fitNp$ints2[coefnames %in% names(fit2N$coefficients)] <- summary(fit2N)$coefficients[,4] <= 0.05
## ALL INT STRENGTHS AND NUM INTS
fit2N <- lm(N~d.tot+allOut+typ, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2N) # .25
Nfit$degstr <- c(AIC(fit2N), summary(fit2N)$r.squared) # 23772
fitNdf$degstr[coefnames %in% names(fit2N$coefficients)] <- fit2N$coefficients
fitNp$degstr[coefnames %in% names(fit2N$coefficients)] <- summary(fit2N)$coefficients[,4] <= 0.05
## ONLY COMP AND MUT
fit2N <- lm(N~nComp+CompOut+nMut+MutOut+typ, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2N) # .41
Nfit$commut <- c(AIC(fit2N), summary(fit2N)$r.squared)  # 23264
fitNdf$commut[coefnames %in% names(fit2N$coefficients)] <- fit2N$coefficients
fitNp$commut[coefnames %in% names(fit2N$coefficients)] <- summary(fit2N)$coefficients[,4] <= 0.05



###########################
# MODEL RANK ABUNDANCE DIST
fitsdf <- fitmadf
sfitp <- fitmap
sfit <- data.frame(full1 = c(NA, NA), full2 = c(NA, NA), struct = c(NA, NA), ints = c(NA, NA), ints2 = c(NA, NA), degstr = c(NA, NA), commut = c(NA, NA))
## ALL PRED
fit3 <- lm(fz~N2, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit3) # .83
sfit$full1 <- c(AIC(fit3), summary(fit3)$r.squared) # -5350
fitsdf$full1[coefnames %in% names(fit3$coefficients)] <- fit3$coefficients
sfitp$full1[coefnames %in% names(fit3$coefficients)] <- summary(fit3)$coefficients[,4] <= 0.05
## WITHOUT COVARYING INT STRENGTHS (IN vs OUT)
fit3 <- lm((fz)~K2+bet.w+d.tot+cc.w+apl.w.mu+nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut+typ, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit3) # .83
sfit$full2 <- c(AIC(fit3), summary(fit3)$r.squared)# -5341
fitsdf$full2[coefnames %in% names(fit3$coefficients)] <- fit3$coefficients
sfitp$full2[coefnames %in% names(fit3$coefficients)] <- summary(fit3)$coefficients[,4] <= 0.05
## ONLY NETWORK STRUCTURE
fit3 <- lm(log(fz)~bet.w+d.tot+cc.w+mod.uw+typ, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit3) # .80
fit3.2 <- lm(log(fz)~K2+bet.uw+d.tot+cc.uw+apl.uw.mu+typ, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit3.2)
sfit$struct <- c(AIC(fit3), summary(fit3)$r.squared) # -5030
fitsdf$struct[coefnames %in% names(fit3$coefficients)] <- fit3$coefficients
sfitp$struct[coefnames %in% names(fit3$coefficients)] <- summary(fit3)$coefficients[,4] <= 0.05
## ONLY INTERACTIONS
fit3 <- lm(fz~nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit3) # .55
sfit$ints <- c(AIC(fit3), summary(fit3)$r.squared)#  -3397
fitsdf$ints[coefnames %in% names(fit3$coefficients)] <- fit3$coefficients
sfitp$ints[coefnames %in% names(fit3$coefficients)] <- summary(fit3)$coefficients[,4] <= 0.05
## ONLY INTERACTIONS W/O COVARYING INT STRENGTHS
fit3 <- lm(fz~nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit3) # .52
sfit$ints2 <- c(AIC(fit3), summary(fit3)$r.squared) # -3287
fitsdf$ints2[coefnames %in% names(fit3$coefficients)] <- fit3$coefficients
sfitp$ints2[coefnames %in% names(fit3$coefficients)] <- summary(fit3)$coefficients[,4] <= 0.05
## ALL INT STRENGTHS AND NUM INTS
fit3 <- lm(fz~d.tot+allOut+typ, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit3) # .22
sfit$degstr <- c(AIC(fit3), summary(fit3)$r.squared) # -2320
fitsdf$degstr[coefnames %in% names(fit3$coefficients)] <- fit3$coefficients
sfitp$degstr[coefnames %in% names(fit3$coefficients)] <- summary(fit3)$coefficients[,4] <= 0.05
## ONLY COMP AND MUT
fit3 <- lm(fz~nComp+CompOut+nMut+MutOut, data = alldat, x = F, y = F, model = F, na.action = "na.fail")
summary(fit3) # .43
sfit$commut <- c(AIC(fit3), summary(fit3)$r.squared) # -2935
fitsdf$commut[coefnames %in% names(fit3$coefficients)] <- fit3$coefficients
sfitp$commut[coefnames %in% names(fit3$coefficients)] <- summary(fit3)$coefficients[,4] <= 0.05

# fit2ma and fit2N second fits better than first
fit4 <- (lm(alldat$fz~fit2ma$fitted.values+fit2N$fitted.values+alldat$typ)) 
summary(fit4)
AIC(fit4)
###########################################################################
###########################################################################


fz1 <- lapply(lapply(dat, "[[", 2), function(x) fzmod(get_abundvec(x$ab/sum(x$ab), 2000)))
fz1 <- do.call(rbind, fz1)
fz2 <- lapply(lapply(dat2, "[[", 2), function(x) fzmod(get_abundvec(x$ab/sum(x$ab), 2000)))
fz2 <- do.call(rbind, fz2)
fz4 <- lapply(lapply(dat4, "[[", 2), function(x) fzmod(get_abundvec(x$ab/sum(x$ab), 2000)))
fz4 <- do.call(rbind, fz4)
fz5 <- lapply(lapply(dat5, "[[", 2), function(x) fzmod(get_abundvec(x$ab/sum(x$ab), 2000)))
fz5 <- do.call(rbind, fz5)
fz6.1 <- lapply(lapply(dat6, "[[", 2), function(x) fzmod((x$ab/sum(x$ab))))
fz6.1 <- do.call(rbind, fz6.1)

points(rbindlist(list(fz1[,1:2], fz2[,1:2], fz4[,1:2], fz5[,1:2], fz6[,1:2])), col = "darkgrey", pch = 20)


ct <- c()
pv <-c()
conns <- unique(round(t(sapply(dat6, "[[", 3))[,1], 2))[1:23]
for(i in 1:length(conns)){
  corr <- cor.test(df1.6$nAmens[round(t(sapply(dat6, "[[", 3))[,1], 2) == conns[i]], rbfz6$s[round(t(sapply(dat6, "[[", 3))[,1], 2) == conns[i]])
  ct[i] <- corr$estimate
  pv[i] <- corr$p.value
  #plot(df1.6$nComp[round(t(sapply(dat6, "[[", 3))[,1], 2) == conns[i]], rbfz6$s[round(t(sapply(dat6, "[[", 3))[,1], 2) == conns[i]])
}

cbind(ct,pv <= 0.05)


###########################################################################
###########################################################################


#crt <- rpart(fz~K2+bet.w+d.in+d.out+cc.w+apl.w.mu+nComp+CompIn+CompOut+nMut+MutIn+MutOut+nPred+PredIn+PredOut+nAmens+AmensIn+AmensOut+nComm+CommIn+CommOut+typ, dat = alldat2[alldat2$r2 > 0.8])
crt <- rpart(fz~N+K2+bet.w+bet.uw+d.in+d.out+cc.w+cc.uw+apl.uw.mu+apl.w.mu+nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut+typ, dat = alldat[alldat2$r2 > 0.8 & alldat$N > 130,])
#prp(crt, extra = 1, main = "S")
#plotcp(crt)
prp(prune(crt, cp = crt$cptable[which.min(crt$cptable[,"xerror"]),"CP"]), extra = 1)
par(mfrow = c(1,2))
rsq.rpart(crt)
rsq.rpart(prune(crt, cp = crt$cptable[which.min(crt$cptable[,"xerror"]),"CP"]))
par(mfrow = c(1,1))

rf <- randomForest(fz~N+K2+bet.w+bet.uw+d.in+d.out+cc.w+cc.uw+apl.uw.mu+apl.w.mu+nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut+typ, dat = alldat2[alldat2$r2 > 0.8], importance = TRUE, mtry = 15, ntree = 2000)
varImpPlot(rf)

crt <- rpart(r2~N+K2+bet.w+bet.uw+d.in+d.out+cc.w+cc.uw+apl.uw.mu+apl.w.mu+nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut+typ, dat = alldat2)
prp(crt, extra = 1)
plotcp(crt)
prp(prune(crt, cp = crt$cptable[which.min(crt$cptable[,"xerror"]),"CP"]))
rsq.rpart(crt)

crt <- rpart(ma~N+K2+bet.w+bet.uw+d.in+d.out+cc.w+cc.uw+apl.uw.mu+apl.w.mu+nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut+typ, dat = alldat2[alldat2$r2 > 0.8], method = "anova")
prp(crt, extra = 1, main = "Max Ab")
plotcp(crt)
prp(prune(crt, cp = crt$cptable[which.min(crt$cptable[,"xerror"]),"CP"]))

rf <- randomForest(ma~N+K2+bet.w+bet.uw+d.in+d.out+cc.w+cc.uw+apl.uw.mu+apl.w.mu+nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut+typ, dat = alldat2[alldat2$r2 > 0.8], importance = TRUE, mtry = 15, ntree = 2000)
varImpPlot(rf)


crt <- rpart((N > 130)~K2+bet.w+bet.uw+d.in+d.out+cc.w+cc.uw+apl.uw.mu+apl.w.mu+nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut+typ, dat = alldat[alldat2$r2 > 0.8], method= "class")
prp(crt, extra = 1, main = "N")
plotcp(crt)
prp(prune(crt, cp = crt$cptable[which.min(crt$cptable[,"xerror"]),"CP"]), extra = 1)
rf <- randomForest(N~K2+bet.w+bet.uw+d.in+d.out+cc.w+cc.uw+apl.uw.mu+apl.w.mu+nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut+typ, dat = alldat2[alldat2$r2 > 0.8], importance = TRUE, mtry = 15, ntree = 2000)
varImpPlot(rf)


plot(alldat2[alldat2$r2 > 0.5,]$N, alldat2[alldat2$r2 > 0.5,]$bet.uw)


########################################################################################
########################################################################################
ga <- lapply(lapply(dat3, "[[", 2), function(x) (get_abundvec(x$ab/sum(x$ab), 2000)))
rbfz <- sapply(ga, function(x) fitpoilog(x)@fullcoef)
rbfz2 <- lapply(ga, function(x) radfit(x))
rbfz6 <- sapply(ga2, function(x) fzmod(x)$s)
rbfz6.2 <- sapply(ga2, function(x) fitpoilog(x)@fullcoef)

cmplx <- sqrt(c(sapply(ga, length),sapply(ga2, length)) * c(sapply(dat3, "[[", 3)[1,],sapply(dat6, "[[", 3)[1,]))

sum((ga[[1]]-mean(ga[[1]]))^3/length(ga[[1]]))/sd(ga[[1]])^3

sk1 <- sapply(ga, function(x) sum((x-mean(x))^3/length(x))/sd(x)^3)
sk2 <- sapply(ga, function(x) sum((x-mean(x))^3/length(x))/sd(x)^3)
sk3 <- sapply(ga, function(x) sum((x-mean(x))^3/length(x))/sd(x)^3)
sk4 <- sapply(ga, function(x) sum((x-mean(x))^3/length(x))/sd(x)^3)
sk5 <- sapply(ga, function(x) sum((x-mean(x))^3/length(x))/sd(x)^3)
sk6 <- sapply(ga, function(x) sum((x-mean(x))^3/length(x))/sd(x)^3)


ku <- sapply(ga, function(x) sum((x-mean(x))^4/length(x))/sd(x)^4 - 3)

idat <- (reshape2::melt(apply(dplyr::select(alldat2, nComp, nMut, nPred, nComm, nAmens), 1, function(x) x/sum(x))))
ggplot(idat, aes(x = Var2, y = value)) + geom_bar(aes(fill = Var1), stat = "identity") + scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2"))

d1 <- c()
for(i in 50:300){
  d1[i-49] <- median(allfitRe$s[allfitRe$N == i]) - median(alldat2$fz[alldat2$N == i]) 
}

m1 <- c()
m2 <- c()
for(i in 1:300){
  m1[i] <- median(allfitRe$s[allfitRe$N == i])
  m2[i] <- median(alldat2$fz[alldat2$N == i])
}
plot(m1, ylim = c(0,5), typ = "l")
points(m2, typ = "l")
