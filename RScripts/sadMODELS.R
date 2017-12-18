get_fz <- function(x, reads = 2000){
  ab1 <- apply(x, 1, function(x) unlist(x)[-c(601,602)])
  ab1 <- apply(ab1, 2, get_abundvec, N = reads)
  
  ab1 <- lapply(ab1, fzmod)
  
  return(do.call(rbind, ab1))
}

get_fz1 <- function(x, reads = 2000){
  ab1 <- apply(x, 1, function(x) unlist(x)[-c(601,602)])
  
  ab2 <- list()
  for(x in 1:ncol(ab1)){
    fz <- list()
    for(i in 1:100){
      ga <- get_abundvec(ab1[,x], N = reads)
      fz[[i]] <- fzmod(ga)
    }
    print(x)
    ab2[[x]] <- colMeans(do.call(rbind, fz))
  }
  
  
  return(do.call(rbind, ab2))
}


#ad1 <- readRDS("~/Documents/ImData2/alldat2.rds")
#ad2 <- readRDS("~/Documents/ImData2/alldat3.rds")

ad1 <- readRDS("D:/UVA/ImData2/alldat2.rds")
ad2 <- readRDS("D:/UVA/ImData2/alldat3.rds")


wp1 <- ad1[[2]]
wp2 <- ad2[[2]]
fdat1 <- ad1[[1]]
fdat2 <- ad2[[1]]
ad1 <- ad1[[3]]
ad2 <- ad2[[3]]

#rm(ad1)
#rm(ad2)


ab1a <- get_fz(ad1, reads = 10000)
ab2a <- get_fz(ad2, reads = 10000)

saveRDS(ab1a, "D:/microbiome-dynamics/data/fz_fit_1a.rds")
saveRDS(ab2a, "D:/microbiome-dynamics/data/fz_fit_2a.rds")

saveRDS(rbind(ab1a, ab2a), "D:/microbiome-dynamics/data/fz_fit_alla.rds")

ab1b <- get_fz1(ad1, reads = 10000)
saveRDS(ab1b, "D:/microbiome-dynamics/data/fz_fit_1b.rds")
ab2b <- get_fz1(ad2, reads = 10000)
saveRDS(ab2b, "D:/microbiome-dynamics/data/fz_fit_2b.rds")

plot(rbind(ab1b,ab2b)[,1:2])


otu <- read.csv("F:/AmericanGut/amergutOTU.csv")
fzag <- apply(otu[,-1], 2, get_abundvec, N = 10000)
fzag <- lapply(fzag, fzmod)
fzag <- do.call(rbind, fzag)
saveRDS(fzag, "D:/microbiome-dynamics/data/fzag10k.rds")

plot(fzag[,1:2], xlim = c(0, 600), ylim = c(.5, 3.5))
points(ab1[,1:2], col = "green3")
points(ab2, col = "green4")

#ab1 <- rbind(ab1, ab2)
#wp1 <- rbind(wp1, wp2)
#fdat1 <- rbind(fdat1, fdat2)


fit1 <- lm(log10(ab1$s)~cc.w+conn+diam.uw+apl+mod.uw+m.int+c1+ls, data = wp1, x = F, y = F, model = F)
summary(fit1)
fit2 <- lm(log10(ab1$s)~nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut, data = fdat1, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2)
fit2.2 <- lm(log10(ab1$s)~m.int+comp+mut+pred+amens+comm, data = wp1, x = F, y = F, model = F, na.action = "na.fail")
summary(fit2.2)
fit3 <- lm(log10(ab1$s)~bet.uw+d.tot+cc.uw+mod.uw, data = fdat1, x = F, y = F, model = F, na.action = "na.fail")
summary(fit3)
fit3.2 <- lm(log10(ab1$s)~bet.w+d.tot+cc.w+mod.uw+conn+diam.uw+diam.w+apl+mod.uw+m.int+c1+ls, data = data.frame(wp1,fdat1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit3.2)
fit4 <- lm(log10(ab1$s)~bet.w+d.tot+cc.w+mod.uw+nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut+conn+diam.uw+diam.w+apl+mod.uw+m.int+c1+ls, data = data.frame(wp1,fdat1), x = F, y = F, model = F, na.action = "na.fail")
summary(fit4)


library(quantreg)
frq <- rq(s~N, tau = .025,  data = fzag)
abline(frq)
cutoff <- matrix(c(1:450, predict(frq, data.frame(N = 1:450))), ncol = 2)
ir <- c()
for(i in 1:nrow(ab1)){
   ir[i] <- ab1$s[i] > cutoff[cutoff[,1] %in% ab1$N[i], 2]
}
sum(!ir)

frqL <- rq((s)~(N), tau = .975,  data = log10(allf))
abline(frqL)
cutoff2 <- matrix(c(10:600, predict(frqL, data.frame(N = log10(10:600)))), ncol = 2)
irL <- c()
for(i in 1:nrow(do.call(rbind, im0k0)[,1:2])){
  irL[i] <- (do.call(rbind, im0k0)[,2][i]) > 10^(cutoff[cutoff[,1] %in% (do.call(rbind, im0k0)[,1][i]), 2]) & (do.call(rbind, im0k0)[,2][i]) < 10^(cutoff2[cutoff2[,1] %in% (do.call(rbind, im0k0)[,1][i]), 2])
}
sum(!irL)

fit5 <- glm(ir~bet.w+d.tot+cc.w+mod.uw+nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut+conn+diam.uw+diam.w+apl+mod.uw+m.int+c1+ls, data = data.frame(ir, wp1, fdat1), family = "binomial")
summary(fit5)
DAAG::cv.binary(fit5)




plot(fzag[,1:2], xlim = c(0, 600), ylim = c(.5, 3.5))
points(ab1[,1:2], col = ifelse(ir, "green4", "blue2"))
lines(cutoff[,1], cutoff[,2], col = "blue")



confusion.glm <- function(model, des.mat=NULL, response=NULL, cutoff=0.5) {
  if (missing(des.mat)) {
    prediction <- predict(model, type='response') > cutoff
    confusion  <- table(as.logical(model$y), prediction)
  } else {
    if (missing(response) || class(response) != "logical") {
      stop("Must give logical vector as response when des.mat given")
    }
    prediction <- predict(model, des.mat, type='response') > cutoff
    confusion  <- table(response, prediction)
  }
  confusion <- cbind(confusion,
                     c(1 - confusion[1,1] / rowSums(confusion)[1],
                       1 - confusion[2,2] / rowSums(confusion)[2]))
  confusion <- as.data.frame(confusion)
  names(confusion) <- c('FALSE', 'TRUE', 'class.error')
  return(confusion)
}

con1 <- confusion.glm(fit5)

ssmod <- function(conf){
  sens <- conf[2,2]/sum(conf[2,1:2])
  spec <- conf[1,1]/sum(conf[1,1:2])
  err.all <- sum(conf[1,1], conf[2,2])/sum(conf[1:2,1:2])
  
  return(c(sensitivity = sens, specificity = spec, overall_error = 1-err.all))
}

ssmod(confusion.glm(fit5))
ssmod(confusion.glm(fit5.2))

rfir <- randomForest::randomForest(factor(irL)~bet.w+d.tot+cc.w+mod.uw+nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut+conn+diam.uw+diam.w+apl+mod.uw+m.int+c1+ls, data = data.frame(irL, wp1, fdat1), mtry = 20, ntree = 2000, importance = T)
randomForest::varImpPlot(rfir)

cart  <- rpart::rpart(factor(irL)~bet.w+d.tot+cc.w+mod.uw+nComp+CompOut+nMut+MutOut+nPred+PredOut+nAmens+AmensOut+nComm+CommOut+conn+diam.uw+diam.w+apl+mod.uw+m.int+c1+ls, data = data.frame(irL, wp1, fdat1), method = "class")
#rpart::plotcp(cart)
rpart.plot::prp(cart, extra = 1)
ssmod(table((predict(cart, type = "class")), irL))
ssmod(con1)



#######################################################################
