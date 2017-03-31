library(vegan)
library(sads)

############################################################################################################
############################################################################################################
############################################################################################################

# Yingnan's modifief fitting function for the zipf dist RAD
fitzipf_r <- function(x, N, trunc, start.value, upper = 20, ...){
  if (any(x <= 0)) stop ("All x must be positive")
  if(class(x)!="rad") rad.tab <- rad(x)
  else rad.tab <- x
  # y <- rep(rad.tab$rank, rad.tab$abund)
  dots <- list(...)
  if (!missing(trunc)){
    if (min(rad.tab$rank)<=trunc) stop("truncation point should be lower than the lowest rank") #
  }
  if(missing(N)){
    N <- max(rad.tab$rank) #
  }
  if(missing(start.value)){
    p <- rad.tab$abund/sum(rad.tab$abund)
    lzipf <- function(s, N) -s*log(1:N) - log(sum(1/(1:N)^s))
    opt.f <- function(s) sum((log(p) - lzipf(s, length(p)))^2)
    opt <- optimize(opt.f, c(0.5, length(p)))
    sss <- opt$minimum
  }
  else{
    sss <- start.value
  }
  if(missing(trunc)){
    LL <- function(N, s) -sum(rad.tab$abund*dzipf(rad.tab$rank, N=N, s=s, log = TRUE)) #
  }
  else{
    LL <- function(N, s) -sum(rad.tab$abund*dtrunc("zipf", x = rad.tab$rank, coef = list(N = N, s = s), trunc = trunc, log = TRUE)) #
  }
  result <- do.call("mle2", c(list(LL, start = list(s = sss), data = list(x = rad.tab$rank), fixed=list(N=N), method = "Brent", lower = 0, upper = upper), dots))
  if(abs(as.numeric(result@coef) - upper) < 0.001)
    warning("mle equal to upper bound provided. \n Try increase value for the 'upper' argument")
  new("fitrad", result, rad="zipf", distr = "zipf of relative abundance", trunc = ifelse(missing(trunc), NaN, trunc), rad.tab=rad.tab)
}

# Function to get interaction numbers and strengths for a given community
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
  

  iL <- c("amensalism", "commensalism", "competition", "mutualism", "predation")
  num[iL %in% aggregate(df1$inty, list(df1$inty), length)$Group.1] <- aggregate(df1$inty, list(df1$inty), length)$x
  umu[iL %in% aggregate(df1$inty, list(df1$inty), length)$Group.1] <- aggregate(df1$i2, list(df1$inty), mean)$x
  lmu[iL %in% aggregate(df1$inty, list(df1$inty), length)$Group.1] <- aggregate(df1$i1, list(df1$inty), mean)$x
  
  
  df2 <- data.frame(typ = iL, num = num, m1 = umu, m2 = lmu)
  
  return(df2)
}


# My r2 function
get_r2 <- function(o, p){
  1 - sum((o-p)^2)/sum((o - mean(o))^2)
}

# Yingnan's modified r2 function
r2modified <- function(x,y,log=FALSE){
  if(log){
    rm = 1-sum((log10(x)-log10(y))^2)/sum((log10(x)-mean(log10(x)))^2)
  }
  else{
    rm = 1-sum((x-y)^2)/sum((x-mean(x))^2)
  }
  return(rm)  
}

############################################################################################################
############################################################################################################
############################################################################################################
### Fit zipf RAD to HMP dataset

otu2 <- read.csv("~/Desktop/otu_table_psn_v13.csv", row.names = 1)
metadat <- read.csv("~/Desktop/v13_map_uniquebyPSN.csv")

stoolsamp <- which(metadat$HMPbodysubsite == "Stool")
spptab <- colnames(otu2) %in% paste0("X",metadat[stoolsamp,]$SampleID)
otu3 <- otu2[-which(rowSums(otu2[,spptab]) == 0),spptab]


fzotu <- lapply(1:ncol(otu3), function(x) fitzipf_r(otu3[,x][otu3[,x]!=0]/sum(otu3[,x][otu3[,x]!=0])))
otuR2 <- sapply(1:ncol(otu3), function(x) r2modified(sort(otu3[,x][otu3[,x] != 0], decreasing = T), radpred(fzotu[[x]])$abund))

plot(t(sapply(fzotu, function(x) x@fullcoef)))
hist(otuR2)

## Check if sampling effort has an effect on s
### use resampling method to standardize # of reads

##################################################################
##################################################################
##################################################################
### Read in data

get_dat <- function(fpath, connected = TRUE){
  lf1 <- list.files(fpath)
  lf2 <- grep("ge", lf1)
  lf3 <- grep("mat", lf1)
  
  eqmat <- list()
  iconn <- c()
  mdstr <- c()
  for(i in 1:length(lf2)){
    ge1 <- readRDS(paste(fpath, lf1[lf2][[i]], sep = ""))
    if(any(is.na(ge1))){next}
    if(any(is.na(ge1$eqst))){next}
    
    mat1 <- readRDS(paste(filepath1, lf1[lf3][[i]], sep = ""))
    
    iconn[i] <- is.connected(graph.adjacency(abs(sign(mat1[ge1$spp, ge1$spp]))))
    
    eqmat[[i]] <- data.frame(itystr(mat1[ge1$spp, ge1$spp]), web = i, N = sum(ge1$spp))
    eqabs[[i]] <- sort(ge1$eqst, decreasing = T)
    
    mdstr[i] <- mean(diag(eqmat))
    
    cat(round(i/length(lf2)), "--:::--")
  }
  
  eqa <- eqabs[!is.na(eqabs) & !sapply(eqabs, is.null)]
  eqm <- eqmat[!is.na(eqabs) & !sapply(eqabs, is.null)]
  mdstr <- mdstr[!is.na(eqabs) & !sapply(eqabs, is.null)]
  
  if(connected){
    eqa <- eqa[iconn[!is.na(eqabs) & !sapply(eqabs, is.null)]]
    eqm <- eqm[iconn[!is.na(eqabs) & !sapply(eqabs, is.null)]]
    mdstr <- mdstr[iconn[!is.na(eqabs) & !sapply(eqabs, is.null)]]
  }
  
  return(list(eqa = eqa, eqm = eqm, ds = mdstr))
}


## File path location of the data
filepath1 <- "~/Documents/Data/parSAD_data2/"

fzmod <- function(x){
  fz1 <- fitzipf_r(x)
  fc1 <- fz1@fullcoef
  nll <- fz1@minuslogl(N = fz1@fullcoef[1], s = fz1@fullcoef[2])
  r2 <- r2modified(sort(x, decreasing = T), radpred(fz1)$abund)
  
  return(data.frame(t(fc1), nll, r2))
}

# Fit zipf RAD to eq abundances for random communities
## Eliminate NAs
eqabsA <- eqabs[!is.na(eqabs)]
## fitting
simfz <- lapply(eqabsA, function(x){x <- x[x > 0]; fitzipf_r(x)})
plot(t(sapply(simfz, function(x) x@fullcoef)))    
## Get R squared
simR2 <- sapply(1:length(eqabsA), function(x) r2modified(sort(eqabsA[[x]][eqabsA[[x]] > 0],decreasing = T), radpred(simfz[[x]])$abund))
sapply(simfz[which(simR2 > 0.9)], function(x) x@fullcoef)


# Fit zipf RAD to eq abundances for hub-like communities
## fitting to all non-null return communities
simfzhub <- lapply(eqabs2[!sapply(eqabs2, is.null)], function(x){x <- x[x > 0]; fitzipf_r(x)})
## Get R squared
simR2hub <- sapply(1:length(eqabs2[!sapply(eqabs2, is.null)]), function(x) r2modified(sort(eqabs2[!sapply(eqabs2, is.null)][[x]][eqabs2[!sapply(eqabs2, is.null)][[x]] > 0],decreasing = T), radpred(simfzhub[[x]])$abund))
hist(simR2hub)

# Fit zipf RAD to random communities with varying pars
sim2fz <- lapply(eqabs3[!sapply(eqabs3, is.null)], function(x){x <- x[x >0]; fitzipf_r(x)})
plot(t(sapply(sim2fz, function(x) x@fullcoef)))    
## Get R squared
sim2R2 <- sapply(1:length(eqabs3[!sapply(eqabs3, is.null)]), function(x){r2modified(sort(eqabs3[!sapply(eqabs3, is.null)][[x]][eqabs3[!sapply(eqabs3, is.null)][[x]] > 0],decreasing = T), radpred(sim2fz[[x]])$abund)})
hist(sim2R2)

# Compute connectance of equilibrium matrices
## Randoms
conn1 <- sapply(eqmat[!is.na(eqabs)], function(x) sum(x$num)/(x$N[1] *x$N[1]))
## Hubs
conn2 <- sapply(eqmat2[!is.na(eqabs2)], function(x) sum(x$num)/(x$N[1] *x$N[1]))
## Randoms with diff pars
conn3 <- sapply(eqmat3[!sapply(eqabs3, is.null)], function(x) sum(x$num)/(x$N[1] *x$N[1]))
## Plotting fitted pars against connectance 
plot(conn1[!is.na(eqabs)], sapply(simfz, function(x) x@coef))
plot(conn2, sapply(simfzhub, function(x) x@coef)[!is.na(eqabs2)])


# Get vector of fitted par (s) for random dat
s.val <- sapply(simfz, function(x) x@coef)
s.val2 <- sapply(simfzhub, function(x) x@coef)
s.val3 <- sapply(sim2fz, function(x) x@coef)
# Get matrix of all fitted pars for random, hub, and real comms
fitpars <- rbind(cbind(t(sapply(simfz, function(x) x@fullcoef)),typ = 1, r2 = simR2),
                 cbind(t(sapply(simfzhub, function(x) x@fullcoef)),typ = 2, r2 = simR2hub),
                 cbind(t(sapply(fzotu, function(x) x@fullcoef)), typ = 3, r2 = otuR2),
                 cbind(t(sapply(sim2fz, function(x) x@fullcoef)), typ = 4, r2 = sim2R2))
fitpars[,"r2"][fitpars[,"r2"] < 0]  <- 0   ## make any negative rsquared 0

# plot relationship between comm size and s
plot(s~N, col = typ, data = fitpars, pch = 20)
ggplot(data.frame(fitpars), aes(x = N, y = s, col = factor(typ), alpha = r2)) + geom_point() + theme_bw()


##################################################################
##################################################################
##################################################################
### Interactions
#testing objs
# add self interaction mean

eqa <- eqabs3[!sapply(eqabs3, is.null)]
eqm <- eqmat3[!sapply(eqabs3, is.null)]
prepdat <- function(eqa, eqm, svals, sr2){
  Nspp <- sapply(eqa, length)
  conn <- sapply(eqm, function(x) sum(x$num)/(x$N[1] *x$N[1]))
  
  # pull out interaction numbers
  int1 <- t(sapply(eqm, function(x){if(is.null(x)){return(c(NA,NA,NA,NA,NA))};x$num}))
  istrs <- lapply(eqm, function(x){x$m1[2] <- x$m2[2]; x$m2[2] <- 0;return(x)}) ## move main commensal strength with others
  # pull out interaction strengths
  int3 <- t(sapply(istrs, function(x){if(any(is.na(unlist(x)))){return(c(NA,NA,NA,NA,NA))};x$m1}))
  # get positive predation strength 
  p2 <- sapply(istrs, function(x) x$m2[5])
  # put interaction strengths into one matrix
  allint <- cbind(int1, int3, p2)
  colnames(allint) <- c("aN", "coN", "cpN", "mN", "pN", "aS", "coS", "cpS", "mS", "pSn", "pSp")
  # put ints and fitted par into one dataframe
  dat <- data.frame(sV = svals, sR = sr2, allint, Nsp = Nspp, C = conn, D = dstr)
  
  return(dat)
}


pdat <- prepdat(eqa = eqa, eqm = eqm, svals = s.val3, sr2 = sim2R2)
subdat <- pdat[iconn2[!sapply(eqabs3, is.null)],]
fit.init <- (lm(sV~Nsp+C+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = subdat))
fit.init.ints <- (lm(sV~Nsp+C+aN*aS+coN*coS+cpN*cpS+mN*mS+pN*pSn+pSp+pN:pSp, data = pdat[iconn2[!sapply(eqabs3, is.null)],], na.action = "na.fail"))
summary(fit.init)

DAAG::cv.lm(data = subdat, form.lm = fit.init, m = 3)

mse.init <- sum((fitted(fit.init) - pdat$sV[iconn2[!sapply(eqabs3, is.null)]])^2)

head(MuMIn::dredge(fit.init.ints))
predict(fit.init, newdata = pdat[iconn2[!sapply(eqabs3, is.null)],][pdat[iconn2[!sapply(eqabs3, is.null)],],])
pdat$sV[iconn2[!sapply(eqabs3, is.null)]][100]

#fitlinmod <- function(pdat, r2cutoff = 0.5, boot = F){
#  fit.init <- (lm(sV~Nsp+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = pdat))
#  summary(fit.init)
#}

Nspp <- sapply(eqabs[!is.na(eqabs)], length)
# pull out interaction numbers
int1 <- t(sapply(eqmat[!is.na(eqabs)], function(x){if(is.null(x)){return(c(NA,NA,NA,NA,NA))};x$num}))
istrs <- lapply(eqmat, function(x){x$m1[2] <- x$m2[2]; x$m2[2] <- 0;x$m1[4] <- x$m2[4]; x$m2[4] <- 0;return(x)}) ## move main commensal strength with others
# pull out interaction strengths
int3 <- t(sapply(istrs[!is.na(eqabs)], function(x){if(any(is.na(unlist(x)))){return(c(NA,NA,NA,NA,NA))};x$m1}))
# get positive predation strength 
p2 <- sapply(istrs[!is.na(eqabs)], function(x) x$m2[5])
# put interaction strengths into one matrix
allint <- cbind(int1, int3, p2)
colnames(allint) <- c("aN", "coN", "cpN", "mN", "pN", "aS", "coS", "cpS", "mS", "pSn", "pSp")
# put ints and fitted par into one dataframe
dat <- data.frame(sR = simR2, allint, Nsp = Nspp)
# rescale by column means (units of 2 standard deviations)
dat2 <- apply(dat, 2, function(x) (x - mean(x))/(2*sd(x)))
dat2[,1] <- dat$sR
#dat2[,1][dat2[,1] < 0] <- 0 

# fit linear models
## to scaled data
summary(lm(sR~Nsp+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = data.frame(dat2[dat2[,1] > 0,])))
## to original data, but only for those where zipf has good fit
fit.init <- (lm(sR~Nsp+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = dat[dat$sR > .5,]))
summary(fit.init)

## subset original data
mydat <- dat[dat$sR > .5,]

## bootstrapping confidence intervals for linmod coefficients
coefs.bs <- matrix(nrow = 200, ncol = 13)
colnames(coefs.bs) <- names(coefficients(fit.init))
r2val.bs <- c()
for(i in 1:200){
  bs.rows <- sample(1:nrow(mydat), nrow(mydat), replace = T)
  fit.bs <- lm(sR~Nsp+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = mydat[bs.rows,])
  
  coefs.bs[i,] <- coefficients(fit.bs)
  r2val.bs[i] <- summary(fit.bs)$r.squared
}
# get intervals
apply(coefs.bs, 2, function(x) quantile(x, probs = c(0.025, 0.975)))

#MuMIn::dredge(lm(sR~aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = dat[dat$sR > .5,], na.action = "na.fail"))

## Check model for how interactions affect fitting of the zipf (rsquared)
summary(glm(as.numeric(simR2 < 0.5)~allint, family = "binomial"))
DAAG::CVbinary(glm(sR~Nsp+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = dat, family = "binomial"))




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
