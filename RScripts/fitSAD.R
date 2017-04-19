############################################################################################################
############################################################################################################
############################################################################################################

library(vegan)
library(sads)
library(MASS)
library(rpart)
library(ggplot2)

############################################################################################################
############################################################################################################
############################################################################################################
# Yingnan's modified fitting function for the zipf dist RAD
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


get_abundvec <- function(abund, N = 10000){
  r.ab <- abund/sum(abund)
  samp2 <- sample(1:length(abund), N, replace = T, prob = r.ab)
  t1 <- table(samp2)
  return(as.numeric(t1))
}


fzmod <- function(x){
  fz1 <- fitzipf_r(x)
  fc1 <- fz1@fullcoef
  nll <- fz1@minuslogl(N = fz1@fullcoef[1], s = fz1@fullcoef[2])
  r2 <- r2modified(sort(x, decreasing = T), radpred(fz1)$abund)
  
  return(data.frame(t(fc1), nll, r2))
}


fzmod2 <- function(x, rad = "zipf"){
  rad1 <- paste("fit", rad, "_r", sep = "")
  rfit <- get(rad1)
  
  fz1 <- rfit(x)
  fc1 <- fz1@fullcoef
  nll <- fz1@minuslogl(fz1@fullcoef[1], fz1@fullcoef[2])
  r2 <- r2modified(sort(x, decreasing = T), radpred(fz1)$abund)
  
  return(data.frame(t(fc1), nll, r2))
}


prepdat <- function(eqa, eqm, svals, sr2, d, r){
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
  dat <- data.frame(sV = unlist(svals), sR = unlist(sr2), allint, Nsp = Nspp, C = conn, D = sapply(d, mean), r = sapply(r, mean))
  
  return(dat)
}


get_dat <- function(fpath, connected = TRUE){
  lf1 <- list.files(fpath)
  lf2 <- grep("ge", lf1)
  lf3 <- grep("mat", lf1)
  
  eqmat <- list()
  eqabs <- list()
  rs <- list()
  iconn <- c()
  mdstr <- list()
  wrks <- c()
  for(i in 1:length(lf2)){
    ge1 <- readRDS(paste(fpath, lf1[lf2][[i]], sep = ""))
    if(any(is.na(ge1))){next}
    if(any(is.na(ge1$eqst))){next}
    
    mat1 <- readRDS(paste(fpath, lf1[lf3][[i]], sep = ""))
    
    iconn[i] <- is.connected(graph.adjacency(abs(sign(mat1[ge1$spp, ge1$spp]))))
    
    eqmat[[i]] <- data.frame(itystr(mat1[ge1$spp, ge1$spp]), web = i, N = sum(ge1$spp))
    eqabs[[i]] <- sort(ge1$eqst, decreasing = T)
    
    mdstr[[i]] <- (diag(mat1[ge1$spp, ge1$spp]))#mean(diag(mat1[ge1$spp, ge1$spp]))
    wrks[i] <- i
    rs[[i]] <- ge1$grs
    
    if(i%%100 == 0){cat(round(i/length(lf2)*100), "--:::--")}
  }
  
  wrks <- wrks[!is.na(eqabs) & !sapply(eqabs, is.null)]
  eqa <- eqabs[!is.na(eqabs) & !sapply(eqabs, is.null)]
  eqm <- eqmat[!is.na(eqabs) & !sapply(eqabs, is.null)]
  mdstr <- mdstr[!is.na(eqabs) & !sapply(eqabs, is.null)]
  rs <- rs[!is.na(eqabs) & !sapply(eqabs, is.null)]
  
  if(connected){
    eqa <- eqa[iconn[!is.na(eqabs) & !sapply(eqabs, is.null)]]
    eqm <- eqm[iconn[!is.na(eqabs) & !sapply(eqabs, is.null)]]
    mdstr <- mdstr[iconn[!is.na(eqabs) & !sapply(eqabs, is.null)]]
    wrks <- wrks[iconn[!is.na(eqabs) & !sapply(eqabs, is.null)]]
    rs <- rs[iconn[!is.na(eqabs) & !sapply(eqabs, is.null)]]
  }
  
  return(list(eqa = eqa, eqm = eqm, ds = mdstr, rs = rs, wrkd = wrks))
}


get_mat <- function(fpath, nums){
  eqmat <- list()
  lf1 <- list.files(fpath)
  lf2 <- grep("ge", lf1)
  lf3 <- grep("mat", lf1)
  for(i in 1:length(nums)){
    ge1 <- readRDS(paste(fpath, lf1[lf2][[nums[i]]], sep = ""))
    if(any(is.na(ge1))){next}
    if(any(is.na(ge1$eqst))){next}
    
    mat1 <- readRDS(paste(fpath, lf1[lf3][[nums[i]]], sep = ""))
    
    eqmat[[i]] <- mat1[ge1$spp, ge1$spp]
    
    if(i%%100 == 0){cat(ceiling(i/length(nums)*100), "--:::--")}
  }
  
  return(eqmat)
}

get_eqa <- function(fpath, nums){
  eqmat <- list()
  lf1 <- list.files(fpath)
  lf2 <- grep("ge", lf1)
  lf3 <- grep("mat", lf1)
  for(i in 1:length(nums)){
    ge1 <- readRDS(paste(fpath, lf1[lf2][[nums[i]]], sep = ""))
    if(any(is.na(ge1))){next}
    if(any(is.na(ge1$eqst))){next}
    
    eqmat[[i]] <- ge1$eqst
    
    if(i%%100 == 0){cat(ceiling(i/length(nums)*100), "--:::--")}
  }
  
  return(eqmat)
}
############################################################################################################
############################################################################################################
############################################################################################################
### Load HMP dataset

otu2 <- read.csv("~/Desktop/otu_table_psn_v13.csv", row.names = 1)
metadat <- read.csv("~/Desktop/v13_map_uniquebyPSN.csv")

stoolsamp <- which(metadat$HMPbodysubsite == "Stool")
spptab <- colnames(otu2) %in% paste0("X",metadat[stoolsamp,]$SampleID)
otu3 <- otu2[-which(rowSums(otu2[,spptab]) == 0),spptab]
rm(otu2)
rm(metadat)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

otuT <- read.csv("Desktop/Archive/feces_M3_spp.csv")
fzT <- apply(otuT[,-1], 1, function(x) fzmod(sort(x[x!=0])))
fzT <- do.call(rbind, fzT)

otuT2 <- read.csv("Desktop/Archive/feces_F4_spp.csv")
fzT2 <- apply(otuT2[,-1], 1, function(x) fzmod(sort(x[x>0])))
fzT2 <- do.call(rbind, fzT2)

dtgut1 <- read.csv("Desktop/Archive/dtgut1.csv")
dt1 <- apply(dtgut1[complete.cases(dtgut1),-1], 1, function(x) fzmod(sort(x[x>0])))
dt1 <- do.call(rbind, dt1)

dtgut2 <- read.csv("Desktop/Archive/dtgut2.csv")
dt2 <- apply(dtgut2[complete.cases(dtgut2),-1], 1, function(x) fzmod(sort(x[x>0])))
dt2 <- do.call(rbind, dt2)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
### Read in data

# Fit zipf RAD to eq abundances for random communities
filepath1 <- "~/Documents/Data/parSAD_data/"
st1 <- Sys.time()
psd1 <- get_dat(filepath1)
st2 <- Sys.time()
st2-st1

#fzd1 <- t(sapply(psd2$eqa, fzmod))

# Fit zipf RAD to random communities with varying pars
filepath2 <- "~/Documents/Data/parSAD_data2/"
st3 <- Sys.time()
psd2 <- get_dat(filepath2)
st4 <- Sys.time()
st4-st3

#fzd2 <- t(sapply(psd2$eqa, fzmod))
#fzd2b <- t(sapply(psd2$eqa, function(x) fzmod(get_abundvec(x, 100))))

# Fit zipf RAD to eq abundances for hub-like communities

filepath3 <- "~/Documents/Data/parSADhub_data/"
st5 <- Sys.time()
psd3 <- get_dat(filepath3)
st6 <- Sys.time()
st6-st5

#fzd3 <- t(sapply(psd3$eqa, fzmod))

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
## Standardize reads to X
X <- 2000

gav1 <- apply(otu3, 2, get_abundvec, N = X)
gav2 <- apply(otu3, 2, function(x) get_abundvec(x[x>=5], N = X))
gavfz1 <- t(sapply(gav1, fzmod))
gavfz2 <- t(sapply(gav2, fzmod))
s.hmp <- unlist(gavfz1[,"s"])[which(apply(otu3, 2, sum) > X)]
s.hmp2 <- unlist(gavfz2[,"s"])[which(apply(otu3, 2, sum) > X)]
n.hmp <- unlist(gavfz1[,"N"])[which(apply(otu3, 2, sum) > X)]
n.hmp2 <- unlist(gavfz2[,"N"])[which(apply(otu3, 2, sum) > X)]
r2.hmp <- unlist(gavfz1[,"r2"])[which(apply(otu3, 2, sum) > X)]

gavt <-  apply(otuT[,-1], 1, get_abundvec, N = X)
fzT <- lapply(gavt, function(x) fzmod(sort(x)))
fzT <- do.call(rbind, fzT)

gavt2 <-  apply(otuT2[,-1], 1, get_abundvec, N = X)
fzT2 <- lapply(gavt, function(x) fzmod(sort(x)))
fzT2 <- do.call(rbind, fzT2)

gavdt <-  apply(dtgut1[complete.cases(dtgut1),-1], 1, get_abundvec, N = X)
dt1 <- lapply(gavdt, function(x) fzmod(sort(x)))
dt1 <- do.call(rbind, dt1)

gavdt2 <-  apply(dtgut2[complete.cases(dtgut2),-1], 1, get_abundvec, N = X)
dt2 <- lapply(gavdt2, function(x) fzmod(sort(x)))
dt2 <- do.call(rbind, dt2)
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

allfit <- data.frame(s = c(fzT$s, fzT2$s, dt1$s, dt2$s, s.hmp),
                     N = c(fzT$N, fzT2$N, dt1$N, dt2$N, n.hmp),
                     r2 = c(fzT$r2, fzT2$r2, dt1$r2, dt2$r2, r2.hmp),
                     dat = rep(c("M3", "F4", "DT1", "DT2", "HMP"), c(length(fzT$s),length(fzT2$s), length(dt1$s), length(dt2$s), length(s.hmp))))


ggplot(allfit[allfit$r2 > 0.9,], aes(x = N, y = s, col = dat)) + geom_point() + geom_smooth() + theme_bw()
ggplot(allfit, aes(x = r2, y = s, col = dat)) + geom_point() + geom_smooth() + theme_bw()

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
gavsim <- lapply(psd2$eqa, get_abundvec, X)
simfz1 <- lapply(gavsim, function(x) fzmod(sort(x)))
simfz1 <- do.call(rbind, simfz1)

gavsim2 <- lapply(psd1$eqa, function(x) get_abundvec(x[x>0], N = X))
simfz2 <- lapply(gavsim2, function(x) fzmod(sort(x)))
simfz2 <- do.call(rbind, simfz2)

gavsim3 <- lapply(psd3$eqa, function(x) get_abundvec(x[x>0], N = X))
simfz3 <- lapply(gavsim3, function(x) fzmod(sort(x)))
simfz3 <- do.call(rbind, simfz3)

allfit <- data.frame(s = c(fzT$s, fzT2$s, dt1$s, dt2$s, s.hmp, simfz1$s, simfz2$s, simfz3$s),
                     N = c(fzT$N, fzT2$N, dt1$N, dt2$N, n.hmp, simfz1$N, simfz2$N, simfz3$N),
                     r2 = c(fzT$r2, fzT2$r2, dt1$r2, dt2$r2, r2.hmp, simfz1$r2, simfz2$r2, simfz3$r2),
                     dat = rep(c("M3", "F4", "DT1", "DT2", "HMP", "sim", "sim2", "sim3"),
                               c(length(fzT$s),length(fzT2$s),length(dt1$s),length(dt2$s),length(s.hmp),nrow(simfz1), nrow(simfz2), nrow(simfz3))))


ggplot(allfit, aes(x = N, y = s, col = dat)) + geom_point() + geom_smooth() + theme_bw() 
ggplot(allfit, aes(x = r2, y = s, col = dat)) + geom_point() + geom_smooth() + theme_bw() 

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# make sure allfit does not include sim dat
rq <- matrix(nrow = nrow(simfz1), ncol = 2)
inrange <- c()
for(i in 1:nrow(simfz1)){
  if(simfz1[i,]$N > 450){
    rq[i,] <- range(allfit$s[allfit$dat != "sim" & allfit$N %in% 400:600])
  }else{
    rq[i,] <- range(allfit$s[allfit$dat != "sim" & allfit$N %in% (simfz1[i,]$N-5):(simfz1[i,]$N+5)])
  }
  inrange[i] <- simfz1[i,]$s <= rq[i,2] & simfz1[i,]$s >= rq[i,1]
}
sum(inrange)
length(inrange)

inrange2 <- c()
test <- c()
for(i in 1:length(t2k$sV)){
  if(t2k$Nsp[i] > 450){
    rq[i,] <- range(allfit$s[allfit$dat != "sim" & allfit$N %in% 1:100])
  }else{
    test[i] <- sum(allfit$N %in% (t2k$Nsp[i]-5):(t2k$Nsp[i]+5))
    rq[i,] <- range(allfit$s[allfit$dat != "sim" & allfit$N %in% (t2k$Nsp[i]-20):(t2k$Nsp[i]+20)])
  }
  inrange2[i] <- t2k$sV[i] <= rq[i,2] & t2k$sV[i] >= rq[i,1]
}


len1 <- sapply(gavsim,length)[(sapply(gavsim, length) > 50)]
hmp1 <- sapply(gav1[apply(otu3, 2, sum) > 2000], function(x) rev(sort(x)))
ir <- c()
q1 <- c()
for(j in 281:length(len1)){
  len2 <- len1[j]
  fzmat <- matrix(nrow = length(hmp1[sapply(hmp1, length) > len2]), ncol = 2)
  if(nrow(fzmat) < 2){next}
  for(i in 1:nrow(fzmat)){
    fzmat[i,] <- fitzipf_r(head(hmp1[sapply(hmp1, length) > len2][[i]], len2))@fullcoef
  }
  r1 <- range(fzmat[,2])
  s.j <- fitzipf_r(gavsim[[j]])@coef 
  
  q1[j] <- sum(fzmat[,2] < s.j)/nrow(fzmat)
  ir[j] <- s.j <= r1[2] & s.j >= r1[1]
  print(j)
} 



####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
### Interactions
#testing objs
# add self interaction mean

#pdat <- prepdat(eqa = psd2$eqa, eqm = psd2$eqm, svals = fzd2[,"s"], sr2 = fzd2[,"r2"], d = psd2$ds)
pdat2 <- prepdat(eqa = gavsim, eqm = psd2$eqm, svals = simfz1$s, sr2 = simfz1$r2, d = psd2$ds, r = psd2$rs)
pdat2$abs <- sapply(psd2$eqa, sum)

subdat <- apply(pdat2[,-c(1,2)], 2, function(x){(x - mean(x))/sd(x)})
subdat <- data.frame(pdat2[,c(1,2)], subdat)
subdat2 <- subdat[pdat2$Nsp > 50,]
subdat2$ir <- ir
subdat2$Nsp <- subdat$Nsp[pdat2$Nsp > 50]
subdat2$Nsp2 <- pdat2$Nsp[pdat2$Nsp > 50]


#######################################
# Models ##############################
#######################################

# all orig data modeling s value
fit1 <- (lm(sV~Nsp+C+r+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = pdat2))
summary(fit1)
# orig data (good fit) modeling s value
fit2 <- (lm(sV~Nsp+C+r+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = pdat2[pdat2$sR > 0.8,]))
summary(fit2)
# all orig data modeling r2
fit3 <- (lm(sR~Nsp+C+r+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = pdat2))
summary(fit3)

#######################################
#######################################
#######################################

# all data as zscore modeling s value
fit1.1 <- (lm(sV~Nsp+C+r+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = subdat))
summary(fit1.1)
# zscore data (good fit) modeling s value
fit2.1 <- (lm(sV~Nsp+C+r+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = subdat[subdat$sR > 0.8,]))
summary(fit2.1)
# zscore data modeling r2
fit3.1 <- (lm(sR~Nsp+C+r+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = subdat))
summary(fit3.1)

fit4 <- lm(cbind(sV, Nsp, abs)~r+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = subdat)
summary(fit4)

subdat3 <- data.frame(mu = t(sa2)[,1], subdat)
fit5 <- (lm(mu~Nsp+C+r+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = subdat3))
summary(fit5)
data.frame(mod1 = summary(fit1.1)$coefficients[,1], mod2 = summary(fit3.1)$coefficients[,1], p1 = summary(fit1.1)$coefficients[,4] <= 0.05, p2 = summary(fit3.1)$coefficients[,4] <= 0.05)

#######################################
#######################################
#######################################


tfit <- (glm(ir~C+r+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = subdat2[complete.cases(subdat2),], family = "binomial"))
summary(tfit)
tfit2 <- glm(ir~C+r+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = subdat2[complete.cases(subdat2),][subdat2[complete.cases(subdat2),]$sR > 0.9,], family = "binomial")
summary(tfit2)
DAAG::cv.binary(tfit2)

subdat2$ir2 <- (subdat2$sR > 0.75)*subdat2$ir
tfit3 <- glm(ir2~Nsp, data = subdat2[complete.cases(subdat2),], family = "binomial")
summary(tfit3)
DAAG::cv.binary(tfit3)


##################################################
##################################################
pdat1 <- prepdat(eqa = gavsim2, eqm = psd1$eqm, svals = simfz2$s, sr2 = simfz2$r2, d = psd1$ds, r = psd1$rs)
pdat3 <- prepdat(eqa = gavsim3, eqm = psd3$eqm, svals = simfz3$s, sr2 = simfz3$r2, d = psd3$ds, r = psd3$rs)

sdat1 <- apply(pdat1[,-c(1,2)], 2, function(x){(x - mean(x))/sd(x)})
sdat1 <- data.frame(pdat1[,c(1,2)], sdat1)
sdat3 <- apply(pdat3[,-c(1,2)], 2, function(x){(x - mean(x))/sd(x)})
sdat3 <- data.frame(pdat3[,c(1,2)], sdat3)


fit1.p1 <- (lm(sV~Nsp+C+r+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = pdat1))
fit1.p3 <- (lm(sV~Nsp+C+r+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = pdat3))
summary(fit1.p1)
summary(fit1.p3)

lm(sV~Nsp+C+r+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = rbind(pdat1, pdat2, pdat3))
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
### Any effect of copy number variation

new.eqa <- lapply(psd2$eqa, function(x) x*sample(c(1,2), length(x), replace = T, prob = c(.95,.5)))
new.eqa <- lapply(psd2$eqa, function(x) rev(sort(x))*rev(sort(rpois(length(x), .01)+1)))
new.eqa <- lapply(psd2$eqa, function(x){x[which.max(x)] <- max(x)*2;return(x)})
new.eqa <- lapply(psd2$eqa, function(x){x[order(x) %in% 1:5] <- x[order(x) %in% 1:5]*sample(c(1,2), 5, replace = T, prob = c(.8,.2));return(x)})
gs.cnv <- lapply(new.eqa, get_abundvec, X)
sim.cnv <- lapply(gs.cnv, function(x) fzmod(sort(x)))
sim.cnv <- do.call(rbind, sim.cnv)

plot(sim.cnv$s, simfz1$s)
abline(a = 0, b = 1, xpd = F)


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
### A look at stability
allmats <- get_mat(filepath2, psd2$wrkd)
alleqa <- get_eqa(filepath2, psd2$wrkd)

eig <- c()
for(i in 1:length(psd2$eqa)){
  p <- list(alpha = psd2$rs[[i]], m = allmats[[i]], K = 20)
  jf <- jacobian.full(alleqa[[i]], lvmodK, parms = p)
  eig[i] <- max(Re(eigen(jf)$values))
  print(i)
}

plot(eig, simfz1$N)

fitE <- lm(eig~Nsp+C+r+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = subdat)
summary(fitE)
