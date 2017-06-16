############################################################################################################
############################################################################################################
############################################################################################################

library(vegan)
library(sads)
library(MASS)
library(rpart)
library(ggplot2)
library(MuMIn)
library(igraph)

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
  kvals <- c()
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
    kvals[i] <- ge1$eqkv 
    
    if(i%%100 == 0){cat(round(i/length(lf2)*100), "--:::--")}
  }
  
  wrks <- wrks[!is.na(eqabs) & !sapply(eqabs, is.null)]
  eqa <- eqabs[!is.na(eqabs) & !sapply(eqabs, is.null)]
  eqm <- eqmat[!is.na(eqabs) & !sapply(eqabs, is.null)]
  mdstr <- mdstr[!is.na(eqabs) & !sapply(eqabs, is.null)]
  rs <- rs[!is.na(eqabs) & !sapply(eqabs, is.null)]
  kvals <- kvals[!is.na(eqabs) & !sapply(eqabs, is.null)]
  
  if(connected){
    eqa <- eqa[iconn[!is.na(eqabs) & !sapply(eqabs, is.null)]]
    eqm <- eqm[iconn[!is.na(eqabs) & !sapply(eqabs, is.null)]]
    mdstr <- mdstr[iconn[!is.na(eqabs) & !sapply(eqabs, is.null)]]
    wrks <- wrks[iconn[!is.na(eqabs) & !sapply(eqabs, is.null)]]
    rs <- rs[iconn[!is.na(eqabs) & !sapply(eqabs, is.null)]]
    kvals <- kvals[iconn[!is.na(eqabs) & !sapply(eqabs, is.null)]]
  }
  
  return(list(eqa = eqa, eqm = eqm, ds = mdstr, rs = rs, wrkd = wrks, kv = kvals))
}


get_mat <- function(fpath, nums){
  eqmat <- list()
  lf1 <- list.files(fpath)
  lf2 <- grep("ge", lf1)
  lf3 <- grep("mat", lf1)
  for(i in 5551:7000){
    ge1 <- readRDS(paste(fpath, lf1[lf2][[i]], sep = ""))
    if(any(is.na(ge1))){next}
    if(any(is.na(ge1$eqst))){next}
    
    mat1 <- readRDS(paste(fpath, lf1[lf3][[i]], sep = ""))
    
    eqmat[[i]] <- mat1[ge1$spp, ge1$spp]
    
    print(i)
    #if(i%%100 == 0){cat(ceiling(i/length(nums)*100), "--:::--")}
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
fzT <- apply(otuT[,-1], 1, function(x) fzmod(get_abundvec(rev(sort(x[x!=0])), 2000)))
fzT <- do.call(rbind, fzT)

otuT2 <- read.csv("Desktop/Archive/feces_F4_spp.csv")
fzT2 <- apply(otuT2[,-1], 1, function(x) fzmod(get_abundvec(rev(sort(x[x!=0])), 2000)))
fzT2 <- do.call(rbind, fzT2)

dtgut1 <- read.csv("Desktop/Archive/dtgut1.csv")
dt1 <- apply(dtgut1[complete.cases(dtgut1),-1], 1, function(x) fzmod(get_abundvec(rev(sort(x[x!=0])), 2000)))
dt1 <- do.call(rbind, dt1)

dtgut2 <- read.csv("Desktop/Archive/dtgut2.csv")
dt2 <- apply(dtgut2[complete.cases(dtgut2),-1], 1, function(x) fzmod(get_abundvec(rev(sort(x[x!=0])), 2000)))
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
psd1 <- readRDS("~/Documents/Data/psd1.rds")
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

psd1 <- readRDS("~/Documents/Data/psd1.rds")
psd2 <- readRDS("~/Documents/Data/psd2.rds")
psd3 <- readRDS("~/Documents/Data/psd3.rds")
psd4 <- readRDS("~/Documents/Data/psd4.rds")
psd5 <- readRDS("~/Documents/Data/psd5.rds")
#fzd3 <- t(sapply(psd3$eqa, fzmod))
psd6 <- readRDS("~/Documents/Data/vdat.rds")
psd7 <- readRDS("~/Documents/Data/vdat2.rds")
psd8 <- readRDS("~/Documents/Data/vdat3.rds")
psd9 <- readRDS("~/Documents/Data/vdat5.rds")
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
## Standardize reads to X
X <- 2000

gav1 <- apply(otu3, 2, get_abundvec, N = 2000)
#gav2 <- apply(otu3, 2, function(x) get_abundvec(x[x>=5], N = X))
gavfz1 <- t(sapply(gav1, fzmod))
#gavfz2 <- t(sapply(gav1, function(x) fzmod(x[x >= 3])))
s.hmp <- unlist(gavfz1[,"s"])[which(apply(otu3, 2, sum) > X)]
s.hmp2 <- unlist(gavfz2[,"s"])[which(apply(otu3, 2, sum) > X)]
n.hmp <- unlist(gavfz1[,"N"])[which(apply(otu3, 2, sum) > X)]
n.hmp2 <- unlist(gavfz2[,"N"])[which(apply(otu3, 2, sum) > X)]
r2.hmp <- unlist(gavfz1[,"r2"])[which(apply(otu3, 2, sum) > X)]
hmp1 <- sapply(gav1[apply(otu3, 2, sum) >= 2000], function(x) rev(sort(x)))

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

allfit <- data.frame(s = c(fzT$s, fzT2$s, dt1$s, dt2$s, s.hmp, fzag$s),
                     N = c(fzT$N, fzT2$N, dt1$N, dt2$N, n.hmp, fzag$N),
                     r2 = c(fzT$r2, fzT2$r2, dt1$r2, dt2$r2, r2.hmp, fzag$r2),
                     dat = rep(c("M3", "F4", "DT1", "DT2", "HMP", "AG"), c(length(fzT$s),length(fzT2$s), length(dt1$s), length(dt2$s), length(s.hmp), length(fzag$N))))


ggplot(allfit[allfit$r2 > 0.9,], aes(x = log10(N), y = log10(s))) + geom_point(aes(col = dat)) + geom_smooth(aes(col = "black")) + theme_bw()
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

gavsim4 <- lapply(psd4$eqa, function(x) get_abundvec(x[x>0], N = X))
simfz4 <- lapply(gavsim4, function(x) fzmod(sort(x)))
simfz4 <- do.call(rbind, simfz4)

gavsim5 <- lapply(psd5$eqa, function(x) get_abundvec(x[x>0], N = X))
simfz5 <- lapply(gavsim5, function(x) fzmod(sort(x)))
simfz5 <- do.call(rbind, simfz5)

gavsim6 <- lapply(psd6$eqa, function(x) get_abundvec(x[x>0], N = X))
simfz6 <- lapply(gavsim6, function(x) fzmod(sort(x)))
simfz6 <- do.call(rbind, simfz6)

gavsim7 <- lapply(psd7$eqa, function(x) get_abundvec(x[x>0], N = X))
simfz7 <- lapply(gavsim7, function(x) fzmod(sort(x)))
simfz7 <- do.call(rbind, simfz7)

gavsim8 <- lapply(psd8$eqa, function(x) get_abundvec(x[x>0], N = X))
simfz8 <- lapply(gavsim8, function(x) fzmod(sort(x)))
simfz8 <- do.call(rbind, simfz8)

gavsim9 <- lapply(psd9$eqa, function(x) get_abundvec(x[x>0], N = X))
simfz9 <- lapply(gavsim9, function(x) fzmod(sort(x)))
simfz9 <- do.call(rbind, simfz9)


allfit <- data.frame(s = c(fzT$s, fzT2$s, dt1$s, dt2$s, s.hmp, simfz1$s, simfz2$s, simfz3$s, simfz4$s, simfz5$s, simfz6$s),
                     N = c(fzT$N, fzT2$N, dt1$N, dt2$N, n.hmp, simfz1$N, simfz2$N, simfz3$N, simfz4$N, simfz5$N, simfz6$N),
                     r2 = c(fzT$r2, fzT2$r2, dt1$r2, dt2$r2, r2.hmp, simfz1$r2, simfz2$r2, simfz3$r2, simfz4$r2, simfz5$r2, simfz6$r2),
                     dat = rep(c("M3", "F4", "DT1", "DT2", "HMP", "sim", "sim2", "sim3", "sim4", "sim5", "sim6"),
                               c(length(fzT$s),length(fzT2$s),length(dt1$s),length(dt2$s),length(s.hmp),nrow(simfz1), nrow(simfz2), nrow(simfz3), nrow(simfz4), nrow(simfz5), nrow(simfz6))))

#allfitRe <- readRDS("~/Documents/AbundData/allfitREAL.rds")

p1 <- ggplot(allfitRe, aes(x = log10(N), y = log10(s))) + geom_point(aes(alpha = r2, col = dat)) + geom_smooth(method = "lm", col = "black") + 
  labs(x = expression(log[10](N)), y = expression(log[10](s)), color = "Data Source", alpha = expression("RAD fit" ~ R^2)) + 
  theme_bw() + theme(axis.title=element_text(size=18,face="bold"), axis.text = element_text(size = 14))

dfresid <- data.frame(resids = residuals(lm(log10(s)~log10(N), data = allfitRe)))
p2 <- ggplot(dfresid, aes(x = resids, y = ..density..)) + geom_density(fill = "green4") + 
  labs(x = "Residuals", y = "Density") + 
  theme_bw() + theme(axis.title=element_text(size=18,face="bold"), axis.text = element_text(size = 14))

p3 <- ggplot(alldat2, aes(x = log10(N), y = log10(fz))) + geom_point(aes(alpha = r2, col = typ)) + geom_smooth(method = "lm", col = "black") + 
  labs(x = expression(log[10](N)), y = expression(log[10](s)), color = "Network Type", alpha = expression("RAD fit" ~ R^2)) + 
  scale_color_manual(labels = c("BA", "Random"), values = c("#F8766D", "#00BFC4")) + 
  theme_bw() + theme(axis.title=element_text(size=18,face="bold"), axis.text = element_text(size = 14))

dfresidsim <- data.frame(resids = residuals(lm(log10(fz)~log10(N), data = alldat2)))
p4 <- ggplot(dfresidsim, aes(x = resids, y = ..density..)) + geom_density(fill = "green4") + 
  labs(x = "Residuals", y = "Density") + 
  theme_bw() + theme(axis.title=element_text(size=18,face="bold"), axis.text = element_text(size = 14))


ggplot(allfit, aes(x = r2, y = s, col = dat)) + geom_point() + geom_smooth() + theme_bw() 

library(multipanelfigure)

cols <- 3
rows <- 2
figure <- multi_panel_figure(
  width = 300,
  columns = cols,
  height = 200,
  rows = rows)

figure %<>% fill_panel(p1, row = 1, column = 1:2)
figure %<>% fill_panel(p3, row = 2, column = 1:2)
figure %<>% fill_panel(p2, row = 1, column = 3)
figure %<>% fill_panel(p4, row = 2, column = 3)
figure
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
#### Effect of truncation on s value
strt <- Sys.time()
trunc <- seq(10, 400, 10)
hmp.sran <- list()
sim.sran1 <- list()
sim.sran2 <- list()
for(i in 1:length(trunc)){
  hmp.sran[[i]] <- sapply(hmp1[sapply(hmp1, length) >= trunc[i]], function(x) fitzipf_r(x[1:trunc[i]])@coef)
  sim.sran1[[i]] <- sapply(gavsim[sapply(gavsim, length) >= trunc[i]], function(x) fitzipf_r(x[1:trunc[i]])@coef)
  sim.sran2[[i]] <- sapply(gavsim4[sapply(gavsim4, length) >= trunc[i]], function(x) fitzipf_r(x[1:trunc[i]])@coef)
}
fin <- Sys.time()
fin-strt

length(hmp.sran)
hist(hmp.sran[[40]])

plot(sapply(hmp.sran, median)~trunc, ylim = c(0,1.5), typ = "l", lwd = 2)
points(sapply(hmp.sran, quantile, prob = 0.025)~trunc, typ = "l", lwd = 2)
points(sapply(hmp.sran, quantile, prob = 0.975)~trunc, typ = "l", lwd = 2)

points(sapply(sim.sran1, median)~trunc, typ = "l", lty = 2, col = "blue", lwd = 2)
points(sapply(sim.sran1, quantile, prob = 0.025)~trunc, typ = "l", lty = 2, col = "blue", lwd = 2)
points(sapply(sim.sran1, quantile, prob = 0.975)~trunc, typ = "l", lty = 2, col = "blue", lwd = 2)

points(sapply(sim.sran2, median)~trunc, typ = "l", lty = 2, col = "darkgreen", lwd = 2)
points(sapply(sim.sran2, quantile, prob = 0.025)~trunc, typ = "l", lty = 2, col = "darkgreen", lwd = 2)
points(sapply(sim.sran2, quantile, prob = 0.975)~trunc, typ = "l", lty = 2, col = "darkgreen", lwd = 2)

####################################################################################################################################
####################################################################################################################################
#### Check whether real data is within range 
# this function replicates everything in this section
inrange <- function(trunc, gav, realdat){
  realdat <- lapply(realdat, function(x) rev(sort(x)))
  gav <- lapply(gav, function(x) rev(sort(x)))
  
  realrange <- sapply(realdat[sapply(realdat, length) >= trunc], function(x) fitzipf_r(x[1:trunc])@coef)
  simrange <- sapply(gav[sapply(gav, length) >= trunc], function(x) fitzipf_r(x[1:trunc])@coef)
  
  irT <- vector(length = length(gav)) 
  irT[sapply(gav, length) >= trunc] <- simrange > min(realrange) & simrange < max(realrange)
  irT[sapply(gav, length) < trunc] <- NA
  
  return(irT)
}

s50range <- sapply(gav1[sapply(gav1, length) >= 50], function(x) fitzipf_r(x[1:50])@coef)
s50sim <- sapply(gavsim4[sapply(gavsim, length) >= 50], function(x) fitzipf_r(x[1:50])@coef)
s50sim4 <- sapply(gavsim4[sapply(gavsim4, length) >= 50], function(x) fitzipf_r(x[1:50])@coef)
s50sim6 <- sapply(gavsim6[sapply(gavsim6, length) >= 50], function(x) fitzipf_r(x[1:50])@coef)

hist(s50sim4, freq = F, xlim = c(0, 2), col = "grey")
hist(s50range, freq = F, add = T)

sum(s50sim4 > min(s50range) & s50sim4 < max(s50range))
ir50.4 <- vector(length = length(gavsim4)) 
ir50.4[sapply(gavsim4, length) >=50] <- s50sim4 > min(s50range) & s50sim4 < max(s50range)
ir50.4[!sapply(gavsim4, length) >=50] <- NA

ir50.6 <- vector(length = length(gavsim6)) 
ir50.6[sapply(gavsim6, length) >=50] <- s50sim6 > min(s50range) & s50sim6 < max(s50range)
ir50.6[!sapply(gavsim6, length) >=50] <- NA
  
s100range <- sapply(hmp1[sapply(hmp1, length) >= 100], function(x) fitzipf_r(x[1:100])@coef)
s100sim <- sapply(gavsim4[sapply(gavsim, length) >= 100], function(x) fitzipf_r(x[1:100])@coef)
s100sim4 <- sapply(gavsim4[sapply(gavsim4, length) >= 100], function(x) fitzipf_r(x[1:100])@coef)
s100sim6 <- sapply(gavsim6[sapply(gavsim6, length) >= 100], function(x) fitzipf_r(x[1:100])@coef)

hist(s100sim4, freq = F, xlim = c(0, 2), col = "grey")
hist(s100range, freq = F, add = T)

sum(s100sim4 > min(s100range) & s100sim4 < max(s100range))
ir100.4 <- vector(length = length(gavsim4)) 
ir100.4[sapply(gavsim4, length) >=100] <- s100sim4 > min(s100range) & s100sim4 < max(s100range)
ir100.4[!sapply(gavsim4, length) >=100] <- NA

ir100.6 <- vector(length = length(gavsim6)) 
ir100.6[sapply(gavsim6, length) >= 100] <- s100sim6 > min(s100range) & s100sim6 < max(s100range)
ir100.6[!sapply(gavsim6, length) >= 100] <- NA

s200range <- sapply(hmp1[sapply(hmp1, length) >= 200], function(x) fitzipf_r(x[1:200])@coef)
s200sim <- sapply(gavsim4[sapply(gavsim, length) >= 200], function(x) fitzipf_r(x[1:200])@coef)
s200sim4 <- sapply(gavsim4[sapply(gavsim4, length) >= 200], function(x) fitzipf_r(x[1:200])@coef)
s200sim6 <- sapply(gavsim6[sapply(gavsim6, length) >= 200], function(x) fitzipf_r(x[1:200])@coef)

hist(s200sim4, freq = F, xlim = c(0, 1.6), col = "grey")
hist(s200range, freq = F, add = T)

sum(s200sim4 > min(s200range) & s200sim4 < max(s200range))
ir200.4 <- vector(length = length(gavsim4)) 
ir200.4[sapply(gavsim4, length) >=200] <- s200sim4 > min(s200range) & s200sim4 < max(s200range)
ir200.4[!sapply(gavsim4, length) >=200] <- NA

ir200.6 <- vector(length = length(gavsim6)) 
ir200.6[sapply(gavsim6, length) >= 200] <- s200sim6 > min(s200range) & s200sim6 < max(s200range)
ir200.6[!sapply(gavsim6, length) >= 200] <- NA
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
### Interactions
#testing objs
# add self interaction mean

pdat1 <- prepdat(eqa = gavsim2, eqm = psd1$eqm, svals = simfz2$s, sr2 = simfz2$r2, d = psd1$ds, r = psd1$rs)
pdat1$k <- 20
pdat1$ir50 <- inrange(50, gavsim2, hmp1)
pdat1$ir100 <- inrange(100, gavsim2, hmp1)
pdat1$ir200 <- inrange(200, gavsim2, hmp1)

pdat2 <- prepdat(eqa = gavsim, eqm = psd2$eqm, svals = simfz1$s, sr2 = simfz1$r2, d = psd2$ds, r = psd2$rs)
pdat2$k <- 20
pdat2$ir50 <- inrange(50, gavsim, hmp1)
pdat2$ir100 <- inrange(100, gavsim, hmp1)
pdat2$ir200 <- inrange(200, gavsim, hmp1)

pdat3 <- prepdat(eqa = gavsim3, eqm = psd3$eqm, svals = simfz3$s, sr2 = simfz3$r2, d = psd3$ds, r = psd3$rs)
pdat3$k <- 20
pdat3$ir50 <- inrange(50, gavsim3, hmp1)
pdat3$ir100 <- inrange(100, gavsim3, hmp1)
pdat3$ir200 <- inrange(200, gavsim3, hmp1)

pdat4 <- prepdat(eqa = gavsim4, eqm = psd4$eqm, svals = simfz4$s, sr2 = simfz4$r2, d = psd4$ds, r = psd4$rs)
pdat4$k <- psd4$kv
pdat4$ir50 <- inrange(50, gavsim4, hmp1)
pdat4$ir100 <- inrange(100, gavsim4, hmp1)
pdat4$ir200 <- inrange(200, gavsim4, hmp1)

pdat5 <- prepdat(eqa = gavsim5, eqm = psd5$eqm, svals = simfz5$s, sr2 = simfz5$r2, d = psd5$ds, r = psd5$rs)
pdat5$k <- psd5$kv
pdat5$ir50 <- inrange(50, gavsim5, hmp1)
pdat5$ir100 <- inrange(100, gavsim5, hmp1)
pdat5$ir200 <- inrange(200, gavsim5, hmp1)

pdat6 <- prepdat(eqa = gavsim6, eqm = psd6$eqm, svals = simfz6$s, sr2 = simfz6$r2, d = psd6$ds, r = psd6$rs)
pdat6$k <- psd6$kv
pdat6$ir50 <- inrange(50, gavsim6, hmp1)
pdat6$ir100 <- inrange(100, gavsim6, hmp1)
pdat6$ir200 <- inrange(200, gavsim6, hmp1)

pdat7 <- prepdat(eqa = gavsim7, eqm = psd7$eqm, svals = simfz7$s, sr2 = simfz7$r2, d = psd7$ds, r = psd7$rs)
pdat7$k <- psd7$kv

pdat8 <- prepdat(eqa = gavsim8, eqm = psd8$eqm, svals = simfz8$s, sr2 = simfz8$r2, d = psd8$ds, r = psd8$rs)
pdat8$k <- psd8$kv

pdat9 <- prepdat(eqa = gavsim9, eqm = psd9$eqm, svals = simfz9$s, sr2 = simfz9$r2, d = psd9$ds, r = psd9$rs)
pdat9$k <- psd9$kv
pdat9$ir50 <- inrange(50, gavsim9, hmp1)

alldat <- data.table::rbindlist(list(pdat1,pdat2,pdat4,pdat6,pdat7,pdat8,pdat9))

subdat <- apply(abs(pdat2[,-c(1,2)]), 2, function(x){(x - mean(x))/sd(x)})
subdat <- data.frame(pdat2[,c(1,2)], subdat)
subdat2 <- subdat[pdat2$Nsp > 50,]
subdat2$ir <- ir
subdat2$Nsp <- subdat$Nsp[pdat2$Nsp > 50]
subdat2$Nsp2 <- pdat2$Nsp[pdat2$Nsp > 50]

subdata <- apply(abs(pdat3[,-c(1,2)]), 2, function(x){(x - mean(x))/sd(x)})
subdata <- data.frame(pdat3[,c(1,2)], subdata)

x1 <- which(pdat7[,"sR"] > 0.99)
x <- x1[79]
#plot(dzipf(1:pdat7[x,"Nsp"], pdat7[x,"Nsp"], 1), typ = "o", col = "green4")
#points(dzipf(1:pdat7[x,"Nsp"], pdat7[x,"Nsp"], .8), typ = "o", col = "green3")
plot(dzipf(1:pdat7[x,"Nsp"], pdat7[x,"Nsp"], pdat7[x,"sV"]), typ = "o", col = "blue")
points(rev(sort(gavsim7[[x]]))/sum(gavsim7[[x]]))
pdat7[x,"sR"]

plot(aggregate(pdat9$Nsp, list(pdat9$k), quantile, prob = 0.975), ylim = c(0,600), typ = "o")
points(aggregate(pdat9$Nsp, list(pdat9$k), median), typ = "o")
points(aggregate(pdat9$Nsp, list(pdat9$k), quantile, prob = 0.025), typ = "o")
points(aggregate(pdat7$Nsp, list(pdat7$k), max), typ = "o")
points(aggregate(pdat7$Nsp, list(pdat7$k), min), typ = "o")

df1 <- data.frame(N = signif(pdat9$Nsp,1)[pdat9$Nsp > 10], K = pdat9$k[pdat9$Nsp > 10], s = pdat9$sV[pdat9$Nsp > 10])
ggplot(df1, aes(x = (K), y = s)) + geom_point() + geom_smooth() + facet_wrap(~N, scales = "free_y")
ggplot(df1, aes(x = factor(K), y = s)) + geom_boxplot() + facet_wrap(~N, scales = "free_y")

#######################################
# Models ##############################
#######################################

# all orig data modeling s value
fit1 <- (lm(sV~Nsp+r+C+k+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = pdat9, model = F, x = F, y = F))
summary(fit1)
# same but for new sim, with K
fit1a <- (lm(sV~Nsp, data = pdat9))
summary(fit1a)
# same but for new sim, with K
fit1b <- (lm(sV~Nsp+r+k+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = pdat4))
summary(fit1b)
# orig data (good fit) modeling s value
fit2 <- (lm(sV~Nsp+r+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = pdat2[pdat2$sR > 0.8,]))
summary(fit2)
# same but for new sim, with K
fit2a <- (lm(sV~Nsp+r+k+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = pdat3[pdat3$sR > 0.8,]))
summary(fit2a)
# all orig data modeling r2
fit3 <- (lm(sR~Nsp+r+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = pdat2))
summary(fit3)
# all orig data modeling r2, with K
fit3a <- (lm(sR~Nsp+r+k+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = pdat3))
summary(fit3a)

#######################################
#######################################
#######################################
tot_count <- function(x, labs, digits, varlen)
{
  paste(labs, "\n\nn =", x$frame$n)
}

fitir <- glm(ir200~cpN+coN, data = pdat3[complete.cases(pdat3),], family = "binomial")
summary(fitir)

fpart50 <- rpart(sV~r+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = alldat, method = "class")
fpart100 <- rpart(ir100~r+k+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = pdat9, method = "class")
fpart200 <- rpart(ir200~r+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = pdat9, method = "class")
plotcp(fpart50)
rpart.plot::prp(fpart50, uniform = T, node.fun = tot_count, extra = 1)
rpart.plot::prp(fpart100, uniform = T, node.fun = tot_count, extra = 1)
rpart.plot::prp(fpart200, uniform = T, node.fun = tot_count, extra = 1)
?rpart
#plot(fpart, uniform = T)
#text(fpart, cex = 0.75)
library(randomForest)
fit50 <- randomForest(factor(ir50)~r+k+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = alldat[!is.na(alldat$ir50),], ntree = 3000, importance = T, mtry = 15)
fit100 <- randomForest(factor(ir100)~r+k+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = alldat[!is.na(alldat$ir100),], ntree = 3000, importance = T, mtry = 15)

fit200 <- randomForest(factor(ir200)~r+k+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = alldat[!is.na(alldat$ir200),], ntree = 3000, importance = T, mtry = 15)
zs <- sample(which(alldat$ir200[!is.na(alldat$ir200)] != 1), 100)
ons <- which(alldat$ir200[!is.na(alldat$ir200)] != 0)
fit200A <- randomForest(factor(ir200)~r+k+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = alldat[!is.na(alldat$ir200),][c(zs,ons),], ntree = 3000, importance = T, mtry = 10)
varImpPlot(fit200A)
fit200Acart <- rpart(factor(ir200)~r+k+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = alldat[!is.na(alldat$ir200),][c(zs,ons),], method = "class")
rpart.plot::prp((fit200Acart), uniform = T, node.fun = tot_count, extra = 1)
plotcp(fit200Acart)

print(fit50)
print(fit100)
print(fit200)
varImpPlot(fit50)
varImpPlot(fit100)
varImpPlot(fit200)

rfS <- randomForest(maxn~Nsp+r+k+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = pdat9, ntree = 3000, importance = T)
print(rfS)
rfSt <- rpart::rpart(sV~Nsp+r+k+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = pdat9, method = "anova")
rpart.plot::prp(rfSt, uniform = T, extra = 1)
varImpPlot(rfS)
rpart::rsq.rpart(rfSt)
#######################################
#######################################
#######################################




#######################################
#######################################
#######################################
# all data as zscore modeling s value
fit1.1 <- (lm(sV~Nsp+r+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = subdat, na.action = "na.fail", model = F, x = F, y = F))
summary(fit1.1)
# all data as zscore modeling s value
fit1.1a <- (lm(sV~Nsp+r+k+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = subdata, na.action = "na.fail", model = F, x = F, y = F))
summary(fit1.1a)
# zscore data (good fit) modeling s value
fit2.1 <- (lm(sV~Nsp+r+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = subdat[subdat$sR > 0.8,]))
summary(fit2.1)
# zscore data modeling r2
fit3.1 <- (lm(sR~Nsp+r+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = subdat, model = F, x = F, y = F, na.action = "na.fail"))
summary(fit3.1)

fit4 <- lm(cbind(sV, Nsp, abs)~r+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = subdat)
summary(fit4)

subdat3 <- data.frame(mu = t(sa2)[,1], subdat)
fit5 <- (lm(mu~Nsp+C+r+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = subdat3))
summary(fit5)
data.frame(smod = summary(fit1)$coefficients[,1], rmod = summary(fit3)$coefficients[,1], p1 = summary(fit1)$coefficients[,4] <= 0.05, p2 = summary(fit3)$coefficients[,4] <= 0.05)

library(MuMIn)
allmod <- dredge(fit1.1)
allmodR <- dredge(fit3.1)
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

pall <- rbind(pdat1, pdat2[,-18])
sall <- apply(pall[,-c(1,2)], 2, function(x){(x - mean(x))/sd(x)})
sall <- data.frame(pall[,c(1,2)], sall)

fitX <- lm(sV~Nsp+C+r+D+aN+coN+cpN+mN+pN+aS+coS+cpS+mS+pSn+pSp, data = sall) %>% summary()

data.frame(smod = summary(fit1.1)$coefficients[,1], rmod = summary(fit3.1)$coefficients[,1], amod = fitX$coefficients[,1], 
           p1 = summary(fit1.1)$coefficients[,4] <= 0.05, p2 = summary(fit3.1)$coefficients[,4] <= 0.05, p3 = fitX$coefficients[,4] <= 0.05)

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

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
### plot all rads

gav1.1 <- lapply(gavsim, function(x) rev(sort(x)))
max(unlist(gav1.1))
max(sapply(gav1.1, length))

plot(gav1.1[[1]], ylim = c(0, 1030), xlim = c(0, 800))
lapply(gav1.1[-1], function(x) points(x))

gav1.2 <- lapply(gav1, function(x) rev(sort(x)))
lapply(gav1.2, function(x) points(x, pch = 20, col = "blue"))
