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


############################################################################################################
############################################################################################################
############################################################################################################

args <- commandArgs(TRUE)

st1 <- Sys.time()
psd1 <- get_dat(args)
st2 <- Sys.time()
st2-st1

gavsim <- lapply(psd1$eqa, get_abundvec, X)
simfz1 <- lapply(gavsim, function(x) fzmod(sort(x)))
simfz1 <- do.call(rbind, simfz1)

pdat2 <- prepdat(eqa = gavsim, eqm = psd2$eqm, svals = simfz1$s, sr2 = simfz1$r2, d = psd2$ds, r = psd2$rs)
pdat2$k <- psd1$kv

saveRDS(pdat2, file = "~/Abundance/vdat.rds")