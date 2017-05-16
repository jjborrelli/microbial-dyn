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


X <- 2000
otu2 <- read.csv("~/Desktop/otu_table_psn_v13.csv", row.names = 1)
metadat <- read.csv("~/Desktop/v13_map_uniquebyPSN.csv")

stoolsamp <- which(metadat$HMPbodysubsite == "Stool")
spptab <- colnames(otu2) %in% paste0("X",metadat[stoolsamp,]$SampleID)
otu3 <- otu2[-which(rowSums(otu2[,spptab]) == 0),spptab]
rm(otu2)
rm(metadat)

gav1 <- apply(otu3, 2, get_abundvec, N = X)
gavfz1 <- t(sapply(gav1, fzmod))
s.hmp <- unlist(gavfz1[,"s"])[which(apply(otu3, 2, sum) > X)]
n.hmp <- unlist(gavfz1[,"N"])[which(apply(otu3, 2, sum) > X)]
hmp1 <- sapply(gav1[apply(otu3, 2, sum) >= 2000], function(x) rev(sort(x)))


psd9 <- readRDS("~/Desktop/vdat5.rds")
gavsim9 <- lapply(psd9$eqa, function(x) get_abundvec(x[x>0], N = X))
simfz9 <- lapply(gavsim9, function(x) fzmod(sort(x)))
simfz9 <- do.call(rbind, simfz9)

pdat9 <- prepdat(eqa = gavsim9, eqm = psd9$eqm, svals = simfz9$s, sr2 = simfz9$r2, d = psd9$ds, r = psd9$rs)
pdat9$k <- psd9$kv
pdat9$ir50 <- inrange(50, gavsim9, hmp1)
pdat9$ir100 <- inrange(100, gavsim9, hmp1)
pdat9$ir200 <- (inrange(200, gavsim9, hmp1))


######################################################################################################
######################################################################################################
#### Machine Learning with h2o
library(h2o)
library(tidyverse)

localh2o <- h2o.init()

pdat.h2o <- as.h2o(pdat9[pdat9$Nsp >= 200,-c(1,2,14,19,20)])
splits <- h2o.splitFrame(pdat.h2o, ratios = c(0.4, 0.4), seed = 42)

train_unsupervised  <- splits[[1]]
train_supervised  <- splits[[2]]
test <- splits[[3]]

response <- "ir200"
features <- setdiff(colnames(train_unsupervised), response)
strt <- Sys.time()

# Autoencoders
model_nn <- h2o.deeplearning(x = features,
                             training_frame = train_unsupervised,
                             model_id = "model_nn",
                             autoencoder = TRUE,
                             reproducible = TRUE, #slow - turn off for real problems
                             ignore_const_cols = FALSE,
                             seed = 42,
                             hidden = c(10, 2, 10), 
                             epochs = 100,
                             activation = "Tanh")
done <- Sys.time()
done-strt

test_autoenc <- h2o.predict(model_nn, test)

# Dimensionality reduction with hidden layers
train_features <- h2o.deepfeatures(model_nn, train_unsupervised, layer = 2) %>%
  as.data.frame() %>%
  mutate(Class = as.vector(train_unsupervised[, 16]))

ggplot(train_features, aes(x = DF.L2.C1, y = DF.L2.C2, color = Class)) + geom_point()

train_features <- h2o.deepfeatures(model_nn, train_unsupervised, layer = 3) %>%
  as.data.frame() %>%
  mutate(ir200 = as.factor(as.vector(train_unsupervised[, 16]))) %>%
  as.h2o()
features_dim <- setdiff(colnames(train_features), response)

model_nn_dim <- h2o.deeplearning(y = response,
                                 x = features_dim,
                                 training_frame = train_features,
                                 reproducible = TRUE, #slow - turn off for real problems
                                 balance_classes = TRUE,
                                 ignore_const_cols = FALSE,
                                 seed = 42,
                                 hidden = c(10, 2, 10), 
                                 epochs = 100,
                                 activation = "Tanh")

test_dim <- h2o.deepfeatures(model_nn, test, layer = 3)

h2o.predict(model_nn_dim, test_dim) %>%
  as.data.frame() %>%
  mutate(actual = as.vector(test[, 16])) %>%
  group_by(actual, predict) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

# Anomaly detection

anomaly <- h2o.anomaly(model_nn, test) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  mutate(ir200 = as.vector(test[, 16]))

mean_mse <- anomaly %>%
  group_by(ir200) %>%
  summarise(mean = mean(Reconstruction.MSE))

ggplot(anomaly, aes(x = as.numeric(rowname), y = Reconstruction.MSE, color = as.factor(ir200))) +
  geom_point() +
  geom_hline(data = mean_mse, aes(yintercept = mean, color = as.factor(ir200))) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "instance number",
       color = "ir200")

anomaly <- anomaly %>%
  mutate(outlier = ifelse(Reconstruction.MSE > 0.02, "outlier", "no_outlier"))

anomaly %>%
  group_by(ir200, outlier) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 

# Pre-trained supervised model

model_nn_2 <- h2o.deeplearning(y = response,
                               x = features,
                               training_frame = train_supervised,
                               pretrained_autoencoder  = "model_nn",
                               reproducible = TRUE, #slow - turn off for real problems
                               balance_classes = TRUE,
                               ignore_const_cols = FALSE,
                               seed = 42,
                               hidden = c(10, 2, 10), 
                               epochs = 100,
                               activation = "Tanh")

pred <- as.data.frame(h2o.predict(object = model_nn_2, newdata = test)) %>%
  mutate(actual = as.vector(test[, 16]))

pred %>%
  group_by(actual, predict) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 

pred %>%
  ggplot(aes(x = actual, fill = predict)) +
  geom_bar() +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap( ~ actual, scales = "free", ncol = 2)