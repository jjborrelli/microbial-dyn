library(igraph)

# compute frequency of different interaction types for whole matrix
itypes <- function(x){
  i1 <- x[upper.tri(x)]
  i2 <- t(x)[upper.tri(x)] 
  
  comp <- sum(i1 < 0 & i2 < 0)
  mut <- sum(i1 > 0 & i2 > 0)
  pred <- sum(i1 > 0 & i2 < 0 | i1 < 0 & i2 > 0)
  amens <- sum(i1 < 0 & i2  == 0 | i1 == 0 & i2 < 0)
  comm <- sum(i1 > 0 & i2  == 0 | i1 == 0 & i2 > 0)
  
  return(c(comp = comp, mut = mut, pred = pred, amens = amens, comm = comm))
}


efun <- function(x, mats, iter = 1000){
  
  e.val <- c()
  for(i in 1:iter){
    m <- mats[[x]]
    m[m < 0] <- runif(sum(m < 0), -1, 0)
    m[m > 0] <- runif(sum(m > 0), 0, 1)
    diag(m) <- -1
    
    e.val[i] <- max(Re(eigen(m)$values)) 
  }
  return(e.val)
}


N = sample(10:50, 10000, replace = T)
C = sample(seq(.1, .8, .1), 10000, replace = T)


g <- lapply(1:10000, function(x) erdos.renyi.game(N, C, "gnp", directed = T))
a1 <- lapply(g, get.adjacency, sparse = F)
a2 <- lapply(a1, function(x){x[x != 0] <- sample(c(-1,1), sum(x), replace = T);return(x)})

ity <- t(sapply(a2, itypes))

eigs <- lapply(1:length(a2), function(x) efun(x, a2, 1000))
qss <- sapply(eigs, function(x) c(sum(x > 0), sum(x <= 0)))

summary(glm(t(qss)~ity[,1], family = "binomial"))
summary(glm(t(qss)~ity[,2], family = "binomial"))
summary(glm(t(qss)~ity[,3], family = "binomial"))

summary(glm(t(qss)~N+C+ity, family = "binomial"))
summary(glm(t(qss)~ity, family = "binomial"))
