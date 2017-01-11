max(Re(eigen(jacobian.full(eqab[[1]], lvmod, parms = list(alpha = growth[eqcomm[[1]]], m = matuse[[1]])))$values))

eig <- c()
x <- 1
matA <- jacobian.full(eqab[[x]], lvmod, parms = list(alpha = growth[eqcomm[[x]]], m = matuse[[x]]))
for(i in 1:1000){
  mat1 <- matA
  int.n <- sum(abs(sign(mat1)))
  mat1[mat1!=0] <- rnorm(int.n, mat1[mat1!=0], .5)
  eig[i] <- max(Re(eigen(mat1)$values))
}

hist(eig)
sum(eig < 0)




eqab <- lapply(dyn, function(x) x[1000,-1][x[1000,-1] > 0])
espp <- lapply(dyn, function(x) which(x[1000,-1] > 0))

gs1 <- lapply(matuse, function(x) graph.adjacency(abs(sign(x))))

get.eig <- function(eq, mat, growth, dev){
  eig <- c()
  matA <- jacobian.full(eq, lvmod, parms = list(alpha = growth, m = mat))
  for(i in 1:1000){
    mat1 <- matA
    int.n <- sum(abs(sign(mat1)))
    mat1[mat1!=0] <- rnorm(int.n, mat1[mat1!=0], dev)
    eig[i] <- max(Re(eigen(mat1)$values))
  }
  
  return(eig)
}

qss.rem <- c()
x <- 3
for(i in 1:length(eqab[[x]])){
  e1 <- get.eig(eq = eqab[[x]][-i], mat = matuse[[x]][-i,-i], growth = growth[espp[[x]]][-i], dev = .1)
  qss.rem[i] <- sum(e1 < 0)/length(e1)
}

iqss <- c()
for(x in 1:sum(use)){
  e1 <- get.eig(eq = eqab[[x]], mat = matuse[[x]], growth = growth[espp[[x]]], dev = .1)
  iqss[x] <- sum(e1 < 0)/length(e1)
}
boxplot(iqss)

q <- (qss.rem-iqss)
plot(rowSums(!ks3[[x]], na.rm = T)~q)
abline(v = iqss)


qss <- lapply(1:sum(use), function(x) get.eig(eq = eqab[[x]], mat = matuse[[x]], growth = growth[espp[[x]]], dev = .1))
qss1 <- sapply(qss, function(x) sum(x < 0)/length(x))

plot(unlist(lapply(ks1, function(x) x[,4]/nrow(x)))~unlist(lapply(gs1, betweenness)))
