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
