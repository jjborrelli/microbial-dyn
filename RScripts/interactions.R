library(MASS)
m3l2 <- read.csv("Data/Data/L2.csv", row.names = 1)

r.m3l2 <- t(apply(m3l2, 2, function(x) x/sum(x)))
r.m3l2A <- t(apply(m3l2, 2, function(x) x/sum(x)))[,which(apply(r.m3l2, 2, function(x) sum(x != 0)) >300)]

v1 <- log(r.m3l2A[2:332,1]) - log(r.m3l2A[1:331,1])
M <- apply(r.m3l2A, 2, function(x) x - median(x))

mydat <- data.frame(M[1:331,], resp = v1)
colnames(mydat) <- c(LETTERS[1:6], "resp")
mydat$resp <- log(mydat$resp) - log(mydat$A)

stepAIC(lm(resp~A, data = mydat),list(lower = resp~A, upper= resp~A+B+C+D+E+F), direction = "forward", data = mydat)


hlf <- length(1:332)/2
t1 <- sample(2:332, hlf)
t2 <- t1-1

fit <- lm(resp~A, data = mydat[t1,])
error1 <- sum(predict(fit, newdata = mydat[-t1,]) - mydat[-t1,"resp"])

fit1 <- lm(resp~A+F, data = mydat[t1,])
error2 <- sum(predict(fit1, newdata = mydat[-t1,]) - mydat[-t1, "resp"])

((error1-error2)/error2)*100
