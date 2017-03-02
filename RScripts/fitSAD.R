load("~/Downloads/sim.Rdata")


plot(hist(sort(ge.mult2$eqst[[1]], decreasing = T), plot = F)$counts)

igraph::fit_power_law(hist(sort(ge.mult2$eqst[[1]], decreasing = T), plot = F)$counts)

table(ge.mult2$eqst[[1]])

fpl3[[1]]

par(mfrow = c(5,5))
for(i in 1:25){
  plot(sort(ge.mult2$eqst[[i]], decreasing = T))
  text(100,.1, label = sapply(fpl2,"[[",5))
}

dev.off()


plot(as.vector(table(ceiling(ge.mult2$eqst[[1]]/min(ge.mult2$eqst[[1]]))))~as.numeric(names(table(ceiling(ge.mult2$eqst[[1]]/min(ge.mult2$eqst[[1]]))))))
hist(sort(ge.mult2$eqst[[1]], decreasing = T))


length(ge.mult2$eqst[[]])
samp1 <- sample(1:327, 100000, prob = ge.mult2$eqst[[1]]/sum(ge.mult2$eqst[[1]]), replace = T)
plot(as.vector(table(as.vector(table(samp1))))~as.numeric(names((table(as.vector(table(samp1)))))))

abund <- ge.mult3$eqst[[1]]
r.ab <- abund/sum(abund)
N = c(10^3, 5*10^3, 10^4, 5*10^4, 10^5)
par(mfrow = c(5,1))
t3 <- list()
for(i in 1:5){
  samp2 <- sample(1:length(abund), N[i], replace = T, prob = r.ab)
  t1 <- table(samp2)
  t2 <- table(as.vector(t1))
  plot(as.vector(t2)~as.numeric(names(t2)), main = N[i], xlab = "Abundance", ylab = "Freq")
  
  t3[[i]] <- as.numeric(t1)
}

dev.off()
plot(r.ab[1:length(r.ab) %in% as.numeric(names(t1))], as.vector(t1)/sum(as.vector(t1)))
  
lapply(t3, fitsad, "pareto")
