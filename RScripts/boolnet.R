library(BoolNet)
library(arules)
library(arulesViz)

data("yeastTimeSeries")
yeastTimeSeries

bin <- binarizeTimeSeries(yeastTimeSeries)
renet <- reconstructNetwork(bin$binarizedMeasurements, method = "bestfit", maxK = 4)
chooseNetwork(renet, c(1,1,1,1))
plotNetworkWiring(chooseNetwork(renet, c(2,2,2,2)))

eg <- expand.grid(1:2, 1:2, 1:3, 1:2)
cNet <- list()
for(i in 1:nrow(eg)){
  cNet[[i]] <- chooseNetwork(renet, eg[i,])
}

rule1 <- apriori(t(bin$binarizedMeasurements))
summary(rule1)
inspect(head(rule1))
plot(rule1)
inspect(head(rule1, n = 3, by = "lift"))

rules.sorted <- sort(rule1, by="lift")
# find redundant rule
subset.matrix <- is.subset(rules.sorted, rules.sorted)
subset.matrix[lower.tri(subset.matrix, diag=T)] <- NA
redundant <- colSums(subset.matrix, na.rm=T) >= 1

# remove redundant rules
rules.pruned <- rules.sorted[!redundant]
inspect(rules.pruned)



install.packages("bnlearn")
library(bnlearn)
plot(hc(data.frame(t(yeastTimeSeries))))

