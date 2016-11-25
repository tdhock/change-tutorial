load("breakpoint.learning.RData")
library(doParallel)
registerDoParallel()
library(iregnet)

set.name <- "original"
inf.target.mat <- breakpoint.learning$targets[[set.name]]
has.info <- 0 < rowSums(is.finite(inf.target.mat))
target.mat <- inf.target.mat[has.info, ]
test.fold.vec <- sub(".*[.]", "", rownames(target.mat))

for(test.fold in unique(test.fold.vec)){
  is.test <- test.fold.vec == test.fold
  train.target.mat <- target.mat[!is.test,]
  train.feature.mat <- breakpoint.learning$features[rownames(train.target.mat), ]
  fit <- cv.iregnet(train.feature.mat, train.target.mat, family="gaussian")
  for(selected.name in names(fit$selected)){
    selected.i <- fit$selected[[selected.name]]
    pred.mat <- cbind(1, feature.mat) %*% fit$beta[, i, drop=FALSE]
    is.lo <- pred.mat < target.mat[,1]
    is.hi <- target.mat[,2] < pred.mat
    stop("TODO compute test error")
  }
}

