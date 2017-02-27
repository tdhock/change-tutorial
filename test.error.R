source("packages.R")

## Load neuroblastoma data set and pre-computed optimal Gaussian
## segmentation models.
data(neuroblastoma, package="neuroblastoma")
load("Segmentor.models.RData")

## Compute model selection functions, label error, and target
## intervals.
selection <- Segmentor.models$loss[, modelSelection(
  .SD, complexity="n.segments"
  ), by=.(profile.id, chromosome)]
changes <- Segmentor.models$segs[1 < start, ]
errors <- labelError(
  selection, neuroblastoma$annotations, changes,
  change.var="chromStart",
  label.vars=c("min", "max"),
  model.vars="n.segments",
  problem.vars=c("profile.id", "chromosome"))
target.dt <- targetIntervals(
  errors$model.errors, c("profile.id", "chromosome"))
errors$model.errors[, pid.chr := paste0(profile.id, ".", chromosome)]
setkey(errors$model.errors, pid.chr)

## Make sure that target.dt is the same as the pre-computed target
## matrix in the neuroblastomaProcessed data set.
data(neuroblastomaProcessed, package="penaltyLearning")
all.pid.chr <- target.dt[, paste0(profile.id, ".", chromosome)]
processed.targets <- neuroblastomaProcessed$target.mat[all.pid.chr, ]
computed.targets <- target.dt[, cbind(min.log.lambda, max.log.lambda)]
dimnames(computed.targets) <- dimnames(processed.targets)
stopifnot(all.equal(computed.targets, processed.targets))
processed.features <- neuroblastomaProcessed$feature.mat[all.pid.chr, ]

## Each of these functions should return a list which contains at
## least one element, named predict, which is function that takes a N
## x P feature matrix and returns an N-vector of log(penalty) values.
model.fun.list <- list(
  BIC=function(){
    list(predict=function(X.pred){
      X.pred[, "log2.n"]
    })
  }, cghseg.k=function(){
    pred.dt <- data.table(
      target.dt[is.train, ],
      pred.log.lambda=train.feature.mat[, "log.n"])
    roc <- ROChange(errors$model.errors, pred.dt, c("profile.id", "chromosome"))
    min.train.error <- roc$thresholds[threshold=="min.error",]
    learned.thresh <- min.train.error[, (min.thresh+max.thresh)/2]
    list(predict=function(X.pred){
      X.pred[, "log.n"] + learned.thresh
    })
  }, log.n=function(){
    IntervalRegressionUnregularized(
      train.feature.mat[, "log.n", drop=FALSE], train.target.mat)
  })

set.seed(1)
test.error.list <- list()
## 6-fold CV, one for for each chrom.
test.fold.vec <- paste(target.dt$chromosome)
for(test.fold in unique(test.fold.vec)){
  is.train <- test.fold.vec != test.fold
  train.target.mat <- processed.targets[is.train, ]
  train.feature.mat <- processed.features[is.train, ]
  test.feature.mat <- processed.features[!is.train, ]
  for(model.name in names(model.fun.list)){
    model.fun <- model.fun.list[[model.name]]
    model <- model.fun()
    pred.dt <- data.table(
      target.dt[!is.train,],
      pred.log.lambda=as.numeric(model$predict(test.feature.mat)))
    test.roc <- ROChange(
      errors$model.errors, pred.dt, c("profile.id", "chromosome"))
    model.test.error <- with(test.roc, {
      data.table(
        test.fold, model.name, auc,
        error.percent=thresholds[threshold=="predicted", error.percent])
    })
    print(model.test.error)
    test.error.list[[paste(test.fold, model.name)]] <- model.test.error
  }
}
test.error <- do.call(rbind, test.error.list)

save(test.error, file="test.error.RData")
