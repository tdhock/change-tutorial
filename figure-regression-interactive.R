source("packages.R")
source("animint.R")
## Visualizes all of the data points, comparing two models (BIC and
## learned), using three plots: threshold, ROC, regression. In
## contrast figure-regression-interactive-some.R shows the same two
## models for only a few data points, but five plots total. The two
## extra plots show the data and the model selection/error functions.
data(neuroblastoma, package="neuroblastoma")
load("Segmentor.models.RData")
selection <-
  Segmentor.models$loss[, modelSelection(
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
some.profiles <- data.table(
  neuroblastoma$profiles)[target.dt, on=.(profile.id, chromosome)]
feature.dt <-  some.profiles[, list(
  log2.n=log(log(.N)),
  log.mad=log(median(abs(diff(logratio))))
), by=.(profile.id, chromosome)]
train.dt <- data.table(
  feature=feature.dt$log2.n,
  target.dt)

## The log2.n feature is log(log(n = number of data points to segment)).

## The model selection criterion is argmin_K totalSquareLoss_K +
## lambda * K where K is the number of segments in the piecewise
## constant model, and lambda is a non-negative penalty constant.

## The target.mat contains the interval of log(lambda) values for
## which the number of incorrect labels is minimized.

## The BIC model selection criterion is lambda = log(n), where n is
## the number of data points to segment. This implies log(lambda) =
## log(log(n)) = the log2.n feature in all.features.mat.

## We can thus visualize the BIC penalty as a line with slope 1 and
## intercept 0 in a plot of log(lambda) versus log(log(n)).
BIC.df <- data.frame(slope=1, intercept=0, model.name="BIC")
train.dt[, pred.log.lambda := feature ] #for the BIC
train.dt$model.name <- "BIC"
train.dt[, residual := targetIntervalResidual(
  cbind(min.log.lambda, max.log.lambda), pred.log.lambda)]

possible <- train.dt[, list(
  negative=sum(-Inf < min.log.lambda),
  positive=sum(max.log.lambda < Inf),
  total=.N
  )]

fit <- survreg(
  Surv(min.log.lambda, max.log.lambda, type="interval2") ~ feature,
  train.dt, dist="gaussian")
pred.dt <- data.table(train.dt)
pred.dt[, pred.log.lambda := predict(fit)]
pred.dt[, residual := targetIntervalResidual(
  cbind(min.log.lambda, max.log.lambda), pred.log.lambda)]
pred.dt$model.name <- "learned"

model.list <- list(
  BIC=train.dt,
  learned=pred.dt)
roc.list <- list()
auc.list <- list()
pred.dot.list <- list()
auc.polygon.list <- list()
for(model.name in names(model.list)){
  model.dt <- model.list[[model.name]]
  feature.result <- ROChange(errors$model.errors, model.dt, c("profile.id", "chromosome"))
  pred.dot.list[[model.name]] <- data.table(
    model.name, feature.result$thresholds[threshold=="predicted",])
  roc.list[[model.name]] <- data.table(model.name, feature.result$roc)
  auc.list[[model.name]] <- data.table(model.name, auc=feature.result$auc)
  auc.polygon.list[[model.name]] <- data.table(model.name, feature.result$auc.polygon)
}
roc <- do.call(rbind, roc.list)
auc <- do.call(rbind, auc.list)
pred.dot <- do.call(rbind, pred.dot.list)
auc.polygon <- do.call(rbind, auc.polygon.list)

## Zoomed ROC curves.
auc$TPR <- c(0.85, 0.9)
auc$FPR <- 1
ggplot()+
  theme_bw()+
  geom_text(aes(FPR, TPR, color=model.name, label=sprintf("AUC = %.4f", auc)),
            vjust=0,
            data=auc)+
  geom_path(aes(FPR, TPR, color=model.name, group=model.name),
            data=roc)+
  scale_fill_manual(values=c(default="grey", min.error="black"))+
  geom_point(aes(FPR, TPR, color=model.name, fill=threshold),
             data=pred.dot,
             size=3,
             shape=21)+
  coord_equal(xlim=c(0, 0.25), ylim=c(0.75, 1))

roc[, mid.thresh := (min.thresh+max.thresh)/2]
pred.dot[, mid.thresh := (min.thresh+max.thresh)/2]
roc[, threshold := "other"]
## We would like to sample some thresholds evenly in the ROC space, so
## here we compute dist.from.zero, which is the cumulative euclidean
## distance from (0,0) in the ROC space.
roc[, dist.from.zero := cumsum(sqrt({
  diff(c(0, FPR))^2 + diff(c(0, TPR))^2
})), by=model.name]
grid.thresh <- roc[is.finite(mid.thresh), {
  .SD[as.integer(seq(1, .N, l=50)[c(1,5,10,15,20, 25, 30:50)]), ]
}, by=model.name]
grid.thresh <- roc[is.finite(mid.thresh), {
  r <- range(dist.from.zero)
  d.vec <- seq(r[1], r[2], l=50)
  i.vec <- sapply(d.vec, function(d)which.min(abs(dist.from.zero-d)))
  .SD[i.vec,]
}, by=model.name]
some.thresh <- unique(rbind(
  grid.thresh[, names(pred.dot), with=FALSE],
  pred.dot))
setkey(some.thresh, model.name, mid.thresh)
some.thresh[, prev.thresh := {
  d <- diff(mid.thresh)/2
  c(mid.thresh[1]-d[1], mid.thresh[-.N]+d)
  }, by=model.name]
some.thresh[, next.thresh := {
  d <- diff(mid.thresh)/2
  c(mid.thresh[-.N]+d, mid.thresh[.N]+d[.N-1])
}, by=model.name]
some.thresh[, stopifnot(
  identical(prev.thresh[-1], next.thresh[-.N])
  ), by=model.name]
some.thresh[, slope := ifelse(model.name=="BIC", 1, coef(fit)[["feature"]])]
some.thresh[, intercept := {
  ifelse(model.name=="BIC", 0, coef(fit)[["(Intercept)"]])+mid.thresh
}]
ggplot()+
  theme_bw()+
  geom_text(aes(FPR, TPR, label=sprintf("AUC = %.4f", auc)),
            vjust=0,
            data=auc)+
  geom_path(aes(FPR, TPR, 
                color=dist.from.zero, group=model.name),
            data=roc)+
  scale_fill_manual(values=c(default="grey", min.error="black"))+
  geom_point(aes(FPR, TPR, fill=threshold),
             data=pred.dot,
             size=3,
             shape=21)+
  scale_color_gradient(low="black", high="red")

both.pred <- rbind(train.dt, pred.dt)
setkey(both.pred, model.name)
setkey(some.thresh, model.name)
some.thresh.pred <- both.pred[some.thresh, allow.cartesian=TRUE]
some.thresh.pred[, pred.log.lambda := slope * feature + intercept]
some.thresh.pred[, residual := targetIntervalResidual(
  cbind(min.log.lambda, max.log.lambda), pred.log.lambda)]
total.thresh.pred <- unique(some.thresh.pred)[, list(
  total.residual=sum(residual),
  intervals=.N
), by=.(model.name, sign.residual=sign(residual), mid.thresh)]
total.thresh.pred[, top := ifelse(model.name=="BIC", 6.5, -4)]
total.thresh.pred[, left := ifelse(model.name=="BIC", 1.45, 1.8)]
total.thresh.pred[, log.penalty := top-(1-sign.residual)*0.7]
residual.thresh.pred <- some.thresh.pred[residual!=0,]

add_selector <- function(DT){
  DT[, selector := paste0(model.name, ".thresh")]
}
add_selector(some.thresh)
add_selector(pred.dot)
add_selector(total.thresh.pred)
add_selector(residual.thresh.pred)
thresh.colors <- c(predicted="green", min.error="black", other="white")
viz <- animint(
  title="BIC versus learned penalty in neuroblastoma data",
  thresholds=ggplot()+
    theme_bw()+
    guides(fill="none", color="none")+
    geom_line(aes(
      mid.thresh, errors,
      group=model.name,
      color=model.name),
      showSelected="model.name",
      data=some.thresh)+
    geom_tallrect(aes(
      xmin=prev.thresh, xmax=next.thresh,
      color=model.name),
      clickSelects=c(selector="mid.thresh"),
      showSelected="model.name",
      size=3,
      alpha=0.5,
      data=some.thresh)+
    ylab("total incorrect labels")+
    xlab("constant/threshold added to log(penalty)")+
    scale_fill_manual(values=thresh.colors)+
    scale_color_manual(values=c(BIC="red", learned="deepskyblue"))+
    geom_point(aes(
    (max.thresh+min.thresh)/2, errors,
    color=model.name,
    fill=threshold),
    clickSelects=c(selector="mid.thresh"),
    showSelected=c("model.name","threshold"),
    data=pred.dot,
    size=4,
    shape=21),
  roc=ggplot()+
    theme_bw()+
    geom_text(aes(
      FPR, TPR,
      color=model.name,
      label=sprintf("%s AUC = %.4f", model.name, auc)),
      hjust=1,
      data=auc)+
    geom_path(aes(
      FPR, TPR,
      color=model.name,
      group=model.name),
      size=2,
      data=roc)+
    scale_fill_manual(values=thresh.colors, breaks=names(thresh.colors))+
    scale_color_manual(values=c(BIC="red", learned="deepskyblue"))+
    geom_point(aes(
      FPR, TPR),
      showSelected=c(selector="mid.thresh","threshold","model.name"),
      data=some.thresh,
      color="grey",
      size=6,
      shape=21)+
    geom_point(aes(
      FPR, TPR,
      color=model.name, fill=threshold),
      clickSelects=c(selector="mid.thresh"),
      data=some.thresh,
      size=4,
      alpha=0.8,
      shape=21)+
    coord_equal(),
  regression=ggplot()+
    theme_bw()+
    theme_animint(width=800)+
    geom_text(aes(
      left, log.penalty, 
      color=model.name,
      key=paste(model.name, sign.residual),
      label=ifelse(sign.residual==0, sprintf(
        "%d intervals correctly predicted", intervals), sprintf(
          "total too %s = %.1f (%d/%d intervals)",
          ifelse(sign.residual==1, "high", "low"),
          total.residual, intervals,
          possible[, ifelse(sign.residual==1, positive, negative)]))),
      showSelected=c(selector="mid.thresh"),
      data=total.thresh.pred,
      hjust=0)+
    geom_point(aes(
      feature, min.log.lambda, fill=limit),
      data=data.table(
        limit="min", train.dt[is.finite(min.log.lambda),]),
      shape=21)+
    geom_point(aes(
      feature, max.log.lambda, fill=limit),
      data=data.table(
        limit="max", train.dt[is.finite(max.log.lambda),]),
      shape=21)+
    geom_abline(aes(
      slope=slope,
      intercept=intercept, color=model.name),
      showSelected=c(selector="mid.thresh"),
      data=some.thresh)+
    geom_segment(aes(
      feature, pred.log.lambda, xend=feature, color=model.name,
      yend=pred.log.lambda-residual),
      showSelected=c(selector="mid.thresh"),
      size=1,
      data=residual.thresh.pred)+
    scale_fill_manual(values=c(min="black", max="white"))+
    scale_color_manual(values=c(BIC="red", learned="deepskyblue"))+
    scale_x_continuous("feature = log(log(n = number of data points to segment))")+
    scale_y_continuous("<-- more changes     log(penalty)     less changes -->"),
  source="https://github.com/tdhock/change-tutorial/blob/master/figure-regression-interactive.R",
  out.dir="figure-regression-interactive",
  first=list())
pred.thresh.only <- some.thresh[threshold=="predicted"]
for(row.i in 1:nrow(pred.thresh.only)){
  r <- pred.thresh.only[row.i, ]
  viz$first[[paste0(r$model.name, ".thresh")]] <- r$mid.thresh
}
viz
if(FALSE){
  animint2pages(viz, "2023-08-interval-regression-BIC-vs-learned")
  ## https://rcdata.nau.edu/genomic-ml/animint-gallery/2017-05-08-BIC-versus-learned-penalty/index.html
}
