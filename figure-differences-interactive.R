library(animint2)
library(data.table)
data(neuroblastoma, package="neuroblastoma")

nb.dts <- lapply(neuroblastoma, data.table, key=c("profile.id", "chromosome"))

if(file.exists("Segmentor.models.RData")){
  load("Segmentor.models.RData")
}else{
  segs.list <- list()
  loss.list <- list()
  for(problem.i in seq_along(labeled.problem.names)){
    problem.name <- labeled.problem.names[[problem.i]]
    meta <- meta.df[problem.i,]
    cat(sprintf("%4d / %4d problems %s\n", problem.i, length(labeled.problem.names), problem.name))
    pro <- problem.list[[problem.name]]
    fit <- Segmentor3IsBack::Segmentor(pro$logratio, model=2, Kmax=max.segments)
    rss.vec <- rep(NA, max.segments)
    for(n.segments in 1:max.segments){
      end <- fit@breaks[n.segments, 1:n.segments]
      data.before.change <- end[-n.segments]
      data.after.change <- data.before.change+1
      pos.before.change <- as.integer(
      (pro$position[data.before.change]+pro$position[data.after.change])/2)
      start <- c(1, data.after.change)
      chromStart <- c(pro$position[1], pos.before.change)
      chromEnd <- c(pos.before.change, max(pro$position))
      seg.mean.vec <- fit@parameters[n.segments, 1:n.segments]
      data.mean.vec <- rep(seg.mean.vec, end-start+1)
      residual.vec <- pro$logratio - data.mean.vec
      rss.vec[n.segments] <- sum(residual.vec * residual.vec)
      segs.list[[paste(problem.name, n.segments)]] <- data.table(
        meta,
        n.segments,
        start,
        end,
        chromStart,
        chromEnd,
        mean=seg.mean.vec)
    }
    if(FALSE){
      ## The likelihood computed by the Segmentor function is just an
      ## affine transformation of the sum of squared residuals.
      plot(rss.vec, fit@likelihood)
    }
    loss.list[[paste(problem.name, n.segments)]] <- data.table(
      meta,
      n.segments=1:max.segments,
      loss=rss.vec)
  }
  Segmentor.models <- list(
    segs=do.call(rbind, segs.list),
    loss=do.call(rbind, loss.list))
  save(Segmentor.models, file="Segmentor.models.RData")
}

selection <- Segmentor.models$loss[, penaltyLearning::modelSelection(
  .SD, complexity="n.segments"
), by=.(profile.id, chromosome)]
changes <- Segmentor.models$segs[1 < start]
errors <- penaltyLearning::labelError(
  selection, neuroblastoma$annotations, changes,
  change.var="chromStart",
  label.vars=c("min", "max"),
  model.vars="n.segments",
  problem.vars=c("profile.id", "chromosome"))
target.dt <- penaltyLearning::targetIntervals(
  errors$model.errors, c("profile.id", "chromosome"))
labeled.profiles <- nb.dts$profiles[target.dt]
feature.dt <-  labeled.profiles[, list(
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
train.dt[, residual := penaltyLearning::targetIntervalResidual(
  cbind(min.log.lambda, max.log.lambda), pred.log.lambda)]

possible <- train.dt[, list(
  negative=sum(-Inf < min.log.lambda),
  positive=sum(max.log.lambda < Inf),
  total=.N
  )]

fit <- survival::survreg(
  survival::Surv(min.log.lambda, max.log.lambda, type="interval2") ~ feature,
  train.dt, dist="gaussian")
pred.dt <- data.table(train.dt)[
, pred.log.lambda := predict(fit)
][
, residual := penaltyLearning::targetIntervalResidual(
  cbind(min.log.lambda, max.log.lambda), pred.log.lambda)
][
, model.name := "learned"
]

model.list <- list(
  ##learned=pred.dt,
  AIC=data.table(target.dt, pred.log.lambda=log(2)),
  BIC=train.dt)
roc.list <- list()
auc.list <- list()
both.pred.list <- list()
pred.dot.list <- list()
auc.polygon.list <- list()
for(model.name in names(model.list)){
  model.dt <- model.list[[model.name]]
  both.pred.list[[model.name]] <- data.table(model.name, model.dt[, .(
    profile.id, chromosome, pred.log.lambda)])
  feature.result <- penaltyLearning::ROChange(
    errors$model.errors, model.dt, c("profile.id", "chromosome"))
  pred.dot.list[[model.name]] <- data.table(
    model.name, feature.result$thresholds[threshold=="predicted",])
  roc.list[[model.name]] <- data.table(model.name, feature.result$roc)
  auc.list[[model.name]] <- data.table(model.name, auc=feature.result$auc)
  auc.polygon.list[[model.name]] <- data.table(model.name, feature.result$auc.polygon)
}
both.pred <- rbindlist(both.pred.list)
roc <- do.call(rbind, roc.list)
auc <- do.call(rbind, auc.list)
pred.dot <- do.call(rbind, pred.dot.list)
auc.polygon <- do.call(rbind, auc.polygon.list)

## Zoomed ROC curves.
auc$TPR <- c(0.85, 1)
auc$FPR <- c(0.2, 0.05)
ggplot()+
  theme_bw()+
  geom_text(aes(
    FPR, TPR, color=model.name, label=sprintf("AUC = %.4f", auc)),
    vjust=0,
    data=auc)+
  geom_path(aes(
    FPR, TPR, color=model.name, group=model.name),
    data=roc)+
  scale_fill_manual(values=c(default="grey", min.error="black"))+
  geom_point(aes(
    FPR, TPR, color=model.name, fill=threshold),
    data=pred.dot,
    size=3,
    shape=21)+
  coord_equal(xlim=c(0, 0.3), ylim=c(0.4, 1))

roc[, mid.thresh := (min.thresh+max.thresh)/2]
rate.grid.list <- list(
  TPR=seq(0.4, 0.95, by=0.05),
  FPR=c(0.015, 0.02,0.03,0.04,0.061,0.08,0.1))
roc.diff.dt.list <- list()
roc.grid.dt.list <- list()
for(rate in names(rate.grid.list)){
  grid.vec <- rate.grid.list[[rate]]
  grid.dt <- data.table(grid=sort(c(pred.dot[[rate]], grid.vec)))
  set(roc, j="grid", value=roc[[rate]])
  grid.roc <- roc[
  , .SD[grid.dt, roll=Inf, mult="first", on="grid"]
  , by=model.name
  ]#max FP for given TP.
  rate.wide <- dcast(
    grid.roc,
    grid ~ model.name,
    value.var=c("FPR", "TPR","fp","tp","mid.thresh"))
  roc.grid.dt.list[[rate]] <- data.table(rate, grid.roc)
  roc.diff.dt.list[[rate]] <- data.table(rate, rate.wide)
}
roc.grid.dt <- rbindlist(roc.grid.dt.list)
roc.diff.dt <- rbindlist(roc.diff.dt.list)[, `:=`(
  tp_diff=tp_BIC-tp_AIC,
  fp_diff=fp_AIC-fp_BIC
)][, `:=`(
  diff=ifelse(rate=="FPR", tp_diff, fp_diff),
  vjust=ifelse(rate=="FPR", 0, 1),
  hjust=ifelse(rate=="FPR", 1, 0)
)]
err.dt <- data.table(
  errors$label.errors,
  key=c("min.log.lambda","max.log.lambda")
)[, tp := possible.fn-fn]

adj.pred <- roc.grid.dt[, {
  data.table(
    model.name, mid.thresh, key="model.name"
  )[both.pred, .(
    model.name, profile.id, chromosome,
    adj.log.lambda=pred.log.lambda+mid.thresh
  ), on="model.name"]
}, by=.(rate,grid)]
grid.labels <- err.dt[adj.pred, on=.(
  profile.id, chromosome,
  max.log.lambda>adj.log.lambda,
  min.log.lambda<adj.log.lambda
)][
, `:=`(same=all(status==status[1]), sameN=.N), by=.(rate,grid,profile.id,chromosome)
]
grid.labels[, .(# to check, ok.
  fp=sum(fp),
  tp=sum(tp)
), by=.(rate,grid,model.name)]
not.same <- grid.labels[same==FALSE]
not.same.wide <- dcast(
  not.same,
  profile.id+chromosome+rate+grid ~ model.name,
  value.var=c("tp","fp")
)[, `:=`(
  tp_diff = tp_BIC-tp_AIC,
  fp_diff = fp_AIC-fp_BIC
)]
not.same.wide[, .(
  tp_diff=sum(tp_diff),
  fp_diff=sum(fp_diff)
), keyby=.(rate,grid)]
not.same[, .(count=.N), by=.(rate,grid,model.name,status)]
ggplot()+
  theme_bw()+
  geom_text(aes(
    FPR, TPR, color=model.name, label=sprintf("AUC = %.4f", auc)),
    vjust=0,
    data=auc)+
  geom_segment(aes(
    FPR_AIC, TPR_AIC,
    xend=FPR_BIC, yend=TPR_BIC),
    color="grey50",
    data=roc.diff.dt)+
  geom_text(aes(
    FPR_AIC, TPR_BIC, label=diff, hjust=hjust, vjust=vjust),
    data=roc.diff.dt)+
  geom_path(aes(
    FPR, TPR, color=model.name, group=model.name),
    data=roc)+
  scale_fill_manual(values=c(default="grey", min.error="black"))+
  geom_point(aes(
    FPR, TPR, color=model.name, fill=threshold),
    data=pred.dot,
    size=3,
    shape=21)+
  coord_equal(xlim=c(0, 0.3), ylim=c(0.4, 1))

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
  geom_text(aes(
    FPR, TPR, label=sprintf("AUC = %.4f", auc)),
    vjust=0,
    data=auc)+
  geom_path(aes(
    FPR, TPR, 
    color=dist.from.zero, group=model.name),
    data=roc)+
  scale_fill_manual(values=c(default="grey", min.error="black"))+
  geom_point(aes(
    FPR, TPR, fill=threshold),
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
viz <- list(
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
  source="https://github.com/tdhock/change-tutorial/blob/master/figure-differences-interactive.R",
  first=list())
pred.thresh.only <- some.thresh[threshold=="predicted"]
for(row.i in 1:nrow(pred.thresh.only)){
  r <- pred.thresh.only[row.i, ]
  viz$first[[paste0(r$model.name, ".thresh")]] <- r$mid.thresh
}
animint2dir(viz, "figure-differences-interactive")
if(FALSE){
  animint2pages(viz, "2023-11-interval-regression-differences")
}

## TODO binseg vs DP.

## TODO AIC, CV.
