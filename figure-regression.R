source("packages.R")

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

BIC.df <- data.frame(slope=1, intercept=0, model="BIC")
train.dt[, pred.log.lambda := feature ] #for the BIC
train.dt$model <- "BIC"
train.dt[, residual := targetIntervalResidual(cbind(min.log.lambda, max.log.lambda), pred.log.lambda)]
possible <- train.dt[, list(
  negative=sum(-Inf < min.log.lambda),
  positive=sum(max.log.lambda < Inf),
  total=.N
)]
totalResiduals <- function(dt, feature, log.lambda){
  res <- dt[, list(
    total.residual=sum(residual),
    intervals=.N
  ), by=.(model, sign.residual=sign(residual))][order(-sign.residual),]
  res[, feature :=  feature ]
  res[, log.lambda :=  log.lambda]
  res[, possible :=  possible[, c(positive, total, negative)]]
  res[, rate := c("FNR", NA, "FPR")]
  res[, too := c("high", NA, "low")]
  res
}
total.BIC <- totalResiduals(train.dt, 1.4, c(6, 5.5, 5))

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

## Note that this FPR/FNR is defined using the target intervals, so is
## not exactly the same as the numbers that would result from the
## computation using the labels. (but it is pretty much the same).

gg.BIC <- ggplot()+
  geom_text(aes(
    feature, log.lambda, color=model,
    label=ifelse(
      sign.residual==0, 
      sprintf(
        "%d / %d = %.1f%% intervals correctly predicted",
        intervals, possible, 100*intervals/possible),      
      sprintf(
        "total too %s = %.1f (%d intervals / %d possible = %.1f%% %s)",
        too, total.residual, intervals, possible,
        100*intervals/possible, rate))),
    data=total.BIC,
    hjust=0)+
  geom_segment(aes(
    feature, pred.log.lambda, xend=feature, color=model,
    yend=ifelse(residual==0, NA, pred.log.lambda-residual)),
               size=1,
               linetype="dotted",
               data=train.dt)+
  geom_abline(aes(slope=slope, intercept=intercept, color=model), data=BIC.df, size=1)+
  geom_point(aes(feature, ifelse(is.finite(min.log.lambda), min.log.lambda, NA), fill="min"), data=train.dt, shape=21)+
  geom_point(aes(feature, ifelse(is.finite(max.log.lambda), max.log.lambda, NA), fill="max"), data=train.dt, shape=21)+
  scale_fill_manual("limit", values=c(min="black", max="white"))+
  scale_color_manual(values=c(BIC="red", learned="deepskyblue"))+
  scale_x_continuous("feature = log(log(n = number of data points to segment))")+
  scale_y_continuous("<-- more changes     log(penalty)     less changes -->")
print(gg.BIC)

## Bug?
## fit <- survreg(Surv(min.log.lambda, max.log.lambda, type="interval2") ~ feature, train.dt, dist="gaussian")
## Error in coxph.wtest(t(x) %*% (wt * x), c((wt * eta + weights * deriv$dg) %*%  : 
##   NA/NaN/Inf in foreign function call (arg 3)

fit <- with(train.dt, {
  IntervalRegressionUnregularized(
    cbind(feature), cbind(min.log.lambda, max.log.lambda))
})
pred.dt <- data.table(train.dt)
pred.dt[, pred.log.lambda := fit$predict(cbind(feature))]
pred.dt[, residual := targetIntervalResidual(cbind(min.log.lambda, max.log.lambda), pred.log.lambda)]
pred.dt$model <- "learned"
total.learned <- totalResiduals(pred.dt, 1.8, c(-4, -4.5, -5))
total.learned

gg.compare <- gg.BIC+
  geom_text(aes(
    feature, log.lambda, color=model,
    label=ifelse(
      sign.residual==0, 
      sprintf(
        "%d / %d = %.1f%% intervals correctly predicted",
        intervals, possible, 100*intervals/possible),      
      sprintf(
        "total too %s = %.1f (%d intervals / %d possible = %.1f%% %s)",
        too, total.residual, intervals, possible,
        100*intervals/possible, rate))),
    data=total.learned,
    hjust=0)+
  geom_segment(aes(
    feature, pred.log.lambda, xend=feature, color=model,
    yend=ifelse(residual==0, NA, pred.log.lambda-residual)),
               linetype="solid",
               data=pred.dt)+
  geom_line(aes(
    feature, pred.log.lambda, color=model), data=pred.dt)
print(gg.compare)

pdf("figure-regression.pdf", 11, 7)
print(gg.compare)
dev.off()
