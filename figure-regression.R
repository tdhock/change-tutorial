source("packages.R")

load("breakpoint.learning.RData")
## TODO: do this using Segmentor. data!
target.mat <- breakpoint.learning$targets$original
all.features.mat <- breakpoint.learning$features[rownames(target.mat), ]

train.dt <- data.table(
  profile.id=sub("[.].*", "", rownames(all.features.mat)),
  chromosome=sub(".*[.]", "", rownames(all.features.mat)),
  feature=all.features.mat[, "log2.n"],
  target.mat)
BIC.df <- data.frame(slope=1, intercept=0, model="BIC")
train.dt[, pred.log.lambda := feature ] #for the BIC
train.dt$model <- "BIC"
train.dt[, residual := targetIntervalResidual(cbind(min.L, max.L), pred.log.lambda)]
total.BIC <- train.dt[, list(
  total.residual=sum(residual),
  intervals=.N
  ), by=.(model, sign.residual=sign(residual))]
total.BIC
possible <- train.dt[, list(
  negative=sum(-Inf < min.L),
  positive=sum(max.L < Inf)
  )]

gg <- ggplot()+
  geom_text(aes(
    1.4, 6, color=model,
    label=sprintf(
      "total too high = %.1f (%d intervals / %d possible)",
      total.residual, intervals, possible$negative)),
            data=total.BIC[sign.residual==1, ],
            hjust=0)+
  geom_text(aes(
    1.4, 5.5, color=model,
    label=sprintf("%d intervals correctly predicted", intervals)),
            data=total.BIC[sign.residual==0, ],
            hjust=0)+
  geom_text(aes(
    1.4, 5, color=model,
    label=sprintf(
      "total too low = %.1f (%d intervals / %d possible)",
      total.residual, intervals, possible$positive)),
            data=total.BIC[sign.residual==-1, ],
            hjust=0)+
  geom_segment(aes(
    feature, pred.log.lambda, xend=feature, color=model,
    yend=ifelse(residual==0, NA, pred.log.lambda-residual)),
               size=1,
               linetype="dotted",
               data=train.dt)+
  geom_abline(aes(slope=slope, intercept=intercept, color=model), data=BIC.df, size=1)+
  geom_point(aes(feature, ifelse(is.finite(min.L), min.L, NA), fill="min"), data=train.dt, shape=21)+
  geom_point(aes(feature, ifelse(is.finite(max.L), max.L, NA), fill="max"), data=train.dt, shape=21)+
  scale_fill_manual("limit", values=c(min="black", max="white"))+
  scale_color_manual(values=c(BIC="red", learned="deepskyblue"))+
  scale_x_continuous("feature = log(log(n = number of data points to segment))")+
  scale_y_continuous("<-- more changes     log(penalty)     less changes -->")
print(gg)

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

gg <- ggplot()+
  geom_text(aes(
    1.4, 6, color=model,
    label=sprintf(
      "total too high = %.1f (%d intervals / %d possible)",
      total.residual, intervals, possible$negative)),
            data=total.BIC[sign.residual==1, ],
            hjust=0)+
  geom_text(aes(
    1.4, 5.5, color=model,
    label=sprintf("%d intervals correctly predicted", intervals)),
            data=total.BIC[sign.residual==0, ],
            hjust=0)+
  geom_text(aes(
    1.4, 5, color=model,
    label=sprintf(
      "total too low = %.1f (%d intervals / %d possible)",
      total.residual, intervals, possible$positive)),
            data=total.BIC[sign.residual==-1, ],
            hjust=0)+
  geom_segment(aes(
    feature, pred.log.lambda, xend=feature, color=model,
    yend=ifelse(residual==0, NA, pred.log.lambda-residual)),
               size=1,
               linetype="dotted",
               data=train.dt)+
  geom_abline(aes(slope=slope, intercept=intercept, color=model), data=BIC.df, size=1)+
  geom_point(aes(feature, ifelse(is.finite(min.L), min.L, NA), fill="min"), data=train.dt, shape=21)+
  geom_point(aes(feature, ifelse(is.finite(max.L), max.L, NA), fill="max"), data=train.dt, shape=21)+
  scale_fill_manual("limit", values=c(min="black", max="white"))+
  scale_color_manual(values=c(BIC="red", learned="deepskyblue"))+
  scale_x_continuous("feature = log(log(n = number of data points to segment))")+
  scale_y_continuous("<-- more changes     log(penalty)     less changes -->")
print(gg)

## Bug?
## fit <- survreg(Surv(min.L, max.L, type="interval2") ~ feature, train.dt, dist="gaussian")
## Error in coxph.wtest(t(x) %*% (wt * x), c((wt * eta + weights * deriv$dg) %*%  : 
##   NA/NaN/Inf in foreign function call (arg 3)

fit <- with(train.dt, {
  IntervalRegressionUnregularized(
    cbind(feature), cbind(min.L, max.L))
})
pred.dt <- data.table(train.dt)
pred.dt[, pred.log.lambda := fit$predict(cbind(feature))]
pred.dt[, residual := targetIntervalResidual(cbind(min.L, max.L), pred.log.lambda)]
pred.dt$model <- "learned"
total.learned <- pred.dt[, list(
  total.residual=sum(residual),
  intervals=.N
  ), by=.(model, sign.residual=sign(residual))]
total.learned

gg+
  geom_text(aes(
    1.8, -4, color=model,
    label=sprintf(
      "total too high = %.1f (%d intervals / %d possible)",
      total.residual, intervals, possible$negative)),
            data=total.learned[sign.residual==1, ],
            hjust=0)+
  geom_text(aes(
    1.8, -4.5, color=model,
    label=sprintf("%d intervals correctly predicted", intervals)),
            data=total.learned[sign.residual==0, ],
            hjust=0)+
  geom_text(aes(
    1.8, -5, color=model,
    label=sprintf(
      "total too low = %.1f (%d intervals / %d possible)",
      total.residual, intervals, possible$positive)),
            data=total.learned[sign.residual==-1, ],
            hjust=0)+
  geom_segment(aes(
    feature, pred.log.lambda, xend=feature, color=model,
    yend=ifelse(residual==0, NA, pred.log.lambda-residual)),
               linetype="solid",
               data=pred.dt)+
  geom_line(aes(
    feature, pred.log.lambda, color=model), data=pred.dt)

