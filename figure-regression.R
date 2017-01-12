source("packages.R")

load("breakpoint.learning.RData")

target.mat <- breakpoint.learning$targets$original
all.features.mat <- breakpoint.learning$features[rownames(target.mat), ]

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
BIC.df <- data.frame(slope=1, intercept=0, model="BIC")
train.df <- data.frame(
  feature=all.features.mat[, "log2.n"],
  target.mat)

gg <- ggplot()+
  geom_point(aes(feature, ifelse(is.finite(min.L), min.L, NA), fill="min"), data=train.df, shape=21)+
  geom_point(aes(feature, ifelse(is.finite(max.L), max.L, NA), fill="max"), data=train.df, shape=21)+
  scale_fill_manual("limit", values=c(min="black", max="white"))+
  geom_abline(aes(slope=slope, intercept=intercept, color=model), data=BIC.df, size=2)+
  scale_color_manual(values=c(BIC="red", learned="deepskyblue"))+
  scale_x_continuous("feature = log(log(n = number of data points to segment))")+
  scale_y_continuous("<-- more changes     log(penalty)     less changes -->")
print(gg)

## Bug?
fit <- survreg(Surv(min.L, max.L, type="interval2") ~ feature, train.df, dist="gaussian")
## Error in coxph.wtest(t(x) %*% (wt * x), c((wt * eta + weights * deriv$dg) %*%  : 
##   NA/NaN/Inf in foreign function call (arg 3)

intercept.vec <- seq(-3, 1, by=0.1)
for(pid.chr in rownames(target.mat)){
  breakpoint.learning$modelSelection$original[[pid.chr]]
}

