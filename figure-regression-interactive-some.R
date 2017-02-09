source("packages.R")

load("breakpoint.learning.RData")
load("Segmentor.models.RData")
ids.str <- paste(c(1, 4, 6, 8, 10, 11))
all.target.mat <- breakpoint.learning$targets$original
target.mat <- all.target.mat[sub("[.].*", "", rownames(all.target.mat)) %in% ids.str, ]
all.features.mat <- breakpoint.learning$features[rownames(target.mat), ]


data(neuroblastoma, package="neuroblastoma")
selection <- Segmentor.models$loss[, modelSelection(
  .SD, complexity="n.segments"
  ), by=.(profile.id, chromosome)]
changes <- Segmentor.models$segs[1 < start, ]
some.anns <- data.table(neuroblastoma$annotations)[profile.id %in% ids.str,]
errors <- labelError(
  selection, some.anns, changes,
  change.var="chromStart",
  label.vars=c("min", "max"),
  model.vars="n.segments",
  problem.vars=c("profile.id", "chromosome"))

model.fun.list <- list(
  BIC=function(train.dt){
    data.table(slope=1, intercept=0)
  }, learned=function(train.dt){
    target.mat <- train.dt[, cbind(min.L, max.L)]
    feature.mat <- train.dt[, cbind(feature)]
    fit <- IntervalRegressionUnregularized(feature.mat, target.mat)
    slope <- with(fit, param.mat["feature",]/sd.vec)
    data.table(
      slope,
      intercept=fit$param.mat["(Intercept)",]-slope*fit$mean.vec
      )
  })
model.label.vec <- c(
  BIC="BIC penalty, feature=log(log(n))",
  learned="Learned penalty, feature=log(var.est)")
feature.name.vec <- c(
  BIC="log2.n",
  learned="log.mad")
train.dt.list <- list()
totals.list <- list()
reg.line.list <- list()
for(model.name in names(feature.name.vec)){
  model.label <- model.label.vec[[model.name]]
  feature.name <- feature.name.vec[[model.name]]
  model.train.dt <- data.table(
    model.name, model.label,
    profile.id=sub("[.].*", "", rownames(all.features.mat)),
    chromosome=sub(".*[.]", "", rownames(all.features.mat)),
    feature=all.features.mat[, feature.name],
    target.mat)
  model.fun <- model.fun.list[[model.name]]
  reg.dt <- data.table(model.fun(model.train.dt), model.name, model.label)
  reg.line.list[[model.name]] <- reg.dt
  model.train.dt[, pred.log.lambda := reg.dt[, feature*slope+intercept] ]
  model.train.dt[, residual := {
    targetIntervalResidual(cbind(min.L, max.L), pred.log.lambda)
  }]
  possible <- model.train.dt[, list(
    positive=sum(-Inf < min.L),
    negative=sum(max.L < Inf)
    )]
  min.feature <- min(model.train.dt$feature)
  totals.list[[model.name]] <- model.train.dt[, list(
    feature=min.feature,
    total.residual=sum(residual),
    intervals=.N
    ), by=.(model.name, model.label, sign.residual=sign(residual))]
  train.dt.list[[model.name]] <- model.train.dt
}
train.dt <- do.call(rbind, train.dt.list)
totals <- do.call(rbind, totals.list)
reg.line <- do.call(rbind, reg.line.list)
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ model.label, scales="free")+
  geom_text(aes(
    feature, 4, 
    label=sprintf(
      "total too high = %.1f (%d intervals / %d possible)",
      total.residual, intervals, possible$negative)),
            data=totals[sign.residual==1, ],
            hjust=0)+
  geom_text(aes(
    feature, 3.5, 
    label=sprintf("%d intervals correctly predicted", intervals)),
            data=totals[sign.residual==0, ],
            hjust=0)+
  geom_segment(aes(
    feature, pred.log.lambda, xend=feature, 
    yend=ifelse(residual==0, NA, pred.log.lambda-residual)),
               size=1,
               linetype="dotted",
               data=train.dt)+
  geom_abline(aes(slope=slope, intercept=intercept), data=reg.line, size=1)+
  geom_point(aes(feature, ifelse(is.finite(min.L), min.L, NA), fill="min"), data=train.dt, shape=21)+
  geom_point(aes(feature, ifelse(is.finite(max.L), max.L, NA), fill="max"), data=train.dt, shape=21)+
  scale_fill_manual("limit", values=c(min="black", max="white"))+
  scale_x_continuous("feature")+
  scale_y_continuous("<-- more changes     log(penalty)     less changes -->")
print(gg)



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
  pred.dot.list[[model.name]] <- data.table(model.name, feature.result$thresholds)
  roc.list[[model.name]] <- data.table(model.name, feature.result$roc)
  auc.list[[model.name]] <- data.table(model.name, auc=feature.result$auc)
  auc.polygon.list[[model.name]] <- data.table(model.name, feature.result$auc.polygon)
}
roc <- do.call(rbind, roc.list)
auc <- do.call(rbind, auc.list)
pred.dot <- do.call(rbind, pred.dot.list)
auc.polygon <- do.call(rbind, auc.polygon.list)

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
some.thresh[, slope := ifelse(model.name=="BIC", 1, learned.slope)]
some.thresh[, intercept := {
  ifelse(model.name=="BIC", 0, learned.intercept)+mid.thresh
}]

both.pred <- rbind(train.dt, pred.dt)
setkey(both.pred, model.name)
setkey(some.thresh, model.name)
some.thresh.pred <- both.pred[some.thresh, allow.cartesian=TRUE]
some.thresh.pred[, pred.log.lambda := slope * feature + intercept]
some.thresh.pred[, residual := targetIntervalResidual(cbind(min.L, max.L), pred.log.lambda)]
total.thresh.pred <- some.thresh.pred[, list(
  total.residual=sum(residual),
  intervals=.N
  ), by=.(model.name, sign.residual=sign(residual), mid.thresh)]
total.thresh.pred
residual.thresh.pred <- some.thresh.pred[residual!=0,]

auc$TPR <- c(0.85, 0.9)
auc$FPR <- 0.2
bic.left <- 1.45
bic.top <- 6.5
bic.space <- 0.7
learned.left <- 1.8
learned.top <- -4
thresh.colors <- c(predicted="green", min.error="black", other="white")
viz <- list(
  title="BIC versus learned penalty in neuroblastoma data",
  thresholds=ggplot()+
    theme_bw()+
    guides(fill="none", color="none")+
    geom_line(aes(mid.thresh, errors, group=model.name, color=model.name,
                  showSelected=model.name), data=some.thresh)+
    geom_tallrect(aes(xmin=prev.thresh, xmax=next.thresh,
                      clickSelects.variable=paste0(model.name, ".thresh"),
                      clickSelects.value=mid.thresh,
                      color=model.name,
                      showSelected=model.name),
                  size=3,
                  alpha=0.5,
                  data=some.thresh)+
    ylab("total incorrect labels")+
    xlab("constant/threshold added to log(penalty)")+
    scale_fill_manual(values=thresh.colors)+
    scale_color_manual(values=c(BIC="red", learned="deepskyblue"))+
    geom_point(aes((max.thresh+min.thresh)/2, errors,
                   clickSelects.variable=paste0(model.name, ".thresh"),
                   clickSelects.value=mid.thresh,
                   color=model.name,
                   showSelected=model.name,
                   showSelected2=threshold,
                   fill=threshold),
               data=pred.dot,
               size=4,
               shape=21),
  roc=ggplot()+
    theme_bw()+
    geom_text(aes(FPR, TPR, color=model.name, label=sprintf("AUC = %.4f", auc)),
              hjust=0,
              data=auc)+
    geom_path(aes(FPR, TPR, color=model.name, group=model.name),
              size=2,
              data=roc)+
    scale_fill_manual(values=thresh.colors, breaks=names(thresh.colors))+
    scale_color_manual(values=c(BIC="red", learned="deepskyblue"))+
    geom_point(aes(FPR, TPR,
                   showSelected.variable=paste0(model.name, ".thresh"),
                   showSelected.value=mid.thresh,
                   showSelected2=threshold,
                   showSelected=model.name),
               data=some.thresh,
               color="grey",
               size=6,
               shape=21)+
    geom_point(aes(FPR, TPR,
                   clickSelects.variable=paste0(model.name, ".thresh"),
                   clickSelects.value=mid.thresh,
                   color=model.name, fill=threshold),
               data=some.thresh,
               size=4,
               alpha=0.8,
               shape=21)+
    coord_equal(),
  regression=ggplot()+
    theme_bw()+
    theme_animint(width=800)+
    geom_text(aes(
      bic.left, bic.top, color=model.name,
      showSelected.variable=paste0(model.name, ".thresh"),
      showSelected.value=mid.thresh,
      label=sprintf(
        "total too high = %.1f (%d intervals / %d possible)",
        total.residual, intervals, possible$negative)),
              data=total.thresh.pred[sign.residual==1 & model.name=="BIC", ],
              hjust=0)+
    geom_text(aes(
      bic.left, bic.top-bic.space, color=model.name,
      showSelected.variable=paste0(model.name, ".thresh"),
      showSelected.value=mid.thresh,
      label=sprintf("%d intervals correctly predicted", intervals)),
              data=total.thresh.pred[sign.residual==0 & model.name=="BIC", ],
              hjust=0)+
    geom_text(aes(
      bic.left, bic.top-bic.space*2, color=model.name,
      showSelected.variable=paste0(model.name, ".thresh"),
      showSelected.value=mid.thresh,
      label=sprintf(
        "total too low = %.1f (%d intervals / %d possible)",
        total.residual, intervals, possible$positive)),
              data=total.thresh.pred[sign.residual==-1 & model.name=="BIC", ],
              hjust=0)+
    geom_text(aes(
      learned.left, learned.top, color=model.name,
      showSelected.variable=paste0(model.name, ".thresh"),
      showSelected.value=mid.thresh,
      label=sprintf(
        "total too high = %.1f (%d intervals / %d possible)",
        total.residual, intervals, possible$negative)),
              data=total.thresh.pred[sign.residual==1 & model.name=="learned", ],
              hjust=0)+
    geom_text(aes(
      learned.left, learned.top-bic.space, color=model.name,
      showSelected.variable=paste0(model.name, ".thresh"),
      showSelected.value=mid.thresh,
      label=sprintf("%d intervals correctly predicted", intervals)),
              data=total.thresh.pred[sign.residual==0 & model.name=="learned", ],
              hjust=0)+
    geom_text(aes(
      learned.left, learned.top-bic.space*2, color=model.name,
      showSelected.variable=paste0(model.name, ".thresh"),
      showSelected.value=mid.thresh,
      label=sprintf(
        "total too low = %.1f (%d intervals / %d possible)",
        total.residual, intervals, possible$positive)),
              data=total.thresh.pred[sign.residual==-1 & model.name=="learned", ],
              hjust=0)+
    geom_point(aes(feature, min.L, fill=limit),
               data=data.table(limit="min", train.dt[is.finite(min.L),]),
               shape=21)+
    geom_point(aes(feature, max.L, fill=limit),
               data=data.table(limit="max", train.dt[is.finite(max.L),]),
               shape=21)+
    geom_abline(aes(slope=slope,
                    showSelected.variable=paste0(model.name, ".thresh"),
                    showSelected.value=mid.thresh,
                    intercept=intercept, color=model.name),
                data=some.thresh)+
    geom_segment(aes(
      feature, pred.log.lambda, xend=feature, color=model.name,
      showSelected.variable=paste0(model.name, ".thresh"),
      showSelected.value=mid.thresh,
      yend=pred.log.lambda-residual),
                 size=1,
                 data=residual.thresh.pred)+
    scale_fill_manual(values=c(min="black", max="white"))+
    scale_color_manual(values=c(BIC="red", learned="deepskyblue"))+
    scale_x_continuous("feature = log(log(n = number of data points to segment))")+
    scale_y_continuous("<-- more changes     log(penalty)     less changes -->"),
  first=list())
pred.thresh.only <- some.thresh[threshold=="predicted",]
for(row.i in 1:nrow(pred.thresh.only)){
  r <- pred.thresh.only[row.i, ]
  viz$first[[paste0(r$model.name, ".thresh")]] <- r$mid.thresh
}
animint2dir(viz, "figure-regression-interactive-some")

## TODO: facet_grid in regression plot, with second panel showing
## learned model in log(n) space. or maybe a single panel with
## showSelected=model, key=profile (dots, regression line,
## residuals/margin move to new positions when model changes).

## TODO: plot data with showSelected=profile, model and labels/errors
## with showSelected=model.

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


