source("packages.R")

data(neuroblastoma, package="neuroblastoma")
load("Segmentor.models.RData")

ids.str <- paste(c(1, 4, 6, 8, 10, 11))
selection <-
  Segmentor.models$loss[profile.id %in% ids.str, modelSelection(
    .SD, complexity="n.segments"
  ), by=.(profile.id, chromosome)]
some.anns <- data.table(
  neuroblastoma$annotations)[profile.id %in% ids.str,]
changes <- Segmentor.models$segs[
  1 < start & profile.id %in% ids.str & chromosome %in% some.anns$chromosome, ]
errors <- labelError(
  selection, some.anns, changes,
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

model.fun.list <- list(
  BIC=function(train.dt){
    data.table(slope=1, intercept=0)
  }, learned=function(train.dt){
    target.mat <- train.dt[, cbind(min.log.lambda, max.log.lambda)]
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
roc.list <- list()
auc.list <- list()
pred.dot.list <- list()
auc.polygon.list <- list()
for(model.name in names(feature.name.vec)){
  model.label <- model.label.vec[[model.name]]
  feature.name <- feature.name.vec[[model.name]]
  model.train.dt <- data.table(
    model.name, model.label,
    target.dt,
    feature=feature.dt[[feature.name]])
  model.fun <- model.fun.list[[model.name]]
  reg.dt <- data.table(model.fun(model.train.dt), model.name, model.label)
  reg.line.list[[model.name]] <- reg.dt
  model.train.dt[, pred.log.lambda := reg.dt[, feature*slope+intercept] ]
  ## This residual is for the non-interactive plot.
  model.train.dt[, residual := {
    targetIntervalResidual(
      cbind(min.log.lambda, max.log.lambda), pred.log.lambda)
  }]
  possible <- model.train.dt[, list(
    positive=sum(-Inf < min.log.lambda),
    negative=sum(max.log.lambda < Inf)
    )]
  min.feature <- min(model.train.dt$feature)
  totals.list[[model.name]] <- model.train.dt[, list(
    feature=min.feature,
    total.residual=sum(residual),
    intervals=.N
    ), by=.(model.name, model.label, sign.residual=sign(residual))]
  train.dt.list[[model.name]] <- model.train.dt
  feature.result <- ROChange(
    errors$model.errors, model.train.dt, c("profile.id", "chromosome"))
  pred.dot.list[[model.name]] <- data.table(
    model.name, model.label, feature.result$thresholds)
  feature.result$roc[, threshold := {
    ifelse(
      min.thresh < 0 & 0 < max.thresh,
      "predicted", ifelse(
        errors==min(errors),
        "min.error", "other"))
  }]
  roc.list[[model.name]] <- data.table(
    model.name, model.label, feature.result$roc)
  auc.list[[model.name]] <- data.table(
    model.name, model.label, auc=feature.result$auc)
  auc.polygon.list[[model.name]] <- data.table(
    model.name, model.label, feature.result$auc.polygon)
}
train.dt <- do.call(rbind, train.dt.list)
totals <- do.call(rbind, totals.list)
reg.line <- do.call(rbind, reg.line.list)
roc <- do.call(rbind, roc.list)
auc <- do.call(rbind, auc.list)
pred.dot <- do.call(rbind, pred.dot.list)
auc.polygon <- do.call(rbind, auc.polygon.list)
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
  geom_point(aes(feature, ifelse(is.finite(min.log.lambda), min.log.lambda, NA), fill="min"), data=train.dt, shape=21)+
  geom_point(aes(feature, ifelse(is.finite(max.log.lambda), max.log.lambda, NA), fill="max"), data=train.dt, shape=21)+
  scale_fill_manual("limit", values=c(min="black", max="white"))+
  scale_x_continuous("feature")+
  scale_y_continuous("<-- more changes     log(penalty)     less changes -->")
print(gg)

w <- roc[is.finite(min.thresh), abs(mean(diff(min.thresh)))]
roc[, prev.thresh := ifelse(min.thresh==-Inf, max.thresh-w, min.thresh)]
roc[, next.thresh := ifelse(max.thresh==Inf, min.thresh+w, max.thresh)]
roc[, mid.thresh := (prev.thresh+next.thresh)/2]
pred.dot[, mid.thresh := (min.thresh+max.thresh)/2]
some.thresh.pred <- train.dt[roc, on=.(
  model.name, model.label), allow.cartesian=TRUE]
some.thresh.pred[, pred.plus.thresh := pred.log.lambda + mid.thresh]
selection.thresh <- selection[some.thresh.pred, on=.(
  profile.id, chromosome,
  min.log.lambda <= pred.plus.thresh,
  max.log.lambda >= pred.plus.thresh)]
stopifnot(nrow(selection.thresh)==nrow(some.thresh.pred))
labels.thresh <- errors$label.errors[selection.thresh, on=.(
  profile.id, chromosome, n.segments)]
stopifnot(nrow(selection.thresh)==nrow(labels.thresh))
changes.thresh <- changes[selection.thresh, on=.(
  profile.id, chromosome, n.segments), nomatch=0L]


some.thresh.pred[, residual := {
  targetIntervalResidual(cbind(min.log.lambda, max.log.lambda), pred.plus.thresh)
}]
total.thresh.pred <- some.thresh.pred[, list(
  total.residual=sum(residual),
  intervals=.N
  ), by=.(model.name, model.label, sign.residual=sign(residual), mid.thresh)]
total.thresh.pred
residual.thresh.pred <- some.thresh.pred[residual!=0,]
space.dt <- rbind(
  BIC=data.table(model.name="BIC", left=1.55, top=4.5, space=0.5),
  learned=data.table(model.name="learned", left=-3.3, top=4.5, space=0.5))
setkey(space.dt, model.name)
setkey(total.thresh.pred, model.name)
total.labels <- total.thresh.pred[space.dt]
total.labels[, log.penalty := top-(1-sign.residual)*space]
all.reg.lines <- reg.line[roc, on=.(model.name, model.label)]

model.colors <- c(
  learned="blue",
  BIC="red"
  )
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  ##facet_grid(. ~ panel.label, scales="free")+
  scale_color_manual(values=model.colors)+
  geom_segment(aes(prev.thresh, errors,
                   color=model.name,
                   xend=next.thresh, yend=errors),
               data=roc)

auc[, TPR := c(0.85, 0.9)]
auc[, FPR := 0.2]
thresh.colors <- c(predicted="green", "min.error"="black", other="white")
some.profiles[, profile.chrom := paste(profile.id, chromosome)]
some.anns[, profile.chrom := paste(profile.id, chromosome)]
train.dt[, profile.chrom := paste(profile.id, chromosome)]
labels.thresh[, profile.chrom := paste(profile.id, chromosome)]
selection.thresh[, profile.chrom := paste(profile.id, chromosome)]
changes.thresh[, profile.chrom := paste(profile.id, chromosome)]
selection.melt <- melt(
  errors$model.errors, measure.vars=c("errors", "n.segments"))
selection.melt[, profile.chrom := paste(profile.id, chromosome)]
breakpoint.colors <- c(
  "breakpoint"="#a445ee",
  "normal"="#f6f4bf")
viz <- list(
  title="BIC versus learned penalty in neuroblastoma data",
  profile=ggplot()+
    ggtitle("Selected problem")+
    theme_bw()+
    theme_animint(height=300, width=400)+
    ylab("logratio (approximate copy number)")+
    xlab("position on chromosome")+
    scale_fill_manual("label", values=breakpoint.colors)+
    geom_tallrect(aes(xmin=min/1e6, xmax=max/1e6,
                      showSelected=profile.chrom,
                      fill=annotation),
                  alpha=0.3,
                  data=some.anns)+
    scale_linetype_manual("error type", values=c(
      correct=0,
      "false negative"=3,
      "false positive"=1))+
    scale_color_manual(values=model.colors)+
    scale_size_manual(values=c(
      BIC=4,
      learned=1.5))+
    guides(color="none", size="none")+
    geom_tallrect(aes(xmin=min/1e6, xmax=max/1e6,
                      linetype=status,
                      size=model.name,
                      color=model.name,
                      showSelected.variable=paste0(model.name, ".thresh"),
                      showSelected.value=mid.thresh,
                      showSelected=profile.chrom),
                  fill=NA,
                  data=labels.thresh)+
    geom_segment(aes(chromStart/1e6, 1,
                     xend=chromStart/1e6, yend=-1,
                     showSelected.variable=paste0(model.name, ".thresh"),
                     showSelected.value=mid.thresh,
                     showSelected=profile.chrom,
                     color=model.name, size=model.name),
               data=changes.thresh)+
    geom_point(aes(position/1e6, logratio,
                   showSelected=profile.chrom),
               fill=NA,
               data=some.profiles),
  thresholds=ggplot()+
    ggtitle("Total error, select threshold")+
    theme_bw()+
    theme_animint(height=300, width=300)+
    guides(fill="none", color="none")+
    geom_tallrect(aes(xmin=prev.thresh, xmax=next.thresh,
                      clickSelects.variable=paste0(model.name, ".thresh"),
                      clickSelects.value=mid.thresh,
                      showSelected=model.name),
                  alpha=0.2,
                  data=roc)+
    geom_vline(aes(xintercept=mid.thresh, color=model.name,
                   showSelected.variable=paste0(model.name, ".thresh"),
                   showSelected.value=mid.thresh,
                   showSelected=model.name),
               data=roc)+
    geom_segment(aes(prev.thresh, errors,
                     color=model.name,
                     clickSelects=model.name,
                     xend=next.thresh, yend=errors),
                 alpha=0.8,
                 size=4,
                 data=roc)+
    ylab("total incorrect labels")+
    xlab("constant/threshold added to log(penalty)")+
    scale_fill_manual(values=thresh.colors)+
    scale_color_manual(values=model.colors)+
    geom_point(aes((max.thresh+min.thresh)/2, errors,
                   clickSelects.variable=paste0(model.name, ".thresh"),
                   clickSelects.value=mid.thresh,
                   color=model.name,
                   showSelected=model.name,
                   showSelected2=threshold,
                   fill=threshold),
               data=roc[threshold!="other",],
               size=4,
               shape=21),
  roc=ggplot()+
    ggtitle("ROC curves, select threshold")+
    theme_bw()+
    theme_animint(height=300, width=300)+
    ##guides(color="none")+
    geom_text(aes(FPR, TPR,
                  clickSelects=model.name,
                  color=model.name, label=sprintf("AUC = %.4f", auc)),
              hjust=0,
              data=auc)+
    geom_path(aes(FPR, TPR, clickSelects=model.name,
                  color=model.name, group=model.name),
              size=4,
              alpha=0.8,
              data=roc)+
    scale_fill_manual(values=thresh.colors, breaks=names(thresh.colors))+
    scale_color_manual(values=model.colors, breaks=names(model.colors))+
    ## geom_point(aes(FPR, TPR,
    ##                showSelected.variable=paste0(model.name, ".thresh"),
    ##                showSelected.value=mid.thresh,
    ##                showSelected2=threshold,
    ##                showSelected=model.name),
    ##            data=roc,
    ##            color="grey",
    ##            size=6,
    ##            shape=21)+
    geom_point(aes(FPR, TPR,
                   showSelected=model.name,
                   clickSelects.variable=paste0(model.name, ".thresh"),
                   clickSelects.value=mid.thresh,
                   color=model.name, fill=threshold),
               data=roc,
               size=4,
               alpha=0.7,
               shape=21)+
    coord_equal(),
  regression=ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(. ~ model.label, scales="free")+
    theme_animint(height=300, width=500)+
    geom_text(aes(
      left, log.penalty, 
      color=model.name,
      showSelected.variable=paste0(model.name, ".thresh"),
      showSelected.value=mid.thresh,
      label=ifelse(sign.residual==0, sprintf(
        "%d intervals correctly predicted", intervals), sprintf(
        "total too %s = %.1f (%d/%d intervals)",
          ifelse(sign.residual==1, "high", "low"),
          total.residual, intervals,
          possible[, ifelse(sign.residual==1, negative, positive)]))),
              data=total.labels,
              hjust=0)+
    geom_abline(aes(slope=slope,
                    ##showSelected=model.name,
                    color=model.name,
                    showSelected.variable=paste0(model.name, ".thresh"),
                    showSelected.value=mid.thresh,
                    intercept=intercept+mid.thresh),
                data=all.reg.lines)+
    geom_segment(aes(
      feature, pred.plus.thresh,
      xend=feature, yend=pred.plus.thresh-residual, 
      ##showSelected=model.name,
      color=model.name,
      showSelected.variable=paste0(model.name, ".thresh"),
      showSelected.value=mid.thresh),
                 size=1,
                 data=residual.thresh.pred)+
    geom_point(aes(feature, min.log.lambda,
                   ##showSelected=model.name,
                   clickSelects=profile.chrom,
                   fill=limit),
               data=data.table(
                 limit="min", train.dt[is.finite(min.log.lambda),]),
               alpha=0.8,
               size=4,
               shape=21)+
    geom_point(aes(feature, max.log.lambda,
                   clickSelects=profile.chrom,
                   ##showSelected=model.name,
                   fill=limit),
               data=data.table(
                 limit="max", train.dt[is.finite(max.log.lambda),]),
               size=4,
               alpha=0.8,
               shape=21)+
    geom_point(aes(feature,
                   pred.log.lambda+mid.thresh,
                   showSelected=profile.chrom,
                   showSelected.variable=paste0(model.name, ".thresh"),
                   showSelected.value=mid.thresh,
                   color=model.name),
               data=selection.thresh)+
    guides(color="none")+
    scale_fill_manual(values=c(min="black", max="white"))+
    scale_color_manual(values=model.colors)+
    scale_x_continuous("feature")+
    scale_y_continuous("<- more changes log(penalty) less changes ->"),
  selection=ggplot()+
    ggtitle("Selected problem error")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    theme_animint(width=300, height=300)+
    facet_grid(variable ~ ., scales="free")+
    xlab("log(penalty)")+
    ylab("")+
    geom_vline(aes(
      xintercept=pred.log.lambda+mid.thresh,
      color=model.name,
      showSelected.variable=paste0(model.name, ".thresh"),
      showSelected.value=mid.thresh,
      showSelected=profile.chrom,
      clickSelects=model.name),
      data=selection.thresh,
      alpha=0.8,
      size=4)+
    guides(color="none")+
    geom_segment(aes(
      ifelse(min.log.lambda==-Inf, min(max.log.lambda)-0.5, min.log.lambda),
      value,
      xend=ifelse(max.log.lambda==Inf, max(min.log.lambda)+0.5, max.log.lambda),
      yend=value,
      showSelected=profile.chrom),
      data=selection.melt)+
    scale_color_manual(values=model.colors),
  first=list(),
  selector.types=list(model.name="single"))
pred.thresh.only <- roc[min.thresh < 0 & 0 < max.thresh,]
for(row.i in 1:nrow(pred.thresh.only)){
  r <- pred.thresh.only[row.i, ]
  viz$first[[paste0(r$model.name, ".thresh")]] <- r$mid.thresh
}
animint2dir(viz, "figure-regression-interactive-some")

