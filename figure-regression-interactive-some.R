source("packages.R")
## The other file figure-regression-interactive.R visualizes all of
## the data points, comparing two models (BIC and learned), using
## three plots: threshold, ROC, regression. In contrast this file
## visualizes the same two models for only a few data points, but five
## plots total. The two extra plots show the data and the model
## selection/error functions.
data(neuroblastoma, package="neuroblastoma")
load("Segmentor.models.RData")
ids.str <- paste(c(
  1, 4, 6, 8, 10, 11,
  ##120,#finite coefs and warns
  13,#finite coefs, no warning
  NULL))
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
    if(TRUE){#buggy as of R-4.4.1
      fit <- survreg(
        Surv(min.log.lambda, max.log.lambda, type="interval2") ~ feature,
        train.dt, dist="gaussian")
      print(weight.vec <- coef(fit))
    }
    fit <- penaltyLearning::IntervalRegressionUnregularized(
      feature.mat, target.mat)
    weight.vec <- coef(fit)[,"0"]
    data.table(
      slope=weight.vec[["feature"]],
      intercept=weight.vec[["(Intercept)"]]
      )
  })
model.label.vec <- c(
  BIC="BIC feature=log(log(n))",
  learned="Learned feature=log(var.est)")
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
    negative=sum(-Inf < min.log.lambda),
    positive=sum(max.log.lambda < Inf)
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
  ##pred.dot.list[[model.name]] <- data.table(model.name, model.label, feature.result$thresholds)
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
##pred.dot <- do.call(rbind, pred.dot.list)
auc.polygon <- do.call(rbind, auc.polygon.list)
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ model.label, scales="free")+
  geom_text(aes(
    feature, 4, 
    label=sprintf(
      "total too high = %.1f (%d intervals / %d possible)",
      total.residual, intervals, possible$positive)),
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

## For reproducing an error.
(out.dt <- model.train.dt[, .(
  id=as.integer(paste(profile.id)), lo=min.log.lambda, hi=max.log.lambda, feature)])
dput(data.frame(out.dt))

w <- roc[is.finite(min.thresh), mean(diff(min.thresh))]
roc[, prev.thresh := ifelse(min.thresh==-Inf, max.thresh-w, min.thresh)]
roc[, next.thresh := ifelse(max.thresh==Inf, min.thresh+w, max.thresh)]
roc[, MID.thresh := (prev.thresh+next.thresh)/2]
round.int <- 1
u.MID <- length(unique(roc$MID.thresh))
while({
  roc[, mid.thresh := round(MID.thresh,round.int)][]
  length(unique(roc$mid.thresh)) != u.MID
}){
  round.int <- round.int+1
}

##pred.dot[, mid.thresh := (min.thresh+max.thresh)/2]
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
residual.thresh.pred <- some.thresh.pred
space.dt <- rbind(
  BIC=data.table(model.name="BIC", left=1.55), 
  learned=data.table(model.name="learned", left=-3.3)
)[, `:=`(
  top=4.7, space=0.6
)][]
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
  geom_segment(aes(
    prev.thresh, errors,
    color=model.name,
    xend=next.thresh, yend=errors),
    data=roc)

auc[, TPR := ifelse(model.name=="BIC", 0.5, 0.4)]
auc[, FPR := 0.2]
thresh.colors <- c(predicted="green", "min.error"="black", other="white")
some.profiles[, profile.chrom := paste(profile.id, chromosome)]
some.anns[, profile.chrom := paste(profile.id, chromosome)]
train.dt[, profile.chrom := paste(profile.id, chromosome)]
labels.thresh[, profile.chrom := paste(profile.id, chromosome)]
selection.thresh[, profile.chrom := paste(profile.id, chromosome)]
changes.thresh[, profile.chrom := paste(profile.id, chromosome)]
selection.melt <- melt(
  data.table(errors$model.errors)[, n.segments := as.numeric(n.segments)],
  measure.vars=c("errors", "n.segments"))
selection.melt[, profile.chrom := paste(profile.id, chromosome)]
b.color <- function(label, color, limit){
  data.table(label, color, limit)
}
break.color.dt <- rbind(
  b.color("breakpoint","#a445ee","max"),
  b.color("normal","orange","min"))
breakpoint.colors <- break.color.dt[, structure(color, names=label)]
min.max.colors <- break.color.dt[, structure(color, names=limit)]
train.long <- nc::capture_melt_single(
  train.dt,
  limit="min|max",
  ".log.lambda"
)[
  is.finite(value)
][
  break.color.dt, on="limit"
]

rect.w <- 0.05
viz <- animint(
  title="BIC versus learned penalty in neuroblastoma data",
  selection=ggplot()+
    ggtitle("Selected problem, label error")+
    theme_bw()+
    theme(panel.margin=grid::unit(1, "lines"))+
    theme_animint(width=300, height=300)+
    facet_grid(variable ~ ., scales="free")+
    xlab("log(penalty)")+
    ylab("")+
    geom_vline(aes(
      xintercept=pred.log.lambda+mid.thresh,
      color=model.name,
      key=model.name),
      help="Predicted threshold shown for each model.",
      showSelected=c(selector="mid.thresh", "profile.chrom"),
      clickSelects="model.name",
      data=selection.thresh[
      , selector := paste0(model.name,'.thresh')
      ],
      alpha=0.8,
      size=4)+
    guides(color="none", fill="none")+
    geom_segment(aes(
      min.log.lambda,
      value,
      xend=max.log.lambda,
      yend=value),
      help="Number of label errors (top) and number of segments (bottom), for the selected data sequence.",
      data=selection.melt,
      showSelected="profile.chrom")+
    geom_point(aes(
      value, 0, fill=label),
      data=data.table(train.long[model.name=="BIC"], variable="errors"),
      size=4,
      help="Dot shows limit of target interval of acceptable penalty values, used as the output when training the regression model.",
      showSelected=c("label", "profile.chrom"))+
    scale_fill_manual(values=breakpoint.colors)+
    scale_color_manual(values=model.colors),
  profile=ggplot()+
    ggtitle("Selected problem, data sequence")+
    theme_bw()+
    theme_animint(height=300, width=650)+
    ylab("Log-ratio (DNA copy number)")+
    xlab("Position on chromosome (mega bases)")+
    scale_fill_manual(values=breakpoint.colors)+
    animint2::geom_tallrect(aes(#not penaltyLearning::geom_tallrect!
      xmin=min/1e6, xmax=max/1e6,
      fill=label),
      alpha=0.3,
      help="Colored rectangle represents a label, which should contain either at least one change (purple breakpoint label), or no change (orange normal label).",
      showSelected="profile.chrom",
      data=some.anns[, label := annotation])+
    scale_linetype_manual("error type", values=c(
      correct=0,
      "false negative"=3,
      "false positive"=1))+
    scale_color_manual(values=model.colors)+
    scale_size_manual(values=c(
      BIC=4,
      learned=1.5))+
    guides(color="none", size="none")+
    animint2::geom_tallrect(aes(
      xmin=min/1e6, xmax=max/1e6,
      linetype=status,
      size=model.name,
      color=model.name),
      help="Border of rectangle shows if there is a label error. False negative occurs for no predicted change-points in a purple breakpoint label, and false positive occurs for at least one predicted change-point in an orange normal label.",
      fill=NA,
      showSelected=c(selector="mid.thresh", "profile.chrom"),
      data=labels.thresh[
      , selector := paste0(model.name,".thresh")
      ])+
    geom_segment(aes(
      chromStart/1e6, 1,
      xend=chromStart/1e6, yend=-1,
      key=paste(model.name, chromStart),
      color=model.name, size=model.name),
      showSelected=c(selector="mid.thresh","profile.chrom"),
      linetype="dashed",
      help="Dashed vertical segments represent predicted change-points.",
      data=changes.thresh[
      , selector := paste0(model.name,".thresh")
      ])+
    geom_point(aes(
      position/1e6, logratio),
      showSelected="profile.chrom",
      fill=NA,
      help="Black dots represent data sequence to segment.",
      data=some.profiles),
  thresholds=ggplot()+
    ggtitle("Total error, select threshold")+
    theme_bw()+
    theme_animint(height=300, width=250)+
    guides(fill="none", color="none")+
    animint2::geom_tallrect(aes(
      xmin=offset-rect.w, xmax=offset+rect.w,
      fill=label),
      size=0.5,
      help="Vertical line shows limit of target interval of acceptable penalty values, with respect to the label in selected data sequence.",
      data=train.long[, offset := value-pred.log.lambda],
      showSelected=c("model.name","label","profile.chrom"))+
    animint2::geom_tallrect(aes(
      xmin=min.thresh, xmax=max.thresh),
      alpha=0.2,
      help="Grey rectangle shows selected threshold.",
      clickSelects=c(selector="mid.thresh"),
      showSelected="model.name",
      data=roc[
      , selector := paste0(model.name,".thresh")
      ])+
    geom_vline(aes(
      xintercept=mid.thresh, color=model.name,
      key=model.name),
      help="Vertical line shows selected threshold.",
      showSelected=c(selector="mid.thresh", "model.name"),
      data=roc)+
    geom_segment(aes(
      min.thresh, errors,
      color=model.name,
      xend=max.thresh, yend=errors),
      clickSelects="model.name",
      help="Segments represent total number of label errors, as a function of constant/threshold added to log(penalty).",
      alpha=0.8,
      size=4,
      data=roc)+
    ylab("total incorrect labels")+
    xlab("Constant added to log(penalty)")+
    scale_fill_manual(values=c(thresh.colors, breakpoint.colors))+
    scale_color_manual(values=model.colors)+
    geom_point(aes(
      x=(max.thresh+min.thresh)/2, errors,
      color=model.name,
      fill=threshold),
      help="Dots represent predicted model, and min error models.",
      showSelected=c("model.name","threshold"),
      clickSelects=c(selector="mid.thresh"),
      data=roc[threshold!="other",],
      size=4,
      shape=21),
  regression=ggplot()+
    theme_bw()+
    theme(
      legend.position="none",
      panel.margin=grid::unit(1, "lines"))+
    facet_grid(. ~ model.label, scales="free")+
    theme_animint(height=300, width=450)+
    geom_text(aes(
      left, log.penalty, 
      color=model.name,
      key=paste(model.name, sign.residual),
      label=ifelse(sign.residual==0, sprintf(
        "%d intervals correctly predicted", intervals), sprintf(
        "total too %s = %.1f (%d/%d labels)",
          ifelse(sign.residual==1, "high", "low"),
          total.residual, intervals,
          possible[, ifelse(sign.residual==1, positive, negative)]))),
      help="Error rates with respect to the target interval of acceptable penalty values.",
      data=total.labels[
      , selector := paste0(model.name,".thresh")
      ],
      showSelected=c(selector="mid.thresh"),
      hjust=0)+
    geom_abline(aes(
      slope=slope,
      key=model.name,
      color=model.name,
      intercept=intercept+mid.thresh),
      showSelected=c(selector="mid.thresh"),
      help="Regression lines, for predicting log(penalty) values.",
      data=all.reg.lines[
      , selector := paste0(model.name,'.thresh')
      ])+
    geom_segment(aes(
      feature, pred.plus.thresh,
      xend=feature, yend=pred.plus.thresh-residual, 
      key=paste(model.name, profile.id, chromosome),
      color=model.name),
      size=1,
      help="Segments represent residuals/errors with respect to the target interval of acceptable penalty values.",
      showSelected=c(selector="mid.thresh"),
      data=residual.thresh.pred[
      , selector := paste0(model.name,'.thresh')
      ])+
    geom_point(aes(
      feature, pred.plus.thresh,
      key=paste(model.name, profile.id, chromosome),
      size=status,
      color=model.name),
      help="A dot is drawn at the predicted log(penalty) value, for each data sequence with a prediction error, with respect to the target interval of acceptable penalty values.",
      showSelected=c(selector="mid.thresh"),
      data=residual.thresh.pred[
      , status := ifelse(residual==0, "correct", "error")
      ])+
    scale_size_manual(values=c(correct=0, error=2))+
    geom_point(aes(
      feature, value,
      fill=label),
      showSelected="label",
      help="Dots represent limits of the target interval of acceptable penalty values.",
      clickSelects="profile.chrom",
      data=train.long,
      alpha=0.8,
      size=4,
      shape=21)+
    guides(color="none")+
    scale_fill_manual(values=breakpoint.colors)+
    scale_color_manual(values=model.colors)+
    scale_x_continuous("feature")+
    scale_y_continuous("log(penalty)", limits=c(-5,5)),
  roc=ggplot()+
    ggtitle("ROC curves, select threshold")+
    theme_bw()+
    theme_animint(height=300, width=250)+
    geom_text(aes(
      FPR, TPR,
      color=model.name,
      label=sprintf("AUC = %.4f", auc)),
      help="AUC = Area Under the ROC Curve.",
      clickSelects="model.name",
      hjust=0,
      data=auc)+
    geom_path(aes(
      FPR, TPR,
      color=model.name,
      group=model.name),
      clickSelects="model.name",
      help="The Receiver Operating Characteristic (ROC) curve is a plot of True Positive Rate versus False Positive Rate, as the prediction threshold is varied.",
      size=4,
      alpha=0.8,
      data=roc)+
    scale_x_continuous(
      "False Positive Rate",
      breaks=seq(0, 1, by=0.2))+
    scale_y_continuous(
      "True Positive Rate",
      breaks=seq(0, 1, by=0.2))+
    scale_fill_manual(values=thresh.colors, breaks=names(thresh.colors))+
    scale_color_manual(values=model.colors, breaks=names(model.colors))+
    geom_text(aes(
      FPR+0.05, TPR,
      label=sprintf(
        "Predicted FPR=%.2f TPR=%.2f",
        FPR, TPR)),
      hjust=0,
      help="FPR and TPR at predicted threshold.",
      data=roc[threshold=="predicted"],
      showSelected=c("model.name","threshold"),
      color=thresh.colors[["predicted"]])+
    geom_point(aes(
      FPR, TPR,
      color=model.name, fill=threshold),
      clickSelects=c(selector="mid.thresh"),
      showSelected="model.name",
      help="Each dot on the ROC curve corresponds to an interval of thresholds, or constants which could be added to the predicted scores.",
      data=roc,
      size=4,
      alpha=0.7,
      shape=21)+
    coord_equal(),
  first=list(profile.chrom="8 2"),
  duration=list(),
  video="https://vimeo.com/1051547245",
  out.dir="figure-regression-interactive-some",
  source="https://github.com/tdhock/change-tutorial/blob/master/figure-regression-interactive-some.R",
  selector.types=list(model.name="single"))
pred.thresh.only <- roc[min.thresh < 0 & 0 < max.thresh,]
for(row.i in 1:nrow(pred.thresh.only)){
  r <- pred.thresh.only[row.i, ]
  viz$first[[paste0(r$model.name, ".thresh")]] <- r$mid.thresh
  viz$duration[[paste0(r$model.name, ".thresh")]] <- 2000
}

viz
if(FALSE){
  animint2pages(viz, "2024-08-bic-learned-details")
  ## https://rcdata.nau.edu/genomic-ml/animint-gallery/2017-02-13-Learned-penalty-function-vs-BIC/index.html
}
