library(animint2)
library(data.table)
data(neuroblastoma, package="neuroblastoma")

nb.dts <- lapply(neuroblastoma, data.table, key=c("profile.id", "chromosome"))

if(file.exists("Segmentor.models.RData")){
  load("Segmentor.models.RData")
}else{
  segs.list <- list()
  loss.list <- list()
  ann.dt <- nb.dts[["annotations"]]
  for(problem.i in 1:nrow(ann.dt)){
    meta <- ann.dt[problem.i, .(profile.id, chromosome)]
    cat(sprintf(
      "%4d / %4d problems\n",
      problem.i, nrow(ann.dt)))
    pro <- nb.dts[["profiles"]][meta]
    max.segments <- 20
    fit <- jointseg::Fpsn(pro$logratio, max.segments)
    rss.vec <- rep(NA, max.segments)
    cum.vec <- c(0, cumsum(pro$logratio))
    for(n.segments in 1:max.segments){
      end <- fit[["t.est"]][n.segments, 1:n.segments]
      data.before.change <- end[-n.segments]
      data.after.change <- data.before.change+1
      pos.before.change <- as.integer(
      (pro$position[data.before.change]+pro$position[data.after.change])/2)
      start <- c(1, data.after.change)
      chromStart <- c(pro$position[1], pos.before.change)
      chromEnd <- c(pos.before.change, max(pro$position))
      seg.mean.vec <- (cum.vec[end+1]-cum.vec[start])/(end-start+1)
      data.mean.vec <- rep(seg.mean.vec, end-start+1)
      residual.vec <- pro$logratio - data.mean.vec
      rss.vec[n.segments] <- sum(residual.vec^2)
      segs.list[[paste(problem.i, n.segments)]] <- data.table(
        meta,
        n.segments,
        start,
        end,
        chromStart,
        chromEnd,
        mean=seg.mean.vec)
    }
    loss.list[[paste(problem.i, n.segments)]] <- data.table(
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
Segmentor.models$changes <- Segmentor.models$segs[1 < start]
errors <- penaltyLearning::labelError(
  selection, neuroblastoma$annotations, Segmentor.models$changes,
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
auc$TPR <- c(0.84, 1)
auc$FPR <- 0.1
auc$hjust <- 1
ggplot()+
  theme_bw()+
  geom_text(aes(
    FPR, TPR,
    color=model.name,
    hjust=hjust,
    label=sprintf("AUC = %.4f", auc)),
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
add_mid <- function(DT){
  DT[, mid.thresh := (min.thresh+max.thresh)/2]
}
add_mid(roc)
add_mid(pred.dot)
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
  grid.roc[pred.dot, mid.thresh := 0, on=c("model.name",rate)]
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
grid.labels <- err.dt[adj.pred, .(
  profile.id, chromosome, n.segments, rate, grid, model.name, tp, fp, status,
  adj.log.lambda=i.adj.log.lambda, min, max
), on=.(
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
    FPR, TPR,
    color=model.name,
    label=sprintf("AUC = %.4f", auc)),
    vjust=0,
    data=auc)+
  geom_segment(aes(
    FPR_AIC, TPR_AIC,
    xend=FPR_BIC, yend=TPR_BIC),
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
  coord_equal(xlim=c(0, 0.1), ylim=c(0.4, 1))

roc.grid.dt[
, slope := ifelse(model.name=="BIC", 1, 0)
][
, intercept := ifelse(model.name=="BIC", 0, log(2))+mid.thresh
]
grid.label.int <- grid.labels[
  target.dt,
  on=.(profile.id,chromosome)
][
  feature.dt,
  on=.(profile.id,chromosome)
][
, residual := penaltyLearning::targetIntervalResidual(
  cbind(min.log.lambda, max.log.lambda), adj.log.lambda
)][
, limit := ifelse(is.finite(min.log.lambda), "min", "max")
]
total.thresh.pred <- grid.label.int[, list(
  total.residual=sum(residual),
  intervals=.N
), by=.(model.name, sign.residual=sign(residual), rate, grid)
][
, top := ifelse(model.name=="AIC", 6.5, -4)
][
, left := ifelse(model.name=="AIC", 1.45, 1.81)
][
, log.penalty := top-(1-sign.residual)*0.7
]
total.thresh.pred[, .(
  intervals=sum(intervals)
), by=.(model.name, rate, grid)]
long_limits <- function(DT){
  not.limit <- names(DT)[names(DT)!="limit"]
  train.long <- nc::capture_melt_single(
    DT[, not.limit, with=FALSE],
    limit="min|max",
    "[.]log[.]lambda",
    value.name="log.penalty"
  )[is.finite(log.penalty)]
}
select_rate <- function(DT){
  if("rate" %in% names(DT)){
    DT[, Rate := sprintf("%s=%.3f", rate, grid)]
  }
  if("profile.id" %in% names(DT)){
    DT[, pid.chr := paste0(profile.id, ".", chromosome)]
  }
}
lapply(nb.dts, select_rate)
select_rate(labeled.profiles)
train.long <- long_limits(train.dt)
select_rate(grid.label.int)
ref.model <- "BIC"#arbitrary.
other.model <- "AIC"
grid.label.dots <- long_limits(grid.label.int[model.name==ref.model])
select_rate(total.thresh.pred)
select_rate(roc.grid.dt)
select_rate(roc.diff.dt)
select_rate(train.long)
pred.dot[, `:=`(
  TPR.text=c(0.77, 0.62),
  hjust=c(0,1),
  FPR.text=FPR+ifelse(model.name=="AIC",0.002,0)
)]
model.colors <- c(BIC="red", AIC="deepskyblue")
pixels <- 450
roc.diff.dt[, other := ifelse(rate=="FPR", "TPR", "FPR")]
get_other <- function(model.name){
  data.table(
    model.name, roc.diff.dt
  )[, label := {
    col.name <- paste0(other,"_",model.name)
    sprintf("%s = %.4f", col.name, get(col.name))
  }, by=.(model.name,other)][]
}
sign2col=c("1"="positive","0"="total","-1"="negative")
sign2rate=c("1"="FN","0"="correct","-1"="FP")
total.thresh.pred[, label := sprintf(
  "total %s = %d / %d possible",
  sign2rate[paste(sign.residual)],
  intervals,
  unlist(as.list(possible)[sign2col[paste(sign.residual)]])
)]
one.y.vec <- seq(0.4, 0.74, by=0.02)
one.x.vec <- seq(0.1, 0.03, by=-0.013)
grid.label.not.same <- grid.label.dots[same==FALSE]
setkey(grid.label.not.same, Rate, limit)
grid.label.not.same[, `:=`(
  label = ifelse(limit=="max", "breakpoint", "normal"),
  y=rep(one.y.vec, l=.N),
  x=rep(one.x.vec, each=length(one.y.vec))[1:.N]
), by=Rate]
label.colors <- c(
  breakpoint="purple",
  normal="orange")
limit2label <- c(max="breakpoint", min="normal")
limit.colors <- structure(label.colors[limit2label], names=names(limit2label))
ggplot()+
  geom_text(aes(
    x, y, label=pid.chr, color=label),
    hjust=1,
    data=grid.label.not.same)+
  scale_color_manual(values=label.colors)+
  facet_wrap("Rate")
grid.label.not.same[
, model.name := ifelse(status=="correct", ref.model, other.model)]
selected.sizes <- selection[
  adj.pred,
  .(rate, grid, profile.id, chromosome, model.name, n.segments),
  on=.(
    profile.id, chromosome,
    min.log.lambda<adj.log.lambda,
    max.log.lambda>adj.log.lambda
  )
]
Segmentor.selected <- lapply(Segmentor.models, function(DT){
  join.dt <- DT[
    selected.sizes, on=.(profile.id, chromosome, n.segments), nomatch=0L]
  select_rate(join.dt)
})
select_rate(adj.pred)
adj.pred[, `:=`(
  adj.lambda = exp(adj.log.lambda),
  norm_position = ifelse(model.name=="AIC", 0, 1),
  hjust=ifelse(model.name=="AIC", 0, 1)
)]
before.norm.list <- c(
  Segmentor.selected,
  list(profiles=labeled.profiles,
       annotations=nb.dts[["annotations"]]))
norm.list <- list(
  position=c("position","min","max","chromStart","chromEnd"),
  logratio=c("logratio","mean"))
norm.dt <- dcast(
  labeled.profiles,
  pid.chr ~ .,
  list(min,max),
  value.var=names(norm.list))
setkey(norm.dt, pid.chr)
after.norm.list <- list()
for(data.name in names(before.norm.list)){
  DT <- before.norm.list[[data.name]]
  join.dt <- norm.dt[DT, on="pid.chr"]
  for(norm.name in names(norm.list)){
    gcol <- function(m)join.dt[[paste0(norm.name,"_",m)]]
    m <- gcol("min")
    M <- gcol("max")
    denominator <- M-m
    for(col.name in norm.list[[norm.name]]){
      if(col.name %in% names(join.dt)){
        j <- paste0("norm_", col.name)
        before <- join.dt[[col.name]]
        value <- (before-m)/denominator
        set(join.dt, j=j, value=value)
      }
    }
  }
  after.norm.list[[data.name]] <- join.dt
}
geom_text_label <- function(L){
  geom_text(aes(
    x, y,
    key=pid.chr,
    label=pid.chr),
    hjust=1,
    showSelected=c("Rate","model.name","limit","same"),
    color=label.colors[[L]],
    color_off=label.colors[[L]],
    clickSelects="pid.chr",
    data=grid.label.not.same[label==L])
}
viz <- animint(
  title="AIC/BIC change-point detection comparison using ROC curves",
  out.dir="figure-differences-interactive",
  roc=ggplot()+
    ggtitle("ROC curves, select comparison")+
    xlab("False Positive Rate")+
    ylab("True Positive Rate")+
    theme_bw()+
    theme_animint(width=pixels, height=pixels)+
    theme(legend.position="none")+
    scale_fill_manual(values=model.colors)+
    scale_color_manual(values=model.colors)+
    geom_text(aes(
      x, y, label=label, hjust=hjust),
      data=rbind(
        data.table(label=c(
          "Black text shows number",
          "of label error differences"),
          x=0, y=c(0.88, 0.86), hjust=0),
        data.table(label=c(
          "Label error differences:",
          "(dot shows correct model)"),
          x=0.1, y=c(0.78, 0.76), hjust=1)))+
    geom_point(aes(
      x, y,
      color=model.name,
      key=pid.chr),
      showSelected=c("Rate","model.name","limit","same"),
      clickSelects="pid.chr",
      data=grid.label.not.same)+
    geom_text_label("breakpoint")+
    geom_text_label("normal")+
    geom_text(aes(
      FPR, TPR,
      color=model.name,
      hjust=hjust,
      label=sprintf("%s AUC = %.4f", model.name, auc)),
      showSelected="model.name",
      vjust=0,
      data=auc)+
    geom_path(aes(
      FPR, TPR,
      color=model.name,
      group=model.name),
      showSelected="model.name",
      data=roc)+
    geom_segment(aes(
      FPR, TPR.text,
      xend=FPR, yend=TPR),
      showSelected="model.name",
      data=pred.dot)+
    geom_text(aes(
      FPR.text, TPR.text+0.02,
      label=paste(errors, "errors"),
      hjust=hjust,
      color=model.name),
      showSelected="model.name",
      data=pred.dot)+
    geom_text(aes(
      FPR, TPR.text,
      label=paste(model.name, "default"),
      hjust=hjust,
      color=model.name),
      showSelected="model.name",
      data=pred.dot)+
    geom_segment(aes(
      FPR_AIC, TPR_AIC,
      xend=FPR_BIC, yend=TPR_BIC),
      size=4,
      alpha=0.7,
      clickSelects="Rate",
      data=roc.diff.dt)+
    geom_point(aes(
      FPR, TPR,
      fill=model.name),
      showSelected="model.name",
      color="black",
      data=pred.dot,
      size=3,
      shape=21)+
    geom_text(aes(
      0, 1,
      key="SelectedRate",
      label=paste("Selected", Rate)),
      hjust=0,
      showSelected="Rate",
      data=roc.diff.dt)+
    geom_text(aes(
      0, 0.97,
      key="AIC",
      color=model.name,
      label=label),
      hjust=0,
      showSelected=c("Rate","model.name"),
      data=get_other("AIC"))+
    geom_text(aes(
      0, 0.94,
      key="BIC",
      color=model.name,
      label=label),
      hjust=0,
      showSelected=c("Rate","model.name"),
      data=get_other("BIC"))+
    geom_text(aes(
      FPR_AIC,
      TPR_BIC-ifelse(vjust==1, 0.02, 0),
      label=diff,
      hjust=hjust),
      clickSelects="Rate",
      data=roc.diff.dt)+
    coord_cartesian(xlim=c(0, 0.1), ylim=c(0.4, 1)),
  regression=ggplot()+
    ggtitle("Predicted penalties and errors in BIC space")+
    theme_bw()+
    theme_animint(width=600, height=pixels)+
    geom_abline(aes(
      slope=slope,
      intercept=intercept,
      key=model.name,
      color=model.name),
      showSelected="Rate",
      data=roc.grid.dt)+
    geom_text(aes(
      left, log.penalty, 
      color=model.name,
      key=paste(model.name, sign.residual),
      label=label),
      data=total.thresh.pred,
      showSelected="Rate",
      hjust=0)+
    geom_segment(aes(
      log2.n, adj.log.lambda,
      alpha=same,
      color=model.name,
      key=paste(pid.chr, model.name),
      size=model.name,
      xend=log2.n,
      yend=adj.log.lambda-residual),
      showSelected=c("Rate","limit"),
      data=grid.label.int)+
    geom_point(aes(
      log2.n, log.penalty,
      fill=limit,
      key=pid.chr,
      alpha=same),
      showSelected="Rate",
      clickSelects="pid.chr",
      color="blue",
      color_off="grey50",
      stroke=2,
      size=3,
      data=grid.label.dots,
      shape=21)+
    scale_alpha_manual(values=c("TRUE"=0.1,"FALSE"=1))+
    scale_fill_manual(values=limit.colors)+
    scale_color_manual(values=model.colors)+
    scale_size_manual(values=c(BIC=2, AIC=1))+
    scale_x_continuous(
      "feature = log(log(n = number of data points to segment))")+
    scale_y_continuous(
      "<-- more changes log(penalty) less changes -->"),
  data=ggplot()+
    ggtitle("Predicted changes for selected profile")+
    theme_bw()+
    theme_animint(
      ##update_axes=c("x","y"),#TODO
      width=1000, height=300)+
    geom_rect(aes(
      xmin=norm_min, xmax=norm_max,
      ymin=0, ymax=1,
      fill=annotation),
      showSelected="pid.chr",
      color="grey",
      alpha=0.5,
      data=after.norm.list[["annotations"]])+
    scale_fill_manual(values=label.colors)+
    geom_point(aes(
      norm_position, norm_logratio),
      showSelected="pid.chr",
      chunk_vars="pid.chr",
      data=after.norm.list[["profiles"]])+
    scale_color_manual(values=model.colors)+
    scale_size_manual(values=c(BIC=5, AIC=2))+
    geom_segment(aes(
      norm_chromStart, norm_mean,
      xend=norm_chromEnd, yend=norm_mean,
      color=model.name,
      key=paste(chromStart, model.name),
      size=model.name),
      showSelected=c("pid.chr","Rate"),
      chunk_vars="pid.chr",
      data=after.norm.list[["segs"]])+
    geom_segment(aes(
      norm_chromStart, 0,
      xend=norm_chromStart, yend=1,
      color=model.name,
      key=paste(chromStart, model.name),
      size=model.name),
      showSelected=c("pid.chr","Rate"),
      chunk_vars="pid.chr",
      data=after.norm.list[["changes"]])+
    geom_text(aes(
      norm_position, 1,
      hjust=hjust,
      color=model.name,
      key=model.name,
      label=sprintf("Adjusted %s penalty = %.4f", model.name, adj.lambda)),
      showSelected=c("pid.chr","Rate"),
      data=adj.pred),
  duration=list(Rate=1000),
  source="https://github.com/tdhock/change-tutorial/blob/master/figure-differences-interactive.R"
)
viz
if(FALSE){
  animint2pages(viz, "2023-11-interval-regression-differences")
}

## TODO binseg(monotonic) vs DP(could be non-monotonic).

## TODO CV.
