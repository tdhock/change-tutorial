source("packages.R")

## Exercise 1: compute label error and plot predicted changepoints for
## another set of profiles and unsupervised penalty.
some.ids <- 140:155
pen.name <- "SIC0"
data(neuroblastoma, package="neuroblastoma")
someProfiles <- function(all.profiles){
  data.table(all.profiles)[paste(profile.id) %in% some.ids]
}
unsupervised.profiles <- someProfiles(neuroblastoma$profiles)
unsupervised.labels <- someProfiles(neuroblastoma$annotations)
unsupervised.models <- unsupervised.profiles[, {
  fit.pelt <- changepoint::cpt.mean(
    logratio, penalty=pen.name, method="PELT")
  end <- fit.pelt@cpts
  before.change <- end[-length(end)]
  after.change <- before.change+1L
  data.table(
    pen.name,
    pen.value=fit.pelt@pen.value,
    changes=list(
    (position[before.change]+position[after.change])/2
    ))
}, by=list(profile.id, chromosome)]
unsupervised.changes <- unsupervised.models[, data.table(
  change=changes[[1]]
), by=list(profile.id, chromosome, pen.name)]
unsupervised.error.list <- penaltyLearning::labelError(
  unsupervised.models, unsupervised.labels, unsupervised.changes,
  problem.vars=c("profile.id", "chromosome"),
  model.vars="pen.name",
  change.var="change")
ggplot()+
  unsupervised.error.list$model.errors[, ggtitle(paste0(
    "Unsupervised ", pen.name, " penalty has ",
    sum(errors),
    " incorrect labels: ",
    sum(fp), " FP + ",
    sum(fn), " FN"))]+
  theme(
    panel.margin=grid::unit(0, "lines"),
    panel.border=element_rect(fill=NA, color="grey50"),
    legend.position="bottom",
    legend.box="horizontal"
  )+
  facet_grid(profile.id ~ chromosome, scales="free", space="free_x")+
  geom_point(aes(position/1e6, logratio),
             data=unsupervised.profiles,
             shape=1)+
  scale_x_continuous(
    "position on chromosome (mega bases)",
    breaks=c(100, 200))+
  scale_y_continuous(
    "logratio (approximate DNA copy number)",
    limits=c(-1,1)*1.1)+
  penaltyLearning::geom_tallrect(aes(
    xmin=min/1e6, xmax=max/1e6, fill=annotation),
    alpha=0.5,
    color=NA,
    data=unsupervised.labels)+
  scale_fill_manual("label", values=c(
    "1change"="#ff7d7d",
    breakpoint="#a445ee",
    normal="#f6c48f"
  ))+
  geom_vline(aes(
    xintercept=change/1e6),
    color="green",
    size=1,
    linetype="dashed",
    data=unsupervised.changes)+
  geom_tallrect(aes(
    xmin=min/1e6, xmax=max/1e6, linetype=status),
    fill=NA,
    data=unsupervised.error.list$label.errors)+
  scale_linetype_manual("error type", values=c(
    correct=0,
    "false negative"=3,
    "false positive"=1))

## Exercise 2: label a profile and visualize which models are
## consistent with your labels. compute the model selection and error
## functions. Change the definition of the labels -- how does the
## target interval change? What happens to the target interval if you
## reduce max.segments?
zoom.profile <- 104     #Exercise: change this value!
zoom.chromosome <- 10 #Exercise: change this value!
label <- function(min, max, annotation){
  data.frame(
    profile.id=paste(zoom.profile),
    chromosome=paste(zoom.chromosome),
    min, max, annotation)
}
zoom.labels <- rbind(# Exercise: change these values!
  label(60e6, 100e6, "normal"),
  ##label(40e6, 60e6, "1change"), 
  label(100e6, 120e6, "breakpoint"))
zoom.pro <- subset(
  neuroblastoma$profiles,
  profile.id==zoom.profile & chromosome==zoom.chromosome)
ggplot()+
  geom_point(aes(position/1e6, logratio),
             data=zoom.pro,
             shape=1)+
  scale_y_continuous(
    "logratio (approximate DNA copy number)")+
  penaltyLearning::geom_tallrect(aes(
    xmin=min/1e6, xmax=max/1e6, fill=annotation),
    alpha=0.5,
    color=NA,
    data=zoom.labels)+
  scale_fill_manual("label", values=c(
    "1change"="#ff7d7d",
    breakpoint="#a445ee",
    normal="#f6c48f"
  ))
max.segments <- 15
fit <- Segmentor3IsBack::Segmentor(
  zoom.pro$logratio, model=2, Kmax=max.segments)
zoom.segs.list <- list()
zoom.loss.vec <- rep(NA, max.segments)
for(n.segments in 1:max.segments){
  end <- fit@breaks[n.segments, 1:n.segments]
  data.before.change <- end[-n.segments]
  data.after.change <- data.before.change+1
  pos.before.change <- as.integer(
  (zoom.pro$position[data.before.change]+
   zoom.pro$position[data.after.change])/2)
  start <- c(1, data.after.change)
  chromStart <- c(zoom.pro$position[1], pos.before.change)
  chromEnd <- c(pos.before.change, max(zoom.pro$position))
  seg.mean.vec <- fit@parameters[n.segments, 1:n.segments]
  zoom.segs.list[[n.segments]] <- data.table(
    profile.id=paste(zoom.profile),
    chromosome=paste(zoom.chromosome),
    n.segments, # model complexity.
    start, # in data points.
    end,
    chromStart, # in bases on chromosome.
    chromEnd,
    mean=seg.mean.vec)
  data.mean.vec <- rep(seg.mean.vec, end-start+1)
  stopifnot(length(data.mean.vec)==nrow(zoom.pro))
  zoom.loss.vec[n.segments] <- sum((zoom.pro$logratio-data.mean.vec)^2)
}
zoom.segs <- do.call(rbind, zoom.segs.list)
zoom.changes <- zoom.segs[1 < start, data.table(
  profile.id, chromosome, n.segments,
  changepoint=chromStart)]
zoom.models <- data.table(
  profile.id=paste(zoom.profile),
  chromosome=paste(zoom.chromosome),
  loss=zoom.loss.vec,
  n.segments=as.numeric(1:max.segments))
zoom.error.list <- penaltyLearning::labelError(
  zoom.models,
  zoom.labels,
  zoom.changes,
  change.var="changepoint",
  problem.vars=c("profile.id", "chromosome"))
ggplot()+
  geom_point(aes(position/1e6, logratio),
             data=zoom.pro,
             shape=1)+
  scale_y_continuous(
    "logratio (approximate DNA copy number)")+
  penaltyLearning::geom_tallrect(aes(
    xmin=min/1e6, xmax=max/1e6, fill=annotation),
    alpha=0.5,
    color=NA,
    data=zoom.labels)+
  scale_fill_manual("label", values=c(
    "1change"="#ff7d7d",
    breakpoint="#a445ee",
    normal="#f6c48f"
  ))+
  theme(
    panel.margin=grid::unit(0, "lines"),
    panel.border=element_rect(fill=NA, color="grey50")
  )+
  facet_grid(n.segments ~ .)+
  geom_vline(aes(
    xintercept=changepoint/1e6),
    data=zoom.changes,
    color="green",
    size=1,
    linetype="dashed")+
  geom_segment(aes(
    chromStart/1e6, mean,
    xend=chromEnd/1e6, yend=mean),
    data=zoom.segs,
    size=1,
    color="green")+
  penaltyLearning::geom_tallrect(aes(
    xmin=min/1e6,
    xmax=max/1e6,
    linetype=status),
    data=zoom.error.list$label.errors,
    fill=NA)+
  scale_linetype_manual("error type", values=c(
    correct=0,
    "false negative"=3,
    "false positive"=1))
zoom.selection <- penaltyLearning::modelSelection(
  zoom.models, complexity="n.segments")
zoom.error.join <- zoom.error.list$model.errors[J(zoom.selection), on=list(
  profile.id, chromosome, n.segments, loss)]
zoom.errors.tall <- data.table::melt(
  zoom.error.join,
  measure.vars=c("n.segments", "errors"))
zoom.target <- penaltyLearning::targetIntervals(
  zoom.error.join,
  problem.vars=c("profile.id", "chromosome"))
zoom.target.tall <- data.table::melt(
  zoom.target,
  measure.vars=c("min.log.lambda", "max.log.lambda"),
  variable.name="limit")[is.finite(value)]
ggplot()+
  geom_segment(aes(
    min.log.lambda, value,
    xend=max.log.lambda, yend=value),
    size=1,
    data=zoom.errors.tall)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(variable ~ ., scales="free")+
  scale_y_continuous("", breaks=0:max.segments)+
  xlab("log(penalty) = log(lambda)")+
  geom_text(aes(
    ifelse(limit=="min.log.lambda", value-1, value+1),
    errors,
    label=paste(
      "false", 
      ifelse(limit=="min.log.lambda", "positives", "negatives"),
      "\ntoo",
      ifelse(limit=="min.log.lambda", "many", "few"),
      "changes"),
    hjust=ifelse(
      limit=="min.log.lambda", 0, 1)),
    data=data.frame(zoom.target.tall, variable="errors"),
    vjust=-1)+
  geom_point(aes(
    value,
    errors,
    fill=limit),
    shape=21,
    size=4,
    data=data.frame(zoom.target.tall, variable="errors"))+
  scale_fill_manual("limit", values=c(
    min.log.lambda="black",
    max.log.lambda="white"))+
  theme(
    legend.position="bottom",
    legend.box="horizontal")

## Exercise 3: univariate regression using survreg. Try changing the
## profiles (some.ids), the input feature (feature.name), and the
## survreg distribution.
some.ids <- 145:155
feature.name <- "n.loglog" 
survreg.dist <- "gaussian" # try gaussian and logistic.
survreg.profiles <- someProfiles(neuroblastoma$profiles)
survreg.labels <- someProfiles(neuroblastoma$annotations)
feature.mat <- penaltyLearning::featureMatrix(
  survreg.profiles,
  problem.vars=c("profile.id", "chromosome"),
  data.var="logratio")
cat("Other possible features:\n")
print(colnames(feature.mat))
cat("Other possible profile ids:\n")
print(sort(unique(neuroblastoma$profiles$profile.id)))
survreg.features.dt <- survreg.profiles[, list(
  feature=feature.mat[paste(profile.id, chromosome), feature.name]
), by=list(profile.id, chromosome)]
setkey(survreg.profiles, profile.id, chromosome)
survreg.segs.list <- list()
survreg.selection.list <- list()
for(problem.i in 1:nrow(survreg.features.dt)){
  meta <- survreg.features.dt[problem.i,]
  pro <- survreg.profiles[meta]
  max.segments <- min(nrow(pro), 10)
  fit <- Segmentor3IsBack::Segmentor(
    pro$logratio, model=2, Kmax=max.segments)
  rss.vec <- rep(NA, max.segments)
  for(n.segments in 1:max.segments){
    end <- fit@breaks[n.segments, 1:n.segments]
    seg.mean.vec <- fit@parameters[n.segments, 1:n.segments]
    data.before.change <- end[-n.segments]
    data.after.change <- data.before.change+1
    pos.before.change <- as.integer(
      (pro$position[data.before.change]+pro$position[data.after.change])/2)
    start <- c(1, data.after.change)
    chromStart <- c(pro$position[1], pos.before.change)
    chromEnd <- c(pos.before.change, max(pro$position))
    data.mean.vec <- rep(seg.mean.vec, end-start+1)
    rss.vec[n.segments] <- sum((pro$logratio-data.mean.vec)^2)
    survreg.segs.list[[paste(problem.i, n.segments)]] <- data.table(
      meta,
      n.segments,
      start,
      end,
      chromStart,
      chromEnd,
      mean=seg.mean.vec)
  }
  loss.dt <- data.table(
    meta,
    n.segments=1:max.segments,
    loss=rss.vec)
  survreg.selection.list[[problem.i]] <- penaltyLearning::modelSelection(
    loss.dt, complexity="n.segments")
}
survreg.selection <- do.call(rbind, survreg.selection.list)
survreg.segs <- do.call(rbind, survreg.segs.list)
survreg.changes <- survreg.segs[1 < start]
survreg.error.list <- penaltyLearning::labelError(
  survreg.selection, survreg.labels, survreg.changes,
  problem.vars=c("profile.id", "chromosome"))
survreg.targets <- penaltyLearning::targetIntervals(
  survreg.error.list$model.errors,
  problem.vars=c("profile.id", "chromosome"))
survreg.targets.tall <- data.table::melt(
  survreg.targets,
  measure.vars=c("min.log.lambda", "max.log.lambda"),
  variable.name="limit",
  value.name="log.lambda")[is.finite(log.lambda)]
survreg.features.tall <- survreg.features.dt[survreg.targets.tall, on=list(
  profile.id, chromosome)]
library(survival)
survreg.features.targets <- survreg.features.dt[survreg.targets, on=list(
  profile.id, chromosome)]
fit.survreg <- survreg(
  Surv(min.log.lambda, max.log.lambda, type="interval2") ~ feature,
  survreg.features.targets, dist=survreg.dist)
survreg.line.dt <- data.table(t(coef(fit.survreg)), model="survreg")
survreg.pred <- data.table(
  survreg.features.dt,
  model="survreg",
  pred.log.lambda=predict(fit.survreg, survreg.features.dt))
survreg.res <- survreg.pred[survreg.features.targets, on=list(
  profile.id, chromosome)]
survreg.res[, residual := 0]
survreg.res[{
  is.finite(min.log.lambda)|is.finite(max.log.lambda)
}, residual := {
  penaltyLearning::targetIntervalResidual(
    cbind(min.log.lambda, max.log.lambda), pred.log.lambda)
}]
ggplot()+
  geom_point(aes(
    feature, log.lambda, fill=limit),
    shape=21,
    data=survreg.features.tall)+
  scale_fill_manual("limit", values=c(
    min.log.lambda="black",
    max.log.lambda="white"),
    breaks=c("max.log.lambda", "min.log.lambda"))+
  scale_x_continuous(
    paste("input feature:", feature.name))+
  scale_y_continuous(
    "output log(penalty)=log(lambda)")+
  geom_abline(aes(
    slope=feature, intercept=`(Intercept)`, color=model),
    size=1,
    data=survreg.line.dt)+
  geom_segment(aes(
    feature, pred.log.lambda,
    xend=feature, yend=pred.log.lambda-residual,
    color=model),
    data=survreg.res)

survreg.pred[, pred.log.penalty := pred.log.lambda]
show.selection <- data.table(survreg.selection)[survreg.pred, on=list(
  profile.id, chromosome,
  min.log.lambda < pred.log.lambda,
  max.log.lambda > pred.log.lambda)]
show.errors <- survreg.error.list$label.errors[{
  show.selection
}, on=list(
     profile.id, chromosome, n.segments), nomatch=0L]
show.changes <- survreg.changes[show.selection, on=list(
  profile.id, chromosome, n.segments), nomatch=0L]
show.segs <- survreg.segs[show.selection, on=list(
  profile.id, chromosome, n.segments)]
ggplot()+
  show.errors[, ggtitle(paste0(
    "Penalty learned with survreg on 1 feature has ",
    sum(fp+fn),
    " incorrect labels: ",
    sum(fp), " FP + ",
    sum(fn), " FN"))]+
  theme(
    panel.margin=grid::unit(0, "lines"),
    panel.border=element_rect(fill=NA, color="grey50"),
    legend.position="bottom",
    legend.box="horizontal"
  )+
  facet_grid(profile.id ~ chromosome, scales="free", space="free_x")+
  geom_point(aes(position/1e6, logratio),
             data=survreg.profiles,
             shape=1)+
  scale_x_continuous(
    "position on chromosome (mega bases)",
    breaks=c(100, 200))+
  scale_y_continuous(
    "logratio (approximate DNA copy number)",
    limits=c(-1,1)*1.1)+
  penaltyLearning::geom_tallrect(aes(
    xmin=min/1e6, xmax=max/1e6, fill=annotation),
    alpha=0.5,
    color=NA,
    data=survreg.labels)+
  scale_fill_manual("label", values=c(
    "1change"="#ff7d7d",
    breakpoint="#a445ee",
    normal="#f6c48f"
  ))+
  penaltyLearning::geom_tallrect(aes(
    xmin=min/1e6,
    xmax=max/1e6,
    linetype=status),
    fill=NA,
    data=show.errors)+
  scale_linetype_manual("error type", values=c(
    correct=0,
    "false negative"=3,
    "false positive"=1))+
  geom_vline(aes(
    xintercept=chromStart/1e6),
    data=show.changes,
    color="green",
    size=1,
    linetype="dashed")+
  geom_segment(aes(
    chromStart/1e6, mean,
    xend=chromEnd/1e6, yend=mean),
    data=show.segs,
    size=1,
    color="green")+
  theme(legend.box="horizontal")
