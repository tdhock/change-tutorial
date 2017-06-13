works_with_R(
  c("3.3.3", "3.4.0"),
  neuroblastoma="1.0",
  future="1.4.0",
  Segmentor3IsBack="2.0",
  changepoint="2.2",
  directlabels="2017.3.31",
  data.table=c("1.10.4", "1.10.5"),
  survival=c("2.41.2", "2.41.3"),
  ggplot2=c("2.1.0", "2.2.1"),
  penaltyLearning="2017.5.8")
future::plan(multiprocess)

## Exercise 1: compute label error and plot predicted changepoints for
## another set of profiles and unsupervised penalty.
some.ids <- 100:105
pen.name <- "SIC0"
data(neuroblastoma, package="neuroblastoma")
someProfiles <- function(all.profiles){
  data.table(all.profiles)[paste(profile.id) %in% some.ids]
}
six.profiles <- someProfiles(neuroblastoma$profiles)
six.labels <- someProfiles(neuroblastoma$annotations)
unsupervised.models <- six.profiles[, {
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
  unsupervised.models, six.labels, unsupervised.changes,
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
             data=six.profiles,
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
    data=six.labels)+
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

## Exercise 2: label a profile.
zoom.profile <- 104     #Exercise: change this value!
zoom.chromosome <- 14 #Exercise: change this value!
label <- function(min, max, annotation){
  data.frame(
    profile.id=paste(zoom.profile),
    chromosome=paste(zoom.chromosome),
    min, max, annotation)
}
zoom.labels <- rbind(# Exercise: change these values!
  label(70e6, 90e6, "1change"), 
  label(20e6, 60e6, "normal"))
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

## Exercise 3: compute a sequence of maximum likelihood segmentation
## models.
max.segments <- 10
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

