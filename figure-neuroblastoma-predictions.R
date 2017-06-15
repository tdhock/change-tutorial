source("packages.R")

data(neuroblastoma)
data(neuroblastomaProcessed)

profiles.dt <- data.table(neuroblastoma$profiles)
all.feature.mat <- penaltyLearning::featureMatrix(
  profiles.dt,
  problem.vars=c("profile.id", "chromosome"),
  data.var="logratio")

train.target.mat <- neuroblastomaProcessed$target.mat
train.pid.chr <- sub("[.]", " ", rownames(train.target.mat))
train.feature.mat <- all.feature.mat[train.pid.chr, ]
fit <- IntervalRegressionCV(train.feature.mat, train.target.mat)
mean.mat <- matrix(
  colMeans(train.feature.mat),
  nrow(all.feature.mat),
  ncol(all.feature.mat),
  byrow=TRUE)
pred.feature.mat <- all.feature.mat
is.bad <- !is.finite(all.feature.mat)
pred.feature.mat[is.bad] <- mean.mat[is.bad]
pred.log.pen.mat <- predict(fit, pred.feature.mat)

pred.change.dt <- profiles.dt[, {
  pid.chr <- paste0(profile.id, " ", chromosome)
  penalty <- exp(pred.log.pen.mat[pid.chr,])
  fit.fpop <- fpop::Fpop(logratio, penalty)
  end <- fit.fpop$t.est
  before.change <- end[-length(end)]
  after.change <- before.change+1L
  list(
    model="predicted",
    change=as.integer((position[before.change]+position[after.change])/2))
}, by=list(profile.id, chromosome)]

label.dt <- data.table(neuroblastoma$annotations)
error.list <- labelError(
  label.dt[, data.table(
    model="predicted", profile.id, chromosome)],
  label.dt, pred.change.dt,
  change.var="change",
  model.vars="model",
  problem.vars=c("profile.id", "chromosome"))

tit <- error.list$model.errors[, sprintf(
  paste(
    "Predicted IntervalRegressionCV penalty,",
    "%d/%d incorrect labels",
    "= %.1f%% train errors,",
    "%d/%d FP = %.1f%%,",
    "%d/%d FN = %.1f%%"),
  sum(errors), .N, sum(errors)/.N*100,
  sum(fp), sum(possible.fp), sum(fp)/sum(possible.fp)*100,
  sum(fn), sum(possible.fn), sum(fn)/sum(possible.fn)*100)]
print(tit)

show.pid.vec <- 1:10
show.pid.vec <- unique(profiles.dt$profile.id)
show.changes <- pred.change.dt[profile.id %in% show.pid.vec]
show.profiles <- profiles.dt[profile.id %in% show.pid.vec]
show.labels <- label.dt[profile.id %in% show.pid.vec]
show.errors <- error.list$label.errors[profile.id %in% show.pid.vec]
gg.supervised <- ggplot()+
  ggtitle(tit)+
  theme_bw()+
  theme(
    panel.margin=grid::unit(0, "lines"),
    legend.position="bottom",
    legend.box="horizontal"
  )+
  facet_grid(profile.id ~ chromosome, scales="free", space="free_x")+
  penaltyLearning::geom_tallrect(aes(
    xmin=min/1e6, xmax=max/1e6, fill=annotation),
    alpha=0.5,
    color=NA,
    data=show.labels)+
  geom_point(aes(position/1e6, logratio),
             data=show.profiles,
             color="grey40",
             shape=1)+
  scale_x_continuous(
    "position on chromosome (mega bases)",
    breaks=c(100, 200))+
  scale_y_continuous(
    "logratio (approximate DNA copy number)",
    breaks=c(-1, 0, 1),
    limits=c(-1,1)*1.3)+
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
    xintercept=change/1e6),
    data=show.changes,
    color="green",
    size=0.7,
    linetype="dashed")
n.profiles <- length(unique(show.profiles$profile.id))
png("figure-neuroblastoma-predictions.png", 1000, 110+30*n.profiles)
print(gg.supervised)
dev.off()
