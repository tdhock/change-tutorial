source("packages.R")

library(penaltyLearning)

load("Segmentor.models.RData")

data(neuroblastoma, package="neuroblastoma")

selection <- Segmentor.models$loss[, modelSelection(
  .SD, complexity="n.segments"
  ), by=.(profile.id, chromosome)]
changes <- Segmentor.models$segs[1 < start, ]
errors <- labelError(
  selection, neuroblastoma$annotations, changes,
  change.var="chromStart",
  label.vars=c("min", "max"),
  model.vars="n.segments",
  problem.vars=c("profile.id", "chromosome"))

profiles <- data.table(neuroblastoma$profiles)
features <- profiles[, list(
  log.n=log(.N),
  log.log.n=log(log(.N))
  ), by=.(profile.id, chromosome)]

roc.list <- list()
auc.list <- list()
pred.dot.list <- list()
auc.polygon.list <- list()
for(feature.name in c("log.n", "log.log.n")){
  features$pred.log.lambda <- features[[feature.name]]
  feature.result <- ROChange(errors$model.errors, features, c("profile.id", "chromosome"))
  pred.dot.list[[feature.name]] <- data.table(feature.name, feature.result$thresholds)
  roc.list[[feature.name]] <- data.table(feature.name, feature.result$roc)
  auc.list[[feature.name]] <- data.table(feature.name, auc=feature.result$auc)
  auc.polygon.list[[feature.name]] <- data.table(feature.name, feature.result$auc.polygon)
}
roc <- do.call(rbind, roc.list)
auc <- do.call(rbind, auc.list)
pred.dot <- do.call(rbind, pred.dot.list)
auc.polygon <- do.call(rbind, auc.polygon.list)

ggplot()+
  theme_bw()+
  geom_segment(aes(min.thresh, errors, xend=max.thresh, yend=errors, color=feature.name),
               data=roc)+
  scale_fill_manual(values=c(predicted="grey", min.error="black"))+
  geom_point(aes((max.thresh+min.thresh)/2, errors, color=feature.name, fill=threshold),
             data=pred.dot,
             size=3,
             shape=21)

auc$TPR <- c(0.85, 0.8)
auc$FPR <- 0.15
ggplot()+
  theme_bw()+
  geom_text(aes(FPR, TPR, color=feature.name, label=sprintf("AUC = %.4f", auc)),
            data=auc)+
  geom_path(aes(FPR, TPR, color=feature.name, group=feature.name),
            size=2,
            alpha=0.25,
            data=roc)+
  geom_path(aes(FPR, TPR, color=feature.name, group=feature.name),
            size=1,
            alpha=0.25,
            data=auc.polygon)+
  scale_fill_manual(values=c(predicted="grey", min.error="black"))+
  geom_point(aes(FPR, TPR, color=feature.name, fill=threshold),
             data=pred.dot,
             size=3,
             shape=21)+
  coord_equal()

ggplot()+
  theme_bw()+
  geom_text(aes(FPR, TPR, color=feature.name, label=sprintf("AUC = %.4f", auc)),
            data=auc)+
  geom_path(aes(FPR, TPR, color=feature.name, group=feature.name),
            data=roc)+
  scale_fill_manual(values=c(default="grey", min.error="black"))+
  geom_point(aes(FPR, TPR, color=feature.name, fill=threshold),
             data=pred.dot,
             size=3,
             shape=21)+
  coord_equal()+
  xlim(0, 0.25)+ylim(0.75, 1)

