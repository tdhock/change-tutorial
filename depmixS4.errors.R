source("packages.R")

library(penaltyLearning)

load("depmixS4.models.RData")

data(neuroblastoma, package="neuroblastoma")

selection <- depmixS4.models$loss[, modelSelection(
  .SD, "neg.logLik", "nstates"
), by=.(profile.id, chromosome)]
some.anns <- subset(
  neuroblastoma$annotations,
  paste(profile.id, chromosome) %in% selection[, paste(profile.id, chromosome)])
errors <- labelError(
  selection, some.anns, depmixS4.models$change,
  change.var="change",
  label.vars=c("min", "max"),
  model.vars="nstates",
  problem.vars=c("profile.id", "chromosome"))

profiles <- data.table(neuroblastoma$profiles)
features <- profiles[paste(profile.id, chromosome) %in% selection[, paste(profile.id, chromosome)], list(
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

pred.dt <- pred.dot[threshold=="min.error" & feature.name=="log.log.n"]
pred.log.penalty <- pred.dt[, (min.thresh+max.thresh)/2]
selection[min.log.lambda < pred.log.penalty & pred.log.penalty < max.log.lambda]

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
  coord_equal(xlim=c(0, 0.5), ylim=c(0.5, 1))

