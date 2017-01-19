## multivariate model.
cvfit <- IntervalRegressionCV(all.features.mat, target.mat, verbose=1)

cv.dt <- data.table(train.dt)
cv.dt[, pred.log.lambda := cvfit$predict(all.features.mat)]

## also plot ROC curves.
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

model.list <- list(
  BIC=train.dt,
  learned=pred.dt,
  cv=cv.dt)
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

ggplot()+
  theme_bw()+
  ylab("total incorrect labels")+
  xlab("constant/threshold added to log(penalty)")+
  geom_segment(aes(min.thresh, errors, xend=max.thresh, yend=errors, color=model.name),
               data=roc)+
  scale_fill_manual(values=c(predicted="grey", min.error="black"))+
  geom_point(aes((max.thresh+min.thresh)/2, errors, color=model.name, fill=threshold),
             data=pred.dot,
             size=3,
             shape=21)

auc$TPR <- c(0.8, 0.85, 0.9)
auc$FPR <- 0.15
ggplot()+
  theme_bw()+
  geom_text(aes(FPR, TPR, color=model.name, label=sprintf("AUC = %.4f", auc)),
            data=auc)+
  geom_path(aes(FPR, TPR, color=model.name, group=model.name),
            size=2,
            alpha=0.25,
            data=roc)+
  geom_path(aes(FPR, TPR, color=model.name, group=model.name),
            size=1,
            alpha=0.25,
            data=auc.polygon)+
  scale_fill_manual(values=c(predicted="grey", min.error="black"))+
  geom_point(aes(FPR, TPR, color=model.name, fill=threshold),
             data=pred.dot,
             size=3,
             shape=21)+
  coord_equal()

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
  coord_equal()+
  xlim(0, 0.25)+ylim(0.75, 1)

