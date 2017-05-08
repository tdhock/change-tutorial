all.pred <- rbind(
  data.table(model="IntervalRegressionCV", six.cv.pred),
  data.table(model="BIC/SIC", six.BIC.dt, pred.log.penalty=NA),
  data.table(model="survreg", six.survreg.pred))
tidy <- all.pred[, {
  model.roc <- penaltyLearning::ROChange(
    six.error.list$model.errors, .SD[chromosome==11],
    problem.vars=c("profile.id", "chromosome"))
  roc.list[[model]] <- data.table(model, model.roc$roc)
  with(model.roc, data.table(
    roc=list(list(roc)),
    model, auc, thresholds[threshold=="predicted"]))
}, by=list(model)]

tall <- tidy[, roc[[1]], by=list(model)]
ggplot()+
  geom_path(aes(
    roc$FPR, roc$TPR, color=model),
    data=tall)
