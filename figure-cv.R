source("packages.R")

load("test.error.RData")

test.error[, accuracy.percent := 100-error.percent]
model.stats <- test.error[, list(
  mean=mean(error.percent),
  sd=sd(error.percent)
  ), by=.(model.name)][order(mean),]
test.error.tall <- melt(
  test.error,
  id.vars=c("test.fold", "model.name"),
  measure.vars=c("auc", "accuracy.percent"))
levs <- rev(model.stats$model.name)
test.error.tall[, model.fac := factor(model.name, levs)]
model.stats[, model.fac := factor(model.name, levs)]
model.stats[, variable := "accuracy.percent"]
best <- model.stats[which.min(mean),]
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(. ~ variable, scales="free")+
  geom_vline(aes(xintercept=100-mean), color="grey", data=best)+
  geom_point(aes(100-mean, model.fac), color="grey", size=4, data=model.stats)+
  geom_point(aes(value, model.fac), shape=1, data=test.error.tall)+
  xlab("")
png("figure-cv.png")
print(gg)
dev.off()
