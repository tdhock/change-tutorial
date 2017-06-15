source("packages.R")

data(neuroblastoma)
data(neuroblastomaProcessed)

profiles.dt <- data.table(neuroblastoma$profiles)[profile.id==1]
##TODO compute full feature mat.

fit <- with(neuroblastomaProcessed, IntervalRegressionCV(feature.mat, target.mat))
pred.mat <- predict(fit, neuroblastomaProcessed$feature.mat)


profiles.dt[, {
  pid.chr <- paste0(profile.id, ".", chromosome)
  penalty <- pred.mat[pid.chr,]
  fit.fpop <- fpop::Fpop(logratio, penalty)
  end <- fit.fpop$t.est
  before.change <- end[-length(end)]
  after.change <- before.change+1L
  as.integer((position[before.change]+position[after.change])/2)
  }, by=list(profile.id, chromosome)]
