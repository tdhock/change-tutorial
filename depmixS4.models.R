source("packages.R")
works_with_R("3.4.0", depmixS4="1.3.4")

max.states <- 5

data(neuroblastoma)

problem.list <- with(neuroblastoma, {
  split(profiles, paste(profiles$profile.id, profiles$chromosome))
})
labeled.problem.names <- with(
  neuroblastoma$annotations, paste(profile.id, chromosome))
meta.df <- with(neuroblastoma$annotations, data.frame(profile.id, chromosome))
problem.i.vec <- seq_along(labeled.problem.names)
##problem.i.vec <- 1:10

change.list <- list()
loss.list <- list()
for(problem.i in problem.i.vec){
  problem.name <- labeled.problem.names[[problem.i]]
  meta <- meta.df[problem.i,]
  pid.chr <- paste(meta, collapse=".")
  cat(sprintf(
    "%4d / %4d problems %s\n",
    problem.i, length(labeled.problem.names), problem.name))
  pro <- problem.list[[problem.name]]
  ## TODO email i.visser@uva.nl with this data set, and ask how to
  ## specify the Gaussian homoscedastic model.
  if(FALSE){
    buggy.data <- pro
    save(buggy.data, file="buggy.data.RData")
  }
  for(nstates in 1:5){
    if(nstates==1){
      sd.est <- sqrt(var(pro$logratio)*(nrow(pro)-1)/nrow(pro))
      mean.est <- mean(pro$logratio)
      data.mean.vec <- rep(mean.est, nrow(pro))
      is.change <- FALSE
      neg.logLik <- -sum(dnorm(pro$logratio, mean.est, sd.est, log=TRUE))
    }else{
      set.seed(1)
      model.spec <- depmix(logratio ~ 1, data=pro, nstates=nstates)
      model.fit <- fit(model.spec, verbose=FALSE)
      is.change <- diff(model.fit@posterior$state) != 0
      i.change.vec <- which(is.change)
      change.list[[paste(pid.chr, nstates)]] <- data.table(
        meta,
        nstates,
        change=as.integer(
        (pro$position[i.change.vec]+pro$position[i.change.vec+1])/2))
      state.mean.vec <- sapply(model.fit@response, function(L){
        L[[1]]@parameters$coefficients
      })
      data.mean.vec <- state.mean.vec[model.fit@posterior$state]
      neg.logLik <- -as.numeric(logLik(model.fit))
    }
    loss.list[[paste(pid.chr, nstates)]] <- data.table(
      meta,
      nstates,
      nchanges=sum(is.change),
      rss=sum((pro$logratio-data.mean.vec)^2),
      neg.logLik)
  }
}
depmixS4.models <- list(
  change=do.call(rbind, change.list),
  loss=do.call(rbind, loss.list))

save(depmixS4.models, file="depmixS4.models.RData")
