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

plan(multiprocess)

depmixS4.models.list <- future_lapply(problem.i.vec, function(problem.i){
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
    buggy.5states <- pro
    save(buggy.5states, file="buggy.5states.RData")
  }
  prob.list.RData <- file.path("depmixS4.models", paste0(pid.chr, ".RData"))
  if(file.exists(prob.list.RData)){
    load(prob.list.RData)
  }else{
    prob.list <- future_lapply(1:5, function(nstates){
      logratio <- NULL
      prob.changes <- if(nstates==1){
        sd.est <- sqrt(var(pro$logratio)*(nrow(pro)-1)/nrow(pro))
        mean.est <- mean(pro$logratio)
        data.mean.vec <- rep(mean.est, nrow(pro))
        is.change <- FALSE
        neg.logLik <- -sum(dnorm(pro$logratio, mean.est, sd.est, log=TRUE))
        integer(0)
      }else{
        set.seed(1)
        model.spec <- depmix(logratio ~ 1, data=pro, nstates=nstates)
        model.fit <- fit(model.spec, verbose=FALSE)
        is.change <- diff(model.fit@posterior$state) != 0
        i.change.vec <- which(is.change)
        state.mean.vec <- sapply(model.fit@response, function(L){
          L[[1]]@parameters$coefficients
        })
        data.mean.vec <- state.mean.vec[model.fit@posterior$state]
        neg.logLik <- -as.numeric(logLik(model.fit))
        as.integer(
          (pro$position[i.change.vec]+pro$position[i.change.vec+1])/2)
      }
      data.table(
        meta,
        nstates,
        changes=list(list(prob.changes)),
        rss=sum((pro$logratio-data.mean.vec)^2),
        neg.logLik)
    })
    dir.create(dirname(prob.list.RData), showWarnings=FALSE, recursive=TRUE)
    save(prob.list, file=prob.list.RData)
  }
  do.call(rbind, prob.list)
})
depmixS4.models <- do.call(rbind, depmixS4.models.list)

save(depmixS4.models, file="depmixS4.models.RData")
