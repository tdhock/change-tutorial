source("packages.R")

max.segments <- 20

data(neuroblastoma)

problem.list <- with(neuroblastoma, {
  split(profiles, paste(profiles$profile.id, profiles$chromosome))
})
labeled.problem.names <- with(neuroblastoma$annotations, paste(profile.id, chromosome))
meta.df <- with(neuroblastoma$annotations, data.frame(profile.id, chromosome))

segs.list <- list()
loss.list <- list()
for(problem.i in seq_along(labeled.problem.names)){
  problem.name <- labeled.problem.names[[problem.i]]
  meta <- meta.df[problem.i,]
  cat(sprintf("%4d / %4d problems %s\n", problem.i, length(labeled.problem.names), problem.name))
  pro <- problem.list[[problem.name]]
  fit <- Segmentor3IsBack::Segmentor(pro$logratio, model=2, Kmax=max.segments)
  rss.vec <- rep(NA, max.segments)
  for(n.segments in 1:max.segments){
    end <- fit@breaks[n.segments, 1:n.segments]
    data.before.change <- end[-n.segments]
    data.after.change <- data.before.change+1
    pos.before.change <- as.integer(
      (pro$position[data.before.change]+pro$position[data.after.change])/2)
    start <- c(1, data.after.change)
    chromStart <- c(pro$position[1], pos.before.change)
    chromEnd <- c(pos.before.change, max(pro$position))
    seg.mean.vec <- fit@parameters[n.segments, 1:n.segments]
    data.mean.vec <- rep(seg.mean.vec, end-start+1)
    residual.vec <- pro$logratio - data.mean.vec
    rss.vec[n.segments] <- sum(residual.vec * residual.vec)
    segs.list[[paste(problem.name, n.segments)]] <- data.table(
      meta,
      n.segments,
      start,
      end,
      chromStart,
      chromEnd,
      mean=seg.mean.vec)
  }
  if(FALSE){
    ## The likelihood computed by the Segmentor function is just an
    ## affine transformation of the sum of squared residuals.
    plot(rss.vec, fit@likelihood)
  }
  loss.list[[paste(problem.name, n.segments)]] <- data.table(
    meta,
    n.segments=1:max.segments,
    loss=rss.vec)
}

Segmentor.models <- list(
  segs=do.call(rbind, segs.list),
  loss=do.call(rbind, loss.list))

save(Segmentor.models, file="Segmentor.models.RData")
