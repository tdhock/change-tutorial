source("packages.R")

load("Segmentor.models.RData")
setkey(Segmentor.models$loss, profile.id, chromosome)
setkey(Segmentor.models$segs, profile.id, chromosome, n.segments)

data(neuroblastoma)

problem.list <- with(neuroblastoma, {
  split(profiles, paste(profiles$profile.id, profiles$chromosome))
})
labeled.problem.names <- with(neuroblastoma$annotations, paste(profile.id, chromosome))
meta.df <- with(neuroblastoma$annotations, data.frame(profile.id, chromosome))

for(problem.i in seq_along(labeled.problem.names)){
  problem.name <- labeled.problem.names[[problem.i]]
  meta <- meta.df[problem.i,]
  cat(sprintf(
    "%4d / %4d problems %s\n",
    problem.i, length(labeled.problem.names), problem.name))
  pro <- problem.list[[problem.name]]
  fit.pelt <- changepoint::cpt.mean(pro$logratio, penalty="SIC0", method="PELT")
  loss <- Segmentor.models$loss[meta]
  selection <- modelSelection(loss, complexity="n.segments")
  penalty <- log(nrow(pro))
  fit.fpop <- fpop::Fpop(pro$logratio, penalty)
  stopifnot(all.equal(penalty, fit.pelt@pen.value))
  n.segs <- data.table(
    selection)[min.lambda < penalty & penalty < max.lambda, n.segments]
  stopifnot(identical(length(fit.pelt@cpts), n.segs))
  segs <- Segmentor.models$segs[J(meta, n.segments=n.segs)]
  stopifnot(identical(segs$end, fit.pelt@cpts))
  stopifnot(identical(segs$end, fit.fpop$t.est))
}

pelt.models <- "SIC0 from changepoint and fpop is the same as our BIC"
save(pelt.models, file="pelt.models.RData")
