source("packages.R")

load("Segmentor.models.RData")

data(neuroblastoma, package="neuroblastoma")

all.labels <- data.table(neuroblastoma$annotations)
setkey(all.labels, profile.id, chromosome, min, max)

all.changes <- Segmentor.models$segs[1 < start, ]
all.changes[, chromStart0 := chromStart]
setkey(all.changes, profile.id, chromosome, chromStart, chromStart0)

changes.in.labels <- foverlaps(all.changes, all.labels, nomatch=0L)
change.counts <- changes.in.labels[, list(
  changes=.N
  ), by=.(profile.id, chromosome, min, max, annotation, n.segments)]

all.loss <- Segmentor.models$loss
setkey(all.loss, profile.id, chromosome)
setkey(all.labels, profile.id, chromosome)
all.loss.labels <- all.loss[all.labels]
stopifnot(nrow(all.loss) == nrow(all.loss.labels))

setkey(all.loss.labels, profile.id, chromosome, min, max, annotation, n.segments)
all.loss.labels[change.counts, changes := change.counts$changes]
all.loss.labels[is.na(changes), changes := 0]
change.labels.dt <- data.table(change.labels)
setkey(change.labels.dt, annotation)
setkey(all.loss.labels, annotation)
label.min.max <- change.labels.dt[all.loss.labels]
stopifnot(nrow(all.loss.labels) == nrow(label.min.max))

label.min.max[, fp := ifelse(max.changes < changes, 1, 0)]
label.min.max[, fn := ifelse(changes < min.changes, 1, 0)]
all.errors <- label.min.max[, list(
  fp=sum(fp),
  possible.fp=sum(possible.fp),
  fn=sum(fn),
  possible.fn=sum(possible.fn),
  errors=sum(fp+fn),
  labels=.N
  ), by=.(profile.id, chromosome, n.segments, loss)]

