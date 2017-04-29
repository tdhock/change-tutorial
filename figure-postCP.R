source("packages.R")
requireGitHub::requireGitHub_package(
  "tdhock",
  "postCP_Improvement/postCP",
  "99f02e241bb66fd2c9cb80e4673c152d75f605da",
  "postCP")

source("packages.R")

## Read clinical data for six patients (relapse or ok, five years
## after treatment).
clinical.limited <- read.csv("clinical-limited.csv")
ids.str <- paste(clinical.limited$profile.id)
relapse.profile <- with(clinical.limited, paste(relapse, profile.id))
names(relapse.profile) <- ids.str

## Consider the subset of profiles and labels for these six patients.
someProfiles <- function(all.profiles){
  some <- subset(all.profiles, profile.id %in% ids.str)
  some$relapse.profile <- relapse.profile[paste(some$profile.id)]
  data.table(some)
}
data(neuroblastoma)
profiles <- someProfiles(neuroblastoma$profiles)
labels <- someProfiles(neuroblastoma$annotations)

## Plot noisy data sets.
ggplot()+
  ggtitle("unsupervised change-point detection = only noisy data series")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(relapse.profile ~ chromosome, scales="free", space="free_x")+
  geom_point(aes(position/1e6, logratio),
             data=profiles,
             shape=1)+
  scale_x_continuous(
    "position on chromosome (mega bases)",
    breaks=c(100, 200))

zoom.pid <- 4
zoom.chr <- 2
zoom.pro <- profiles[profile.id==zoom.pid & chromosome==zoom.chr]

max.segments <- 4
fit <- Segmentor3IsBack::Segmentor(zoom.pro$logratio, model=2, Kmax=max.segments)
zoom.segs.list <- list()
zoom.loss.vec <- rep(NA, max.segments)
for(n.segments in 1:max.segments){
  end <- fit@breaks[n.segments, 1:n.segments]
  data.before.change <- end[-n.segments]
  data.after.change <- data.before.change+1
  pos.before.change <- as.integer(
  (zoom.pro$position[data.before.change]+
   zoom.pro$position[data.after.change])/2)
  start <- c(1, data.after.change)
  chromStart <- c(zoom.pro$position[1], pos.before.change)
  chromEnd <- c(pos.before.change, max(zoom.pro$position))
  seg.mean.vec <- fit@parameters[n.segments, 1:n.segments]
  zoom.segs.list[[n.segments]] <- data.frame(
    profile.id=zoom.pid,
    chromosome=zoom.chr,
    n.segments,
    start,
    end,
    chromStart,
    chromEnd,
    mean=seg.mean.vec,
    row.names=NULL)
  data.mean.vec <- rep(seg.mean.vec, end-start+1)
  stopifnot(length(data.mean.vec)==nrow(zoom.pro))
  sq.res.vec <- (zoom.pro$logratio-data.mean.vec)^2
  zoom.loss.vec[n.segments] <- sum(sq.res.vec)
}

n.segs <- 4
var.est <- zoom.loss.vec[n.segs]/nrow(zoom.pro)
sd.est <- sqrt(var.est)
segs <- data.table(zoom.segs.list[[n.segs]])
changes <- segs[-.N, ]
fit <- postcp(logratio ~ 1, family=gaussian(), data=zoom.pro, bp=changes$end)
str(fit)
prob.dt <- with(fit, data.table(
  prob=as.numeric(post.cp),
  observation=as.integer(row(post.cp)),
  change=as.integer(col(post.cp))))
ggplot()+
  geom_point(aes(seq_along(logratio), logratio),data=data.table(zoom.pro, what="data"))+
  geom_line(aes(observation, prob, color=factor(change)), data=data.table(prob.dt, what="prob"))+
  geom_vline(aes(xintercept=end+0.5, color=factor(seq_along(end))), data=changes)

prob.dt[, prob.norm := prob/sum(prob), by=change]
ggplot()+
  geom_point(aes(seq_along(logratio), logratio),data=data.table(zoom.pro, what="data"))+
  geom_line(aes(observation, prob.norm, color=factor(change)), data=data.table(prob.dt, what="prob"))+
  geom_vline(aes(xintercept=end+0.5, color=factor(seq_along(end))), data=changes)+
  ylim(0, 0.1)

prob.dt[, list(total=sum(prob.norm)), by=change]
changes[, change := 1:.N]
join.dt <- prob.dt[changes, on=list(change)]
prob.dist <- join.dt[, list(
  total.prob=sum(prob.norm)
), by=list(change, dist.from.change=abs(end-observation))]
setkey(prob.dist, change, dist.from.change)
prob.dist[, cum.prob := cumsum(total.prob), by=change]
min.dist <- prob.dist[0.67 < cum.prob, .SD[1,], by=change]
error.bands <- changes[min.dist, on=list(change)]

gg <- ggplot()+
  geom_point(aes(seq_along(logratio), logratio),data=data.table(zoom.pro, what="data"))+
  penaltyLearning::geom_tallrect(aes(
    xmin=end-dist.from.change,
    xmax=end+dist.from.change,
    fill=factor(change)),
    alpha=0.5,
    color=NA,
    data=error.bands)+
  geom_vline(aes(xintercept=end+0.5, color=factor(seq_along(end))), data=changes)
print(gg)

pdf("figure-postCP.pdf")
print(gg)
dev.off()
