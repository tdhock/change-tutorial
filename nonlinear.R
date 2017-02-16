works_with_R(
  "3.3.2",
  "Rdatatable/data.table@9fadbcdfb7109b5d15adce07f826c9b58a86d2cf",
  party="1.0.13",
  partykit="2.0.2",#R CMD INSTALL partykit/pkg/devel/partykit/
  mlt="0.1.3",
  libcoin="0.9.1",
  trtf="0.1.1",#R CMD INSTALL ctm/pkg/trtf/
  rpart="4.1.10",
  randomForestSRC="2.4.1",
  interval="1.1.0.1",
  LTRCtrees="0.1.5",
  mboost="2.2.3",
  ## Above I have tried to use.
  ##bnnSurvival="1.0",
  ipred="0.9.5",
  LogicReg="1.5.9",
  kaps="1.0.2",
  Icens="1.46.0",
  MLEcens="0.1.4",
  npsurv="0.3.4",
  fitdistrplus="1.0.8",
  coarseDataTools="0.6.3",
  logconcens="0.16.4",
  "tdhock/penaltyLearning@82b16a3c204713818cc23517af85b3a92516af87")

ids.str <- paste(c(1, 4, 6, 8, 10, 11))
someProfiles <- function(all.profiles){
  data.table(all.profiles)[profile.id %in% ids.str, ]
}
data(neuroblastoma, package="neuroblastoma")
profiles <- someProfiles(neuroblastoma$profiles)
labels <- someProfiles(neuroblastoma$annotations)
## Plot labels along with noisy data sets.
breakpoint.colors <- c(
  "breakpoint"="#a445ee",
  "normal"="#f6f4bf")
library(ggplot2)
ggplot()+
  ggtitle("supervised change-point detection = data + labels")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(profile.id ~ chromosome, scales="free", space="free_x")+
  geom_tallrect(aes(xmin=min/1e6, xmax=max/1e6, fill=annotation),
                color="grey",
                data=labels)+
  scale_fill_manual("label", values=breakpoint.colors)+
  geom_point(aes(position/1e6, logratio),
             data=profiles,
             shape=1)+
  scale_x_continuous(
    "position on chromosome (mega bases)",
    breaks=c(100, 200))
problem.list <- split(profiles, profiles[, paste(profile.id, chromosome)])
segs.list <- list()
loss.list <- list()
for(problem.i in seq_along(problem.list)){
  problem.name <- names(problem.list)[[problem.i]]
  cat(sprintf(
    "%4d / %4d problems %s\n",
    problem.i, length(problem.list), problem.name))
  pro <- problem.list[[problem.name]]
  meta <- pro[1, .(profile.id, chromosome)]
  max.segments <- min(nrow(pro), 10)
  fit <- Segmentor3IsBack::Segmentor(
    pro$logratio, model=2, Kmax=max.segments)
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
    segs.list[[paste(problem.name, n.segments)]] <- data.table(
      meta,
      n.segments,
      start,
      end,
      chromStart,
      chromEnd,
      mean=seg.mean.vec)
  }
  loss.list[[paste(problem.name, n.segments)]] <- data.table(
    meta,
    n.segments=1:max.segments,
    loss=as.numeric(fit@likelihood))
}
loss <- do.call(rbind, loss.list)
segs <- do.call(rbind, segs.list)
selection <- loss[, {
  penaltyLearning::modelSelection(.SD, "loss", "n.segments")
}, by=.(profile.id, chromosome)]
changes <- segs[1 < start, ]
errors <- penaltyLearning::labelError(
  selection, labels, changes,
  change.var="chromStart",
  label.vars=c("min", "max"),
  problem.vars=c("profile.id", "chromosome"))
all.errors <- data.table(errors$model.errors)
all.errors[, set := ifelse(chromosome=="11", "test", "train")]
target.dt <- targetIntervals(
  all.errors[set=="train", ],
  c("profile.id", "chromosome"))
target.mat <- target.dt[, cbind(min.log.lambda, max.log.lambda)]
rownames(target.mat) <- target.dt[, paste(profile.id, chromosome)]
feature.dt <- profiles[, list(
  log.data=log(.N),
  log.var=log(median(abs(diff(logratio))))
  ), by=.(profile.id, chromosome)]
all.feature.mat <- feature.dt[, cbind(log.data, log.var)]
rownames(all.feature.mat) <- feature.dt[, paste(profile.id, chromosome)]
train.feature.mat <- all.feature.mat[rownames(target.mat), ]

## Try using my data set.
finite.targets <- data.frame(log.lambda=target.mat[is.finite(target.mat)])
m <- ctm(as.basis(~log.lambda, data=finite.targets), todistr="Normal")
train.Surv <- target.dt[, Surv(min.log.lambda, max.log.lambda, type="interval2")]
train.df <- data.frame(log.lambda=train.Surv, train.feature.mat)
mlt.fit <- mlt(m, data=train.df)
tree.fit <- trafotree(
  m, formula = log.lambda ~ ., data=train.df,
  mltargs=list(theta=coef(mlt.fit)))

## penaltyLearning
fit <- IntervalRegressionUnregularized(train.feature.mat, target.mat)

## party package. ctree docs mention "median predicted survival times"
## and give an example using Surv(time, cens) but no mention of
## interval censoring. cforest example also says Surv(time,
## cens). References mention "Survival ensembles" by Hothorn et al
## https://www.ncbi.nlm.nih.gov/pubmed/16344280. Both ctree and
## cforest learn a constant function.
tfit <- ctree(Surv(min.log.lambda, max.log.lambda, type="interval2") ~ log.data + log.var, data=data.frame(target.dt, train.feature.mat))
predict(tfit)
ffit <- cforest(Surv(min.log.lambda, max.log.lambda, type="interval2") ~ log.data + log.var, data=data.frame(target.dt, train.feature.mat))
predict(ffit)
tfit <- ctree(Surv(exp(min.log.lambda), exp(max.log.lambda), type="interval2") ~ log.data + log.var, data=data.frame(target.dt, train.feature.mat))
predict(tfit)
ffit <- cforest(Surv(exp(min.log.lambda), exp(max.log.lambda), type="interval2") ~ log.data + log.var, data=data.frame(target.dt, train.feature.mat))
predict(ffit)

## rpart documentation mentions 'If ‘y’ is a survival object, then
## ‘method = "exp"’ is assumed' but does not give example. Training
## stops with an error (Observation time must be > 0)
pfit <- rpart(Surv(min.log.lambda, max.log.lambda, type="interval2") ~ log.data + log.var, data=data.frame(target.dt, train.feature.mat))
pfit <- rpart(Surv(exp(min.log.lambda), exp(max.log.lambda), type="interval2") ~ log.data + log.var, data=data.frame(target.dt, train.feature.mat))

## randomForestSRC::rfsrc mentions "right-censored (including
## competing risk)" and examples include Surv(days, status), but I got
## an error "for survival families censoring variable must be coded as
## a non-negative integer (perhaps the formula is set incorrectly?)"
rfit <- rfsrc(Surv(min.log.lambda, max.log.lambda, type="interval2") ~ log.data + log.var, data=data.frame(target.dt, train.feature.mat))
rfit <- rfsrc(Surv(exp(min.log.lambda), exp(max.log.lambda), type="interval2") ~ log.data + log.var, data=data.frame(target.dt, train.feature.mat))

## interval::icfit does not have a predict method, and the input
## variables are treated as categorical.
icout <- icfit(Surv(exp(min.log.lambda), exp(max.log.lambda), type="interval2") ~ log.data + log.var, data=data.frame(target.dt, train.feature.mat))
predict(icout)
coef(icout)

## LTRCtrees "fit left-truncated and right censored (LTRC) data."
## Always gives errors or predicts Inf. Examples do not include open
## intervals.
cart.fit <- LTRCART(Surv(min.log.lambda, max.log.lambda, type="interval2") ~ log.data + log.var, data.frame(target.dt, train.feature.mat))
cart.fit <- LTRCART(Surv(min.log.lambda, max.log.lambda, 0) ~ log.data + log.var, data.frame(target.dt, train.feature.mat))
cart.fit <- LTRCART(Surv(min.log.lambda, max.log.lambda, 1) ~ log.data + log.var, data.frame(target.dt, train.feature.mat))
cit.fit <- LTRCIT(Surv(min.log.lambda, max.log.lambda, type="interval2") ~ log.data + log.var, data.frame(target.dt, train.feature.mat))
cit.fit <- LTRCIT(Surv(min.log.lambda, max.log.lambda, rep(0, length(min.log.lambda))) ~ log.data + log.var, data.frame(target.dt, train.feature.mat))
cit.fit <- LTRCIT(Surv(min.log.lambda, max.log.lambda, rep(1, length(min.log.lambda))) ~ log.data + log.var, data.frame(target.dt, train.feature.mat))

## mboost
ipc <- target.dt[, IPCweights(Surv(min.log.lambda, max.log.lambda, type="interval2"))]

## trtf code from Torsten Hothorn.
set.seed(29)
### unconditional example: left / right-censoring only
n <- 10000
y <- rnorm(n, mean = 3, sd = 1)   ### response
i <- rnorm(n, mean = 3, sd = .5)  ### censoring
left <- y < i
d <- data.frame(y = y)
### set-up linear transformation: P(Y <= y) = pnorm(1 * y - 3)
m <- ctm(as.basis(~ y, data = d), todistr = "Normal")
### set-up left- and right-censored observations
### use mlt::R here for simpler starting value computations
d$y <- R(cleft = ifelse(left, -Inf, i), cright = ifelse(left, i, Inf))
mod <- mlt(m, data = d)
coef(mod) ### c(-3, 1) are the "true" values
### conditional example
x <- runif(n)
### two normals, one with mean 2.5 and one with mean 3.5
y <- rnorm(n, mean = c(2.5, 3.5)[(x > .5) + 1], sd = 1)
i <- rnorm(n, mean = 3, sd = .5) ### censoring
left <- y < i
d <- data.frame(y = y, x = x)
### same model
m <- ctm(as.basis(~ y, data = d), todistr = "Normal")
### use Surv here because trafotree can't deal with R() 
### at the moment
d$y <- Surv(time = ifelse(left, -Inf, i), time2 = ifelse(left, i, Inf),
            type = "interval2")
### need to supply starting value because Surv can't
### handle this at the moment
coef(mod)
mod <- mlt(m, data = d, theta = coef(mod))
coef(mod)
### seach for x > .5 split, make sure to supply starting value
tr <- trafotree(m, formula = y ~ x, data = d, 
                mltargs = list(theta = coef(mod)),
                control = ctree_control(stump = TRUE))
### print(tr) is buggy, here is the tree
node_party(tr)
### and here are the coefficients in the two daughther 
### nodes: (-2.5, 1) and (-3.5, 1)
tr$coef
### plot the tree
plot(tr, tp_args = list(type = "density", K = 200))


feature.dt$pred.log.lambda <- fit$predict(all.feature.mat)
test.pred <- feature.dt[chromosome=="11",]
ROChange(all.errors, test.pred, c("profile.id", "chromosome"))
feature.dt[, pred.log.penalty := pred.log.lambda ]
setkey(feature.dt, profile.id, chromosome, pred.log.lambda, pred.log.penalty)
setkey(selection, profile.id, chromosome, min.log.lambda, max.log.lambda)
pred.models <- foverlaps(feature.dt, selection)
setkey(pred.models, profile.id, chromosome, n.segments)
setkey(segs, profile.id, chromosome, n.segments)
pred.segs <- segs[pred.models]
pred.changes <- pred.segs[1 < start, ]
setkey(errors$label.errors, profile.id, chromosome, n.segments)
pred.labels <- errors$label.errors[pred.models, nomatch=0L]
ggplot()+
  ggtitle("data + labels + predicted segment means and changes")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(profile.id ~ chromosome, scales="free", space="free_x")+
  geom_tallrect(aes(
    xmin=min/1e6, xmax=max/1e6, fill=annotation, linetype=status),
                size=1.5,
                data=pred.labels)+
  scale_linetype_manual("error type",
                        values=c(correct=0,
                          "false negative"=3,
                          "false positive"=1))+
  scale_fill_manual("label", values=breakpoint.colors)+
  geom_point(aes(position/1e6, logratio),
             data=profiles,
             shape=1)+
  scale_x_continuous(
    "position on chromosome (mega bases)",
    breaks=c(100, 200))+
  geom_segment(aes(chromStart/1e6, mean, xend=chromEnd/1e6, yend=mean),
               data=pred.segs,
               color="green")+
  geom_vline(aes(xintercept=chromStart/1e6),
             data=pred.changes,
             color="green",
             linetype="dashed")

