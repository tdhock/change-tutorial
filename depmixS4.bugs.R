works_with_R("3.4.0", depmixS4="1.3.4")

## Data set where nstates=1 does not work, Error in lm.wfit(x =
## as.matrix(object@x[!nas, ]), y = as.matrix(object@y[!nas, : missing
## or negative weights not allowed
if(!file.exists("buggy.data.RData")){
  download.file("http://members.cbio.ensmp.fr/~thocking/data/buggy.data.RData", "buggy.data.RData")
}
load("buggy.data.RData")
model.spec <- depmixS4::depmix(logratio ~ 1, data=buggy.data, nstates=1)
set.seed(1)
model.fit <- depmixS4::fit(model.spec)

## Data set where nstates=5 does not work, Error in fb(init = init, A
## = trDens, B = dens, ntimes = ntimes(object), : NA/NaN/Inf in
## foreign function call (arg 10)
if(!file.exists("buggy.5states.RData")){
  download.file("http://members.cbio.ensmp.fr/~thocking/data/buggy.5states.RData", "buggy.5states.RData")
}
load("buggy.5states.RData")
model.spec <- depmixS4::depmix(logratio ~ 1, data=buggy.5states, nstates=5)
set.seed(1)
model.fit <- depmixS4::fit(model.spec)
