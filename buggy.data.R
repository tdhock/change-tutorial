if(!file.exists("buggy.data.RData")){
  download.file("http://members.cbio.ensmp.fr/~thocking/data/buggy.data.RData", "buggy.data.RData")
}
load("buggy.data.RData")
library(depmixS4)
model.spec <- depmix(logratio ~ 1, data=buggy.data, nstates=1)
set.seed(1)
model.fit <- fit(model.spec)
