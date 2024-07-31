library(survival)
full.df <- structure(list(id = c(1L, 1L, 1L, 1L, 1L, 1L, 4L, 4L, 4L, 4L, 
4L, 4L, 6L, 6L, 6L, 6L, 6L, 6L, 8L, 8L, 8L, 8L, 8L, 8L, 10L, 
10L, 10L, 10L, 10L, 10L, 11L, 11L, 11L, 11L, 11L, 11L, 13L, 13L, 
13L, 13L, 13L, 13L, 120L, 120L, 120L, 120L, 120L, 120L), lo = c(0.195727051766757, 
-1.1087553605429, -1.00065642410924, -0.995123241844028, -Inf, 
-1.81065319423307, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -0.694428755856871, 
-2.4887060208859, -3.44756655828744, -3.03711252743929, -3.14523459850271, 
-0.892319064459033, -0.552651313503626, -Inf, -2.43732789951484, 
-1.24650632092422, -Inf, -Inf, -1.50791410241527, -2.24743277679978, 
-2.37206425088304, -3.83417364430955, -2.60448041546997, -1.44599430482559, 
-Inf, -0.456570392651195, -0.249593970642512, -Inf, 0.552800358347374, 
1.03693235764513, 0.366866754402124, -2.52677122267308, -2.48507949525124, 
-1.54283966348088, -1.59313273785538, -2.44940446306835, -0.513882229855861, 
-3.450873922835, -3.02209394810166, -3.97458092881313, -3.57901135468839, 
-2.1845818324618), hi = c(Inf, Inf, Inf, Inf, 2.06739808397873, 
Inf, 3.26518538624204, 1.13643272732406, 2.27284346410842, 1.38565083489553, 
0.492585941448176, 2.13713710441581, Inf, Inf, Inf, Inf, Inf, 
Inf, Inf, 0.937658965079102, Inf, Inf, 3.17164978103177, 3.0679625505633, 
Inf, Inf, Inf, Inf, Inf, Inf, 2.39473287695651, Inf, Inf, 0.730793396188341, 
Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 
Inf), feature = c(-2.65456939094638, -2.42126427679538, -2.0705540906376, 
-2.41236133374335, -2.38707136376868, -2.50995205248869, -2.36922164775071, 
-2.39136308372706, -2.35954108440544, -2.52855213967643, -2.61717785427761, 
-2.32148054113571, -2.90042209374967, -2.81341071676004, -2.90042209374967, 
-2.7646205525906, -2.97592964625781, -2.3434070875143, -2.87790396885268, 
-2.92818004477678, -2.91781583734955, -3.04532040079498, -2.44619735876083, 
-2.36685302923454, -3.17008566069877, -2.7646205525906, -2.54465665419377, 
-3.12356564506388, -3.20645330486964, -2.39689577246529, -2.43836426623075, 
-2.17187107950722, -2.52267471748792, -2.59336839058434, -1.84434871178509, 
-1.9783687469636, -3.01613697015902, -3.0756635490299, -2.93034936935621, 
-2.94970009845407, -2.89011456307748, -2.87110297166178, -3.01593498087151, 
-3.01593498087151, -2.94694210938456, -2.71055333132033, -3.11226608980994, 
-2.68824757380603)), class = "data.frame", row.names = c(NA, 
-48L))
model.fun.list <- list(
  survreg=function(train.df){
    fit <- with(train.df, survreg(
      Surv(lo, hi, type="interval2") ~ feature,
      dist="gaussian"))
    coef(fit)
  },
  penaltyLearning=function(train.df){
    fit <- with(train.df, penaltyLearning::IntervalRegressionUnregularized(
      cbind(feature), cbind(lo,hi)))
    coef(fit)[,"0"]
  })

dot.df.list <- list()
coef.df.list <- list()
text.df.list <- list()
for(extra.id in c(0,13,120)){
  sub.df <- subset(full.df, id %in% c(1, 4, 6, 8, 10, 11, extra.id))
  for(pkg in names(model.fun.list)){
    model.fun <- model.fun.list[[pkg]]
    try.out <- tryCatch({
      weight.vec <- model.fun(sub.df)
    }, warning=function(w)w$message)
    weight.vec <- model.fun(sub.df)
    dot.df.list[[paste(extra.id, pkg)]] <- with(sub.df, data.frame(
      extra.id, pkg,
      rbind(
        data.frame(limit="lo", output=lo, feature),
        data.frame(limit="hi", output=hi, feature))))
    coef.df.list[[paste(extra.id, pkg)]] <- data.frame(
      extra.id, pkg,
      warn.text=if(is.character(try.out))try.out else "",
      intercept=weight.vec[["(Intercept)"]],
      slope=weight.vec[["feature"]])
  }
}
dot.df <- subset(do.call(rbind, dot.df.list), is.finite(output))
coef.df <- do.call(rbind, coef.df.list)
text.df <- do.call(rbind, text.df.list)

library(ggplot2)
ggplot()+
  theme_bw()+
  geom_abline(aes(
    slope=slope, intercept=intercept),
    data=coef.df)+
  facet_grid(pkg ~ extra.id, labeller=label_both)+
  geom_point(aes(
    feature, output, fill=limit),
    shape=21,
    data=dot.df)+
  geom_text(aes(
    -3.2, 3, label=sprintf(
      "slope=%.2f\nintercept=%.2f\n%s", slope, intercept, warn.text)),
    vjust=1,
    hjust=0,
    data=coef.df)
