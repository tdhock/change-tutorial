## Exercise: verify that BIC/SIC we computed is the same as penalty="SIC0" 
pid <- 8
chr <- 2
data.vec <- six.profiles[profile.id==pid & chromosome==chr, logratio]
(fit.pelt <- changepoint::cpt.mean(data.vec, penalty="SIC0", method="PELT"))
fit.fpop <- fpop::Fpop(data.vec, fit.pelt@pen.value)
str(fit.fpop)
six.BIC.dt[profile.id==pid & chromosome==chr, exp(pred.log.lambda)]
six.BIC.changes[profile.id==pid & chromosome==chr]

## Exercise: verify that PELT with a manual penalty computes the same
## model as we learned using survreg.
(one.selection <- six.survreg.selection[profile.id==pid & chromosome==chr])
six.profiles[profile.id==pid & chromosome==chr, changepoint::cpt.mean(
  logratio, penalty="Manual", method="PELT",
  pen.value=exp(one.selection$pred.log.penalty))]
six.survreg.changes[profile.id==pid & chromosome==chr]
