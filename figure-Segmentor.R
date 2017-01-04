source("packages.R")

data(neuroblastoma)

pid <- "4"
chr <- "2"
pro <- subset(neuroblastoma$profiles, profile.id==pid & chromosome==chr)
lab <- subset(neuroblastoma$annotations, profile.id==pid & chromosome==chr)
if(lab$annotation=="breakpoint"){
  max.changes <- Inf 
  min.changes <- 1
}else{ # normal = no breakpoints.
  max.changes <- 0
  min.changes <- 0
}
breakpoint.colors <- c(
  "breakpoint"="#a445ee",
  "normal"="#f6f4bf")
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_tallrect(aes(xmin=min/1e6, xmax=max/1e6, fill=annotation),
                data=lab)+
  scale_fill_manual("label", values=breakpoint.colors)+
  geom_point(aes(position/1e6, logratio),
             data=pro,
             shape=1)+
  scale_x_continuous(
    "position on chromosome (mega bases)",
    breaks=c(100, 200))

max.segments <- 5
fit <- Segmentor3IsBack::Segmentor(pro$logratio, model=2, Kmax=max.segments)
segs.list <- list()
label.error.list <- list()
for(n.segments in 1:max.segments){
  end <- fit@breaks[n.segments, 1:n.segments]
  data.before.change <- end[-n.segments]
  data.after.change <- data.before.change+1
  pos.before.change <- as.integer(
    (pro$position[data.before.change]+pro$position[data.after.change])/2)
  changes.in.lab <- sum(lab$min < pos.before.change & pos.before.change < lab$max)
  status <- ifelse(
    changes.in.lab < min.changes, "false negative", ifelse(
      max.changes < changes.in.lab, "false positive", "correct"))
  label.error.list[[n.segments]] <- data.frame(n.segments, lab, status)
  start <- c(1, data.after.change)
  chromStart <- c(pro$position[1], pos.before.change)
  chromEnd <- c(pos.before.change, max(pro$position))
  segs.list[[n.segments]] <- data.frame(
    n.segments,
    start,
    end,
    chromStart,
    chromEnd,
    mean=fit@parameters[n.segments, 1:n.segments],
    row.names=NULL)
}
segs <- do.call(rbind, segs.list)
label.error <- do.call(rbind, label.error.list)

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_tallrect(aes(xmin=min/1e6, xmax=max/1e6, fill=annotation, linetype=status),
                size=1.5,
                data=label.error)+
  scale_fill_manual("label", values=breakpoint.colors)+
  geom_point(aes(position/1e6, logratio),
             data=pro,
             shape=1)+
  scale_x_continuous(
    "position on chromosome (mega bases)",
    breaks=c(100, 200))+
  scale_linetype_manual("error type",
                        values=c(correct=0,
                          "false negative"=3,
                          "false positive"=1))+
  facet_grid(n.segments ~ ., labeller=label_both)+
  geom_segment(
    aes(
      chromStart/1e6, mean,
      xend=chromEnd/1e6, yend=mean),
    data=segs,
    color="green")+
  geom_vline(aes(xintercept=chromStart/1e6),
    data=subset(segs, 1 < start),
    linetype="dashed",
    color="green")
