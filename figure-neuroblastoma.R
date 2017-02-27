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
  some
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

## Plot labels as well.
breakpoint.colors <- c(
  "breakpoint"="#a445ee",
  "normal"="#f6f4bf")
ggplot()+
  ggtitle("supervised change-point detection = data + labels")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(relapse.profile ~ chromosome, scales="free", space="free_x")+
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

one.pro <- subset(neuroblastoma$profiles, profile.id==4 & chromosome==14)
ggplot()+
  geom_point(aes(position/1e6, logratio),
             data=one.pro,
             shape=1)

label <- function(annotation, min, max){
  data.table(annotation, min, max)
}
one.pro.labels <- rbind(
  label("1change", 70e6, 80e6),
  label("0changes", 20e6, 60e6))
ggplot()+
  geom_tallrect(aes(xmin=min/1e6, xmax=max/1e6, fill=annotation),
                data=one.pro.labels)+
  scale_fill_manual("label", values=change.colors)+
  geom_point(aes(position/1e6, logratio),
             data=one.pro,
             shape=1)
