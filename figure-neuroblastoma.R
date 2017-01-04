source("packages.R")

data(neuroblastoma)
profile.list <- with(neuroblastoma, split(profiles, profiles$profile.id))
label.list <- with(neuroblastoma, split(annotations, annotations$profile.id))

clinical.limited <- read.csv("clinical-limited.csv")
ids.str <- paste(clinical.limited$profile.id)
relapse.profile <- with(clinical.limited, paste(relapse, profile.id))
names(relapse.profile) <- ids.str

profiles <- do.call(rbind, profile.list[ids.str])
profiles$relapse.profile <- relapse.profile[paste(profiles$profile.id)]
labels <- do.call(rbind, label.list[ids.str])
labels$relapse.profile <- relapse.profile[paste(labels$profile.id)]

## Plot noisy data sets.
ggplot()+
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
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(relapse.profile ~ chromosome, scales="free", space="free_x")+
  geom_tallrect(aes(xmin=min/1e6, xmax=max/1e6, fill=annotation),
                data=labels)+
  scale_fill_manual(values=breakpoint.colors)+
  geom_point(aes(position/1e6, logratio),
             data=profiles,
             shape=1)+
  scale_x_continuous(
    "position on chromosome (mega bases)",
    breaks=c(100, 200))
