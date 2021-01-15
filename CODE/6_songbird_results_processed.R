workdir <- "/home/acari/github/RECAST_project/"
setwd(workdir)

library(ggplot2)

# save data.frames for songbird analysis
sns_top10 <- read.csv("OUTPUT/song_sns_top10.csv")
sg_top10 <- read.csv("OUTPUT/song_sg_top10.csv")

sns_top10$group <- "settle/not settle"
sg_top10$group <- "stay/gone"

all_top10 <- rbind(sns_top10, sg_top10)

all_top10_plot <- ggplot(all_top10, aes(reorder(ko, effect_size), effect_size, fill = brite))+
     geom_bar(stat = "identity")+
     coord_flip()+
     theme_classic()+
     scale_fill_brewer(palette="Set1")+
     facet_wrap(~group, ncol = 2)+
     xlab("KO")

svg(filename="FIGURES/all_top10_plot.svg", width=8, height=6, pointsize=12)
all_top10_plot
dev.off()