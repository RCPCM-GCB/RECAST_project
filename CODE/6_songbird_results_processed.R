###################################################################################
# This R code designed to made Figure 5B "The KOs differentially distinguish      # 
# “baskets” categories."                                                          #
# of Olekhnovich et al. (2020) manuscript                                         #
#                                                                                 #
# "Separation of donor and recipient microbial diversity                          #
# allow to determine taxonomic and functional features of                         #
# microbiota restructuring following fecal transplantation"                       #
#                                                                                 #
## E I. Olekhnovich, January 14, 2021                                             #
#                                                                                 #
###################################################################################

# Set work directory
workdir <- "/home/acari/github/RECAST_project/"
setwd(workdir)

# Set libraries
library(ggplot2)

# Read songbird analysis table 
sns_top10 <- read.csv("OUTPUT/song_sns_top10.csv")
sg_top10 <- read.csv("OUTPUT/song_sg_top10.csv")

# Processing tables
sns_top10$group <- "settle/not settle"
sg_top10$group <- "stay/gone"

# Merge Merge Merge
all_top10 <- rbind(sns_top10, sg_top10)

# Make plot
all_top10_plot <- ggplot(all_top10, aes(reorder(ko, effect_size), effect_size, fill = brite))+
     geom_bar(stat = "identity")+
     coord_flip()+
     theme_classic()+
     scale_fill_brewer(palette="Set1")+
     facet_wrap(~group, ncol = 2)+
     xlab("KO")

# Save plot
svg(filename="FIGURES/all_top10_plot.svg", width=8, height=6, pointsize=12)
all_top10_plot
dev.off()