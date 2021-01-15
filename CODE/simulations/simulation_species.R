###################################################################################
# This R code designed to made Figure 2B "The results of the second step of the   #
# RECAST algorithm testing using the simulation datasets (species)" of            #
# Olekhnovich et al. (2020) manuscript                                            #
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
library(Metrics)
library(stringr)
library(reshape2)

# Make the function for calculate main classification quality metrics
metrics_calc <- function(df){
    y <- as.numeric(unlist(df.sbs["class"]))
    predictions <- as.numeric(unlist(df.sbs["scored.class"]))
    precision <- precision(predictions, y)
    recall <- recall(predictions, y)
    F1 <- (2 * precision * recall) / (precision + recall)
    DF <- data.frame(F1, precision, recall)
    return(DF)
}

# Import data
df.idxstats.all <- read.csv("DATA/simulation_species.txt", sep = "\t")

# Processing data
df.idxstats.all$coverage[is.na(df.idxstats.all$coverage)] <- 0
df.idxstats.all$mapped[is.na(df.idxstats.all$mapped)] <- 0
df.idxstats.all$coverage[df.idxstats.all$coverage > 0] <- 1
df.idxstats.all$mapped[df.idxstats.all$mapped > 0] <- 1

colnames(df.idxstats.all)[c(4,5)] <- c('scored.class', 'class')

df.idxstats.all$sample <- as.character(df.idxstats.all$sample)

df.f1.from_donor <- NULL
for (i in unique(df.idxstats.all$sample)){
    
    for (k in unique(df.idxstats.all$group)){
        
        df.sbs <- df.idxstats.all[df.idxstats.all$sample == i & df.idxstats.all$group == k,]
        df.sbs <- metrics_calc(df.sbs)
        df.sbs$sample <- i
        df.sbs$group <- k
        df.f1.from_donor <- rbind(df.sbs, df.f1.from_donor)    
    }
}

df.f1.from_donor$nreads <- as.numeric(gsub("high_|medium_|low_|M", "", df.f1.from_donor$sample))
df.f1.from_donor$complexity <- sapply(str_split(df.f1.from_donor$sample, "_"), function(x) x[1])
df.f1.from_donor <- melt(df.f1.from_donor, id.vars = c("sample", "nreads", "group", "complexity"))
colnames(df.f1.from_donor)[5] <- "metrics"
df.f1.from_donor$complexity <- factor(df.f1.from_donor$complexity, levels = c("high", "medium", "low"))

df.f1.from_donor.sbs <- df.f1.from_donor[df.f1.from_donor$metrics == "F1",]
df.f1.from_donor.sbs <- df.f1.from_donor.sbs[df.f1.from_donor.sbs$group != "came_from_both",]
df.f1.from_donor.sbs <- df.f1.from_donor.sbs[df.f1.from_donor.sbs$group != "came_itself",]

df.f1.from_donor.sbs$group <- factor(df.f1.from_donor.sbs$group, levels = c("came_from_donor", "came_from_before"))
species <- df.f1.from_donor.sbs

# Make the species classification quality plot
from_donor.metrics <- ggplot(df.f1.from_donor.sbs, aes(nreads, value, col = complexity))+
    geom_point()+
    facet_wrap(~group, ncol = 1)+
    stat_smooth()+
    ylim(c(0,1.05))+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_x_continuous(breaks = c(1,5,10,25,50))+
    # ggtitle("From donor")+
    xlab("# reads")+
    scale_color_brewer(palette = "Set1")

# Save plot
svg(filename="f1_metrics_species.svg", width=4.0, height=5.0)
from_donor.metrics
dev.off()