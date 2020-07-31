workdir <- "/home/acari/github/RECAST_project/"

library(ggplot2)
library(tidyr)
library(stringr)
library(gridExtra)
library(reshape2)
library(pheatmap)

setwd(workdir)

# import data tables
df_sorting <- read.csv("DATA/df_reads_count_benchmark_sorting.org", sep = "\t", stringsAsFactors = F)
df_non_sorting <- read.csv("DATA/df_reads_count_benchmark_non_sorting.org", sep = "\t", stringsAsFactors = F)

sample.metadata <- read.csv("DATA/sample.metadata", sep = "\t", stringsAsFactors = F)
sample.metadata.sbs <- sample.metadata[c(1:3,7)]
sample.metadata.sbs <- sample.metadata.sbs[sample.metadata.sbs$Dataset == "VOIGT15",]

base.meta <- NULL
for (subject in unique(sample.metadata.sbs$Subject)){
    
    dSBS <- sample.metadata.sbs[sample.metadata.sbs$Subject == subject,]
    
    base.meta <- rbind(
        data.frame(
            baseline = dSBS$Sample[dSBS$Time == 0],
            Sample = dSBS$Sample[dSBS$Time != 0]
            
        ), 
        base.meta
    )
    
}

# add groups to sorting data frame 
samples <- sapply(str_split(df_sorting$sample, "\\_", n = 3), function(x) x[1])
donors <- sapply(str_split(df_sorting$sample, "\\_", n = 3), function(x) x[2])
groups <- sapply(str_split(df_sorting$sample, "\\_", n = 4), function(x) x[4])

group.df <- data.frame(samples, donor = donors, group = groups)

group.df$group <- as.character(group.df$group)
group.df$group[group.df$group == "s"] <- "settle"
group.df$group[group.df$group == "n"] <- "not_settle"
group.df$group <- gsub("came_", "", group.df$group)

df_sort_gr <- cbind(group.df, df_sorting[-3])

df_sns <- df_sort_gr[df_sort_gr$group %in% c("settle", "not_settle"),]
df_sns_spread <- spread(df_sns, group, n_reads, fill = 0)
df_sns_spread$not_settle <- df_sns_spread$not_settle+1
df_sns_spread$settle <- df_sns_spread$settle+1
df_sns_spread$index <- df_sns_spread$settle/(df_sns_spread$settle+df_sns_spread$not_settle)

df_sg <- df_sort_gr[df_sort_gr$group %in% c("stay", "gone"),]
df_sg_spread <- spread(df_sg, group, n_reads, fill = 0)
df_sg_spread$gone <- df_sg_spread$gone+1
df_sg_spread$stay <- df_sg_spread$stay+1
df_sg_spread$index <- df_sg_spread$stay/(df_sg_spread$stay+df_sg_spread$gone)

density_plot <- ggplot()+
    geom_density(df_sns_spread, mapping = aes(index), col = "white", fill = 'red', alpha = 0.35)+
    geom_density(df_sg_spread, mapping = aes(index), col = "white", fill = 'blue', alpha = 0.35)+
    theme_classic()+
    scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0), limits = c(-0.2, 1.2))+
    ggtitle("Control")
    
svg(filename="FIGURES/benchmark_density_plot.svg", width=3.5, height=2.5, pointsize=12)
density_plot
dev.off()