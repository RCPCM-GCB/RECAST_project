###################################################################################
# This R code designed to made Figure 3B "Density plot describes microbial reads  #
# distributions in settle-not settle and stay-gone “baskets” categories."         #
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
workdir <- "/home/acari/github/RECAST_project/" # please change to your work dorectory 
setwd(workdir)

# Set libraries
library(ggplot2)
library(tidyr)
library(stringr)
library(gridExtra)
library(reshape2)
library(pheatmap)

# Import tables
## Sorting
df_sorting <- read.csv("DATA/df_reads_count_sorting.org", sep = "\t", stringsAsFactors = F)
## Non-sorting
df_non_sorting <- read.csv("DATA/df_reads_count_non_sorting.org", sep = "\t", stringsAsFactors = F)

## Sample metadata
sample.metadata <- read.csv("DATA/sample.metadata", sep = "\t", stringsAsFactors = F)
sample.metadata.sbs <- sample.metadata[c(1:4,6:7)]
sample.metadata.sbs <- sample.metadata.sbs[sample.metadata.sbs$Status == "Allogenic",]

# Extract baseline IDs
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

# Add groups to sorting data frame 
groups <- sapply(str_split(df_sorting$sample, "\\_", n = 3), function(x) x[3])
pattern <- paste0(paste0("_", unique(groups)), collapse = "|")
non_sorting_sample_id <- gsub(pattern, "", df_sorting$sample)

group.df <- data.frame(non_sorting_sample_id, group = groups)
group.df$group <- as.character(group.df$group)
group.df$group[group.df$group == "s"] <- "settle"
group.df$group[group.df$group == "n"] <- "not_settle"
group.df$group <- gsub("came_", "", group.df$group)

df_sort_gr <- cbind(df_sorting, group.df)

# Make baseline reads count data frame 
df_baseline <- merge(base.meta, df_non_sorting[c(3,1,2)], by = 1)[-1]

# Make donor reads count data frame
df_donor <- sample.metadata.sbs[sample.metadata.sbs$Time !=0,][c(1,5)]
df_donor$Donor[df_donor$Donor == "FAT_DON_11"] <- "FAT_DON_11-22-pooled"
df_donor$Donor[df_donor$Donor == "FAT_DON_19"] <- "FAT_DON_19-22-0-0"
df_donor$Donor[df_donor$Donor == "FAT_DON_8"] <- "FAT_DON_8-22-0-0"
df_donor <- merge(df_donor[c(2,1)], df_non_sorting[c(3,1,2)], by = 1)
df_donor <- df_donor[-1]

# Merge Merge Merge
df_sort_gr_2 <- df_sort_gr[c(4,5,1,2)]

colnames(df_sort_gr_2)[1] <- "Sample"
colnames(df_baseline)[3] <- "n_reads_baseline"
colnames(df_donor)[3] <- "n_reads_donor"

df_sort_baseline <- merge(df_sort_gr_2, df_baseline, by = c("Sample", "taxa"), all = T)
df_sort_donor <- merge(df_sort_gr_2, df_donor, by = c("Sample", "taxa"), all = T)

df_sort_donor_sns <- df_sort_donor[df_sort_donor$group %in% c("settle", "not_settle", "from_donor"),]
df_sort_baseline_sg <- df_sort_baseline[df_sort_baseline$group %in% c("stay", "gone", "from_before"),]

df_sort_donor_sns_group <- df_sort_donor_sns[df_sort_donor_sns$group %in% c("settle", "not_settle"),]
df_sort_donor_sns_group.spread <- spread(df_sort_donor_sns_group, group, n_reads, fill = 0)
df_sort_donor_sns_group.spread$not_settle <- df_sort_donor_sns_group.spread$not_settle+1
df_sort_donor_sns_group.spread$settle <- df_sort_donor_sns_group.spread$settle+1

df_sort_donor_sns_group.spread$index <- df_sort_donor_sns_group.spread$settle/(df_sort_donor_sns_group.spread$settle+
                                                                                   df_sort_donor_sns_group.spread$not_settle)

df_sort_baseline_sg_group <- df_sort_baseline_sg[df_sort_baseline_sg$group %in% c("stay", "gone"),]
df_sort_baseline_sg_group.spread <- spread(df_sort_baseline_sg_group, group, n_reads, fill = 0)
df_sort_baseline_sg_group.spread$gone <- df_sort_baseline_sg_group.spread$gone+1
df_sort_baseline_sg_group.spread$stay <- df_sort_baseline_sg_group.spread$stay+1

df_sort_baseline_sg_group.spread$index <- df_sort_baseline_sg_group.spread$stay/(df_sort_baseline_sg_group.spread$stay+
                                                                                     df_sort_baseline_sg_group.spread$gone)

# Make density plot
density_plot <- ggplot()+
    geom_density(df_sort_donor_sns_group.spread, mapping = aes(index), col = "white", fill = 'red', alpha = 0.35)+
    geom_density(df_sort_baseline_sg_group.spread, mapping = aes(index), col = "white", fill = 'blue', alpha = 0.35)+
    theme_classic()+
    scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0), limits = c(-0.2, 1.2))+
    ggtitle("FMT")+
    ylim(c(0,10))

# Save plot
svg(filename="FIGURES/density_plot_case.svg", width=3.5, height=2.5, pointsize=12)
density_plot
dev.off()