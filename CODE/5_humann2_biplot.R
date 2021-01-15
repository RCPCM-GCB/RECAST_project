###################################################################################
# This R code designed to made Figure 5A "Non-metric multidimensional scaling     #
# biplot obtained by using KOs distribution profiles of “baskets” categories      #
# and Bray-Curtis dissimilarity"                                                  #
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
library(stringr)
library(vegan)
library(ggplot2)
library(reshape2)

# Init KEGG orthology groups relative abundance cross samples data
df.ko <- read.csv("DATA/df.sorting.humann2.ko", sep = "\t", row.names = 1)

# Init group data
groups <- sapply(str_split(rownames(df.ko), "\\_", n = 3), function(x) x[3])
pattern <- paste0(paste0("_", unique(groups)), collapse = "|")
non_sorting_sample_id <- gsub(pattern, "", rownames(df.ko))

group.df <- data.frame(non_sorting_sample_id, sorting_id = rownames(df.ko), group = groups)
group.df$group <- as.character(group.df$group)
group.df$group[group.df$group == "s"] <- "settle"
group.df$group[group.df$group == "n"] <- "not_settle"
group.df$group <- gsub("came_", "", group.df$group)
group.df$sorting_id <- as.character(group.df$sorting_id)

ind_sbs <- as.numeric(rownames(group.df[!group.df$group %in% c("from_donor", "from_both", "from_before", "itself"),]))

# Init metadata
sample.metadata <- read.csv("DATA/sample.metadata", sep = "\t", stringsAsFactors = F)

group.metadata <- merge(group.df, sample.metadata, by = 1)[-1]
group.metadata <- group.metadata[-c(4,5)]
group.metadata <- group.metadata[c(1,3:6,2)]

# Filtering taxonomic data.frame
num.zeros <- colSums(df.ko==0)
save.zeros  <-  (num.zeros<(nrow(df.ko)*0.80))
df.ko.no_zeros  <-  df.ko[,(save.zeros==TRUE),]

# Calculate Bray-Curtis dissimilarity 
bray.mds <- metaMDS(df.ko.no_zeros[ind_sbs,])
bray.mds.points <- as.data.frame(bray.mds$points)

bray.mds.points$groups <- group.df[ind_sbs,]$group

# Make biplot
bray.mds.points$group <- factor(bray.mds.points$group, 
                                levels = c("settle", "not_settle", 
                                           "stay", "gone", "itself")
)

bray_biplot <- ggplot(bray.mds.points, aes(MDS1, MDS2, col = group))+
    geom_point(size = 1.8)+
    theme_bw()+
    # theme(legend.position = "bottom")+
    scale_color_brewer(palette="Set1")

# Save biplot
svg(filename="FIGURES/biplot_humann2.svg", width=4.5, height=3.3, pointsize=12)
bray_biplot
dev.off()