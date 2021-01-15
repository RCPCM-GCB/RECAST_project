###################################################################################
# This R code designed to made Figure 3A "Non-metric multidimensional scaling     #
# biplot obtained using species distribution profiles of “baskets” categories     #
# and Bray-Curtis dissimilarity" of Olekhnovich et al. (2020) manuscript          #
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
library(stringr)
library(zCompositions)
library(pheatmap)
library(vegan)
library(ggplot2)

# Import tables
## Case
df.org <- read.table("DATA/df.sorting.mp2.org", sep = "\t", stringsAsFactors = F, row.names = NULL)
colnames(df.org)[1] <- "Samples"
colnames(df.org) <- gsub("\\.", "\\|", colnames(df.org))

## Controls
df.org.benchmark <- read.csv("DATA/df.sorting_benchmark.mp2.org", sep = "\t", stringsAsFactors = F, row.names = NULL)
colnames(df.org.benchmark)[1] <- "Samples"
colnames(df.org.benchmark) <- gsub("\\.", "\\|", colnames(df.org.benchmark))

# Merge Case and Control tables to taxonomic data.frame
df.org.all <- merge(df.org, df.org.benchmark, all = TRUE)
df.org.all[is.na(df.org.all)] <- 0
df.org.all <- df.org.all[-1]

# Make Group.df metadata table
groups <- sapply(str_split(df.org$Samples, "\\_", n = 3), function(x) x[3])
pattern <- paste0("_",unique(groups), collapse = "|")
non_sorting_sample_id <- gsub(pattern, "", df.org$Samples)

group.df.fmt <- data.frame(
    Sample = df.org$Samples, 
    non_sorting_sample_id, 
    group = groups, 
    dataset = "FMT"
)

groups.benchmark <- sapply(str_split(df.org.benchmark$Samples, "\\_", n = 4), function(x) x[4])
pattern.benchmark <- paste0(paste0("_donor", "_", unique(groups.benchmark)), collapse = "|")
non_sorting_sample_id_benchmark <- gsub(pattern.benchmark, "", df.org.benchmark$Samples)
non_sorting_sample_id_benchmark <- gsub("_daisy|_bugkiller|_scavenger|_peacemaker|_tigress", "", non_sorting_sample_id_benchmark)

group.df.benchmark <- data.frame(
    Sample = df.org.benchmark$Samples, 
    non_sorting_sample_id = non_sorting_sample_id_benchmark, 
    group = groups.benchmark, 
    dataset = "Control"
)

group.df <- rbind(group.df.fmt, group.df.benchmark)

group.df$group <- as.character(group.df$group)
group.df$group[group.df$group == "s"] <- "settle"
group.df$group[group.df$group == "n"] <- "not_settle"
group.df$group <- gsub("came_", "", group.df$group)

ind_sbs <- rownames(group.df[!group.df$group %in% c("from_donor", "from_both", "from_before", "itself"),])

# Make taxonomy metadata tables
colnames(df.org.all) <- gsub("k__|p__|c__|o__|f__|g__|s__", "", colnames(df.org.all))
tax.kingdom <- sapply(str_split(colnames(df.org.all), "\\|"), function(x) x[1])
tax.phyla <- sapply(str_split(colnames(df.org.all), "\\|"), function(x) x[2])
tax.class <- sapply(str_split(colnames(df.org.all), "\\|"), function(x) x[3])
tax.order <- sapply(str_split(colnames(df.org.all), "\\|"), function(x) x[4])
tax.family <- sapply(str_split(colnames(df.org.all), "\\|"), function(x) x[5])
tax.genus <- sapply(str_split(colnames(df.org.all), "\\|"), function(x) x[6])
tax.species <- sapply(str_split(colnames(df.org.all), "\\|"), function(x) x[7])

species_id <- paste0(rep("sp_"), 1:ncol(df.org.all))

taxonomy.df <- data.frame(species_id = species_id, kingdom = tax.kingdom, phyla = tax.phyla, 
                          class = tax.class, order = tax.order, family = tax.family, genus= tax.genus, 
                          species = tax.species
)

colnames(df.org.all) <- species_id

# Filtering taxonomic data.frame 
# (include the species present in at least 20% of the “baskets” categories)
num.zeros <- colSums(df.org.all==0)
save.zeros  <-  (num.zeros<round(nrow(df.org.all)*0.8))
df.org.no_zeros  <-  df.org.all[,(save.zeros==TRUE)]

# Make Bray-Curtis MDS biplot
df.org.no_zeros.sbs <- df.org.no_zeros[ind_sbs,]

bray.mds <- metaMDS(df.org.no_zeros.sbs)
bray.mds.points <- as.data.frame(bray.mds$points)
bray.mds.points <- cbind(group.df[ind_sbs,], bray.mds.points)

bray.mds.points$group <- factor(bray.mds.points$group, levels = c("settle", "not_settle", "stay", "gone"))
bray.mds.points$dataset <- factor(bray.mds.points$dataset, levels = c("FMT", "Control"))

bray_biplot <- ggplot(bray.mds.points, aes(MDS1, MDS2, col = group, shape = dataset))+
    geom_point(size = 2.5)+
    xlim(c(-1.5,1.5))+
    ylim(c(-1.5,1.5))+
    scale_color_brewer(palette = "Set1")+
    theme_bw()+
    theme(legend.position = "right")+
    scale_shape_manual(values = c(19, 0))

# Save Bray-Curtis MDS biplot
svg(filename="FIGURES/biplot_mp2.svg", width=4.5, height=3.3, pointsize=12)
bray_biplot
dev.off()