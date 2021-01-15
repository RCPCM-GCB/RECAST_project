###################################################################################
# This R code designed to made Figure 5C (part 1) "Distribution of antibiotic     #
# and nisin resistance genes in “basket” categories over time."                   #
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
library(zCompositions)
library(pheatmap)
library(vegan)
library(ggplot2)
library(tidyr)
library(reshape2)
library(dplyr)

source("CODE/mini_compositions_lib.R")

# Import tables
## Metadata
sample.metadata <- read.csv("DATA/sample.metadata", sep = "\t")

## MEGARes data
### Case
df.sort <- read.csv("DATA/df.megares.mech.sorting.txt", sep = "\t", stringsAsFactors = F, row.names = NULL)
df.sort <- spread(df.sort, Mechanism, Hits, fill = 0)

### Control
df.sort_benchmark <- read.csv("DATA/df.megares.mech.sorting_benchmark.txt", sep = "\t", stringsAsFactors = F, row.names = NULL)
df.sort_benchmark <- spread(df.sort_benchmark, Mechanism, Hits, fill = 0)

# Merge Merge Merge 
df.sort.all <- merge(df.sort, df.sort_benchmark, all = TRUE)
df.sort.all[is.na(df.sort.all)] <- 0

# Make groups table
groups <- sapply(str_split(df.sort$Sample, "\\_", n = 3), function(x) x[3])
pattern <- paste0("_",unique(groups), collapse = "|")
non_sorting_sample_id <- gsub(pattern, "", df.sort$Sample)

group.df.fmt <- data.frame(
    Sample = df.sort$Sample, 
    non_sorting_sample_id, 
    group = groups, 
    dataset = "FMT"
)

groups.benchmark <- sapply(str_split(df.sort_benchmark$Sample, "\\_", n = 4), function(x) x[4])
pattern.benchmark <- paste0(paste0("_donor", "_", unique(groups.benchmark)), collapse = "|")
non_sorting_sample_id_benchmark <- gsub(pattern.benchmark, "", df.sort_benchmark$Sample)
non_sorting_sample_id_benchmark <- gsub("_daisy|_bugkiller|_scavenger|_peacemaker|_tigress", "", non_sorting_sample_id_benchmark)

group.df.benchmark <- data.frame(
    Sample = df.sort_benchmark$Sample, 
    non_sorting_sample_id = non_sorting_sample_id_benchmark, 
    group = groups.benchmark, 
    dataset = "Control"
)

group.df <- rbind(group.df.fmt, group.df.benchmark)

group.df$group <- as.character(group.df$group)
group.df$group[group.df$group == "s"] <- "settle"
group.df$group[group.df$group == "n"] <- "not_settle"
group.df$group <- gsub("came_", "", group.df$group)

write.table(group.df, "OUTPUT/group.df.txt", sep = "\t", quote = F)

ind_sbs <- as.numeric(rownames(group.df[!group.df$group %in% c("from_donor", "from_both", "from_before", "itself"),]))

rownames(df.sort.all) <- df.sort.all$Sample
df.sort.all <- df.sort.all[-1]
df.sort.all <- df.sort.all[as.character(group.df$Sample),]

# Filtering taxonomic data.frame
num.zeros <- colSums(df.sort.all == 0)
save.zeros  <-  (num.zeros<round(nrow(df.sort.all)*0.8))
df.sort.no_zeros  <-  df.sort.all[,(save.zeros==TRUE)]

# write.table(df.sort.no_zeros, "OUTPUT/df_megares.txt", sep = "\t", quote = F)

# Make Bray-Curtis MDS biplot
df.sort.no_zeros.sbs <- df.sort.no_zeros[ind_sbs,]
df.sort.no_zeros.sbs <- df.sort.no_zeros.sbs[rowSums(df.sort.no_zeros.sbs) > 0,]

bray.mds <- metaMDS(df.sort.no_zeros.sbs, distance = "bray", k = 2)

bray.mds.points <- as.data.frame(bray.mds$points)
bray.mds.points <- cbind(group.df[group.df$Sample %in% rownames(bray.mds.points),], bray.mds.points)
bray.mds.points <- bray.mds.points[-1]

bray.mds.points$group <- factor(bray.mds.points$group, levels = c("settle", "not_settle", "stay", "gone"))

# bray_biplot <- ggplot(bray.mds.points, aes(MDS1, MDS2, col = group, shape = dataset))+
#     geom_point(size = 2)+
#     # xlim(c(-1.5,1.5))+
#     # ylim(c(-1.5,1.5))+
#     scale_color_brewer(palette = "Set1")+
#     theme_bw()+
#     theme(legend.position = "right")+
#     scale_shape_manual(values = c(19,0))
# 
# svg(filename="FIGURES/biplot_megares.svg", width=4.5, height=3.3, pointsize=12)
# bray_biplot
# dev.off()

df.sort.no_zeros.sbs_2 <- cbind(group.df[group.df$Sample %in% rownames(df.sort.no_zeros.sbs),][-1], df.sort.no_zeros.sbs)

df.sort.melt <- melt(df.sort.no_zeros.sbs_2[-9])
colnames(df.sort.melt)[c(4,5)] <- c("ARGs_group", "Hits")

df.sort.melt$group <- factor(df.sort.melt$group, levels = c("settle", "not_settle", "stay", "gone"))
df.sort.melt$non_sorting_sample_id <- as.character(df.sort.melt$non_sorting_sample_id)

df.sort.melt$Dataset <- NA
df.sort.melt$Dataset[which(!is.na(str_extract(df.sort.melt$non_sorting_sample_id, "FAT")))] <- "VRIEZE12"
df.sort.melt$Dataset[which(!is.na(str_extract(df.sort.melt$non_sorting_sample_id, "R0")))] <- "LEE17"
df.sort.melt$Dataset[which(!is.na(str_extract(df.sort.melt$non_sorting_sample_id, "SPB")))] <- "SPB18"
df.sort.melt$Dataset[is.na(df.sort.melt$Dataset)] <- "VOIGT15"

df.sort.melt_2 <- df.sort.melt %>% 
                    group_by(group, dataset) %>% 
                    summarise(mean = mean(Hits), sd = sd(Hits))
df.sort.melt_2 <- as.data.frame(df.sort.melt_2)

# megares_barplot <- ggplot(df.sort.melt_2, aes(group, mean, fill = group))+
#     geom_bar(color="black", position = position_dodge(), stat = "identity", width = 0.55)+
#     geom_errorbar(aes(ymin = mean, ymax = mean+sd), width=.2,
#                       position=position_dodge(.9))+
#     facet_wrap(~dataset, ncol = 2)+
#     theme_classic()+
#     scale_fill_brewer(palette = "Set1")+
#     theme(legend.position = "bottom")+
#     xlab("Groups")+
#     ylab("Hits")+
#     ylim(c(0,20000))
# 
# svg(filename="FIGURES/barplot_megares.svg", width=4.5, height=3.3, pointsize=12)
# megares_barplot
# dev.off()

df.sort.melt_3 <- df.sort.melt %>% 
    group_by(group, dataset, ARGs_group) %>% 
    summarise(mean = mean(Hits), sd = sd(Hits))
df.sort.melt_3 <- as.data.frame(df.sort.melt_3)

# megares_barplot_2 <- ggplot(df.sort.melt_3, aes(ARGs_group, mean, fill = group))+
#     geom_bar(color="black", position = position_dodge(), stat = "identity", width = 0.75)+
#     facet_wrap(~dataset, ncol = 1)+
#     theme_classic()+
#     scale_fill_brewer(palette = "Set1")+
#     theme(legend.position = "bottom")+
#     xlab("Groups of ARGs")+
#     ylab("Mean of hits")
# 
# svg(filename="FIGURES/barplot_megares_by_groups.svg", width=6.5, height=4.5, pointsize=12)
# megares_barplot_2
# dev.off()

# Wilcoxon testing summary hits
settle_not_settle_FMT <- wilcox.test(df.sort.melt$Hits[df.sort.melt$group == "settle" & df.sort.melt$dataset == "FMT"], 
    df.sort.melt$Hits[df.sort.melt$group == "not_settle" & df.sort.melt$dataset == "FMT"], 
    alternative = "greater")$p.val

stay_gone_FMT <- wilcox.test(df.sort.melt$Hits[df.sort.melt$group == "stay" & df.sort.melt$dataset == "FMT"], 
            df.sort.melt$Hits[df.sort.melt$group == "gone" & df.sort.melt$dataset == "FMT"], 
            alternative = "greater")$p.val

settle_not_settle_control <- wilcox.test(df.sort.melt$Hits[df.sort.melt$group == "settle" & df.sort.melt$dataset == "Control"], 
            df.sort.melt$Hits[df.sort.melt$group == "not_settle" & df.sort.melt$dataset == "Control"], 
            alternative = "greater")$p.val

stay_gone_control <- wilcox.test(df.sort.melt$Hits[df.sort.melt$group == "stay" & df.sort.melt$dataset == "Control"], 
            df.sort.melt$Hits[df.sort.melt$group == "gone" & df.sort.melt$dataset == "Control"], 
            alternative = "greater")$p.val

wilcox.megares.results <- data.frame(group = c("settle/not_settle", "stay/gone", "settle/not_settle", "stay/gone"), 
                                     dataset = c("FMT", "FMT", "Control", "Control"), 
           p.val = c(settle_not_settle_FMT, stay_gone_FMT, settle_not_settle_control, stay_gone_control))

wilcox.megares.results$p.adj <- p.adjust(wilcox.megares.results$p.val, method = "fdr")

# Wilcox testing hits over groups

arg_group <- as.character(unique(df.sort.melt$ARGs_group))

sns.wilcox.fmt <- NULL
for (i in arg_group){
    p.val <- wilcox.test(df.sort.melt$Hits[df.sort.melt$group == "settle" & df.sort.melt$dataset == "FMT" & df.sort.melt$ARGs_group == i], 
                df.sort.melt$Hits[df.sort.melt$group == "not_settle" & df.sort.melt$dataset == "FMT" & df.sort.melt$ARGs_group == i], 
                alternative = "greater")$p.val
    sns.wilcox.fmt <- rbind(sns.wilcox.fmt, data.frame(arg_group = i, dataset = "FMT", group = "settle/not_settle", p.val))
}

sns.wilcox.control <- NULL
for (i in arg_group){
    p.val <- wilcox.test(df.sort.melt$Hits[df.sort.melt$group == "settle" & df.sort.melt$dataset == "Control" & df.sort.melt$ARGs_group == i], 
                         df.sort.melt$Hits[df.sort.melt$group == "not_settle" & df.sort.melt$dataset == "Control" & df.sort.melt$ARGs_group == i], 
                         alternative = "greater")$p.val
    sns.wilcox.control <- rbind(sns.wilcox.control, data.frame(arg_group = i, dataset = "Control", group = "settle/not_settle", p.val))
}

sg.wilcox.FMT <- NULL
for (i in arg_group){
    p.val <- wilcox.test(df.sort.melt$Hits[df.sort.melt$group == "stay" & df.sort.melt$dataset == "FMT" & df.sort.melt$ARGs_group == i], 
                         df.sort.melt$Hits[df.sort.melt$group == "gone" & df.sort.melt$dataset == "FMT" & df.sort.melt$ARGs_group == i], 
                         alternative = "greater")$p.val
    sg.wilcox.FMT <- rbind(sg.wilcox.FMT, data.frame(arg_group = i, dataset = "FMT", group = "stay/gone", p.val))
}

sg.wilcox.control <- NULL
for (i in arg_group){
    p.val <- wilcox.test(df.sort.melt$Hits[df.sort.melt$group == "stay" & df.sort.melt$dataset == "Control" & df.sort.melt$ARGs_group == i], 
                         df.sort.melt$Hits[df.sort.melt$group == "gone" & df.sort.melt$dataset == "Control" & df.sort.melt$ARGs_group == i], 
                         alternative = "greater")$p.val
    sg.wilcox.control <- rbind(sg.wilcox.control, data.frame(arg_group = i, dataset = "Control", group = "stay/gone", p.val))
}

wilcox.args.all <- rbind(sns.wilcox.fmt, sg.wilcox.FMT, sns.wilcox.control, sg.wilcox.control)
wilcox.args.all$p.adj <- p.adjust(wilcox.args.all$p.val, method = "fdr")

wilcox.args.all.sbs <- wilcox.args.all[wilcox.args.all$p.adj < 0.01,]

#################
ind_sbs_2 <- as.numeric(rownames(group.df[group.df$group %in% c("from_donor", "from_both", "from_before", "itself"),]))
df_line2 <- df.sort.no_zeros[ind_sbs_2,]

df_line2 <- cbind(group.df[group.df$Sample %in% rownames(df_line2),][-1], df_line2)
df_line2.melt <- melt(df_line2)
df_line2.melt$non_sorting_sample_id <- as.character(df_line2.melt$non_sorting_sample_id)
df_line2.melt$non_sorting_sample_id <- gsub("_donor_came_from_donor", "", df_line2.melt$non_sorting_sample_id)

df_line2.melt <- merge(df_line2.melt, sample.metadata, by = 1)[-c(8:10)]
head(df_line2.melt)
colnames(df_line2.melt)[c(4,5)] <- c("arg", "hits")

df_line2.melt_2 <- df_line2.melt %>% 
    group_by(group, arg, Dataset, Time) %>%
    summarise(hits = sum(hits))
df_line2.melt_2 <- as.data.frame(df_line2.melt_2)

df_line2.melt_3 <- df_line2.melt %>% 
    group_by(arg, Dataset, Time) %>%
    summarise(sum_hits = sum(hits))
df_line2.melt_3 <- as.data.frame(df_line2.melt_3)

df_line2.melt_4 <- merge(df_line2.melt_2, df_line2.melt_3, by = c("arg", "Dataset", "Time"))
df_line2.melt_4$relab <- (df_line2.melt_4$hits/df_line2.melt_4$sum_hits)*100
    
df_line2.melt_4$Dataset <- factor(df_line2.melt_4$Dataset, 
                          levels = c("SPB18", "LEE17", "VRIEZE12", "VOIGT15"))
df_line2.melt_4.sbs <- df_line2.melt_4[df_line2.melt_4$Time < 400,]    

df_dp_2.sbs_2 <- df_line2.melt_4.sbs
df_dp_2.sbs_2$group <- factor(df_dp_2.sbs_2$group, levels = c("from_donor", "from_both", "from_before", "itself"))
df_dp_2.sbs_2 <- df_dp_2.sbs_2[df_dp_2.sbs_2$Time != 45 &  df_dp_2.sbs_2$Time != 75,]

tetracyclines <- ggplot(df_dp_2.sbs_2[df_dp_2.sbs_2$arg == "Tetracyclines",], aes(Time, relab, group = group, fill = group))+
    geom_area()+
    facet_wrap(~Dataset, ncol = 5)+
    ylab("mean % of reads")+
    xlab("log Time, days")+
    theme_bw()+
    scale_x_log10(breaks=c(sort(unique(df_dp_2.sbs_2$Time))))+
    theme(legend.position = "none")+
    # ggtitle("Tetracyclines")+
    scale_fill_manual(values = c("#E41A1C", "#FF7F00", "#4DAF4A", "#984EA3"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.x  = element_text(angle=90, vjust=0.5, size=6.5))

svg(filename="FIGURES/lineplots/lineplot_tetracyclines.svg", width=6, height=1.5, pointsize=12)
tetracyclines
dev.off()

mls <- ggplot(df_dp_2.sbs_2[df_dp_2.sbs_2$arg == "MLS",], aes(Time, relab, group = group, fill = group))+
    geom_area()+
    facet_wrap(~Dataset, ncol = 5)+
    ylab("mean % of reads")+
    xlab("log Time, days")+
    theme_bw()+
    scale_x_log10(breaks=c(sort(unique(df_dp_2.sbs_2$Time))))+
    theme(legend.position = "none")+
    # ggtitle("MLS")+
    scale_fill_manual(values = c("#E41A1C", "#FF7F00", "#4DAF4A", "#984EA3"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.x  = element_text(angle=90, vjust=0.5, size=6.5))

svg(filename="FIGURES/lineplots/lineplot_mls.svg", width=6, height=1.5, pointsize=12)
mls
dev.off()

aminoglycosides <- ggplot(df_dp_2.sbs_2[df_dp_2.sbs_2$arg == "Aminoglycosides",], aes(Time, relab, group = group, fill = group))+
    geom_area()+
    facet_wrap(~Dataset, ncol = 5)+
    ylab("mean % of reads")+
    xlab("log Time, days")+
    theme_bw()+
    scale_x_log10(breaks=c(sort(unique(df_dp_2.sbs_2$Time))))+
    theme(legend.position = "none")+
    # ggtitle("Aminoglycosides")+
    scale_fill_manual(values = c("#E41A1C", "#FF7F00", "#4DAF4A", "#984EA3"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.x  = element_text(angle=90, vjust=0.5, size=6.5))

svg(filename="FIGURES/lineplots/lineplot_aminoglycosides.svg", width=6, height=1.5, pointsize=12)
aminoglycosides
dev.off()

#################
ind_sbs_3 <- as.numeric(rownames(group.df[group.df$group %in% c("settle", "not_settle"),]))

df_line3 <- df.sort.no_zeros[ind_sbs_3,]
df_line3 <- cbind(group.df[group.df$Sample %in% rownames(df_line3),][-1], df_line3)
df_line3.melt <- melt(df_line3)

df_line3.melt$non_sorting_sample_id <- as.character(df_line3.melt$non_sorting_sample_id)
df_line3.melt$non_sorting_sample_id <- gsub("_donor_came_from_donor", "", df_line3.melt$non_sorting_sample_id)

df_line3.melt <- merge(df_line3.melt, sample.metadata, by = 1)[-c(3,8:10)]
colnames(df_line3.melt)[c(3,4)] <- c("arg", "hits")

df_line3.melt_2 <- df_line3.melt %>% 
    group_by(group, arg, Dataset, Time) %>%
    summarise(hits = sum(hits))
df_line3.melt_2 <- as.data.frame(df_line3.melt_2)

df_line3.melt_3 <- df_line3.melt %>% 
    group_by(arg, Dataset, Time) %>%
    summarise(sum_hits = sum(hits))
df_line3.melt_3 <- as.data.frame(df_line3.melt_3)

df_line3.melt_4 <- merge(df_line3.melt_2, df_line3.melt_3, by = c("arg", "Dataset", "Time"))
df_line3.melt_4$relab <- (df_line3.melt_4$hits/df_line3.melt_4$sum_hits)*100

df_line3.melt_4$Dataset <- factor(df_line3.melt_4$Dataset, 
                                  levels = c("SPB18", "LEE17", "VRIEZE12", "VOIGT15"))
df_line3.melt_4.sbs <- df_line3.melt_4[df_line3.melt_4$Time < 400,]    

df_line3.melt_4.sbs$group <- factor(df_line3.melt_4.sbs$group, levels = c("settle", "not_settle"))
df_line3.melt_4.sbs <- df_line3.melt_4.sbs[df_line3.melt_4.sbs$Time != 45 &  df_line3.melt_4.sbs$Time != 75,]

tetracyclines_sns <- ggplot(df_line3.melt_4.sbs[df_line3.melt_4.sbs$arg == "Tetracyclines",], aes(Time, relab, group = group, fill = group))+
    geom_area()+
    facet_wrap(~Dataset, ncol = 5)+
    ylab("mean % of reads")+
    xlab("log Time, days")+
    theme_bw()+
    scale_x_log10(breaks=c(sort(unique(df_dp_2.sbs_2$Time))))+
    theme(legend.position = "none")+
    scale_fill_brewer(palette="Set1")+
    # ggtitle("Tetracyclines")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.x  = element_text(angle=90, vjust=0.5, size=6.5))

svg(filename="FIGURES/lineplots/lineplot_tetracyclines_sns.svg", width=6, height=1.5, pointsize=12)
tetracyclines_sns
dev.off()

mls_sns <- ggplot(df_line3.melt_4.sbs[df_line3.melt_4.sbs$arg == "MLS",], aes(Time, relab, group = group, fill = group))+
    geom_area()+
    facet_wrap(~Dataset, ncol = 5)+
    ylab("mean % of reads")+
    xlab("log Time, days")+
    theme_bw()+
    scale_x_log10(breaks=c(sort(unique(df_dp_2.sbs_2$Time))))+
    theme(legend.position = "none")+
    scale_fill_brewer(palette="Set1")+
    # ggtitle("MLS")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.x  = element_text(angle=90, vjust=0.5, size=6.5))

svg(filename="FIGURES/lineplots/lineplot_mls_sns.svg", width=6, height=1.5, pointsize=12)
mls_sns
dev.off()

aminoglycosides_sns <- ggplot(df_line3.melt_4.sbs[df_line3.melt_4.sbs$arg == "Aminoglycosides",], aes(Time, relab, group = group, fill = group))+
    geom_area()+
    facet_wrap(~Dataset, ncol = 5)+
    ylab("mean % of reads")+
    xlab("log Time, days")+
    theme_bw()+
    scale_x_log10(breaks=c(sort(unique(df_dp_2.sbs_2$Time))))+
    theme(legend.position = "none")+
    scale_fill_brewer(palette="Set1")+
    # ggtitle("Aminoglycosides")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.x  = element_text(angle=90, vjust=0.5, size=6.5))

svg(filename="FIGURES/lineplots/lineplot_aminoglycosides_sns.svg", width=6, height=1.5, pointsize=12)
aminoglycosides_sns
dev.off()

#### 
ind_sbs_4 <- as.numeric(rownames(group.df[group.df$group %in% c("stay", "gone"),]))

df_line4 <- df.sort.no_zeros[ind_sbs_4,]
df_line4 <- cbind(group.df[group.df$Sample %in% rownames(df_line4),][-1], df_line4)
df_line4.melt <- melt(df_line4)

df_line4.melt <- merge(df_line4.melt, sample.metadata, by = 1)[-c(3,8:10)]
colnames(df_line4.melt)[c(3,4)] <- c("arg", "hits")

df_line4.melt_2 <- df_line4.melt %>% 
    group_by(group, arg, Dataset, Time) %>%
    summarise(hits = sum(hits))
df_line4.melt_2 <- as.data.frame(df_line4.melt_2)

df_line4.melt_3 <- df_line4.melt %>% 
    group_by(arg, Dataset, Time) %>%
    summarise(sum_hits = sum(hits))
df_line4.melt_3 <- as.data.frame(df_line4.melt_3)

df_line4.melt_4 <- merge(df_line4.melt_2, df_line4.melt_3, by = c("arg", "Dataset", "Time"))
df_line4.melt_4$relab <- (df_line4.melt_4$hits/df_line4.melt_4$sum_hits)*100

df_line4.melt_4$Dataset <- factor(df_line4.melt_4$Dataset, 
                                  levels = c("SPB18", "LEE17", "VRIEZE12", "VOIGT15"))
df_line4.melt_4.sbs <- df_line4.melt_4[df_line4.melt_4$Time < 400,]    

df_line4.melt_4.sbs$group <- factor(df_line4.melt_4.sbs$group, levels = c("stay", "gone"))
df_line4.melt_4.sbs <- df_line4.melt_4.sbs[df_line4.melt_4.sbs$Time != 45 &  df_line4.melt_4.sbs$Time != 75,]

tetracyclines_sg <- ggplot(df_line4.melt_4.sbs[df_line4.melt_4.sbs$arg == "Tetracyclines",], aes(Time, relab, group = group, fill = group))+
    geom_area()+
    facet_wrap(~Dataset, ncol = 5)+
    ylab("mean % of reads")+
    xlab("log Time, days")+
    theme_bw()+
    scale_x_log10(breaks=c(sort(unique(df_dp_2.sbs_2$Time))))+
    theme(legend.position = "none")+
    scale_fill_manual(values = c("#4DAF4A", "#984EA3"))+
    # ggtitle("Tetracyclines")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.x  = element_text(angle=90, vjust=0.5, size=6.5))

svg(filename="FIGURES/lineplots/lineplot_tetracyclines_sg.svg", width=6, height=1.5, pointsize=12)
tetracyclines_sg
dev.off()

mls_sg <- ggplot(df_line4.melt_4.sbs[df_line4.melt_4.sbs$arg == "MLS",], aes(Time, relab, group = group, fill = group))+
    geom_area()+
    facet_wrap(~Dataset, ncol = 5)+
    ylab("mean % of reads")+
    xlab("log Time, days")+
    theme_bw()+
    scale_x_log10(breaks=c(sort(unique(df_dp_2.sbs_2$Time))))+
    theme(legend.position = "none")+
    scale_fill_manual(values = c("#4DAF4A", "#984EA3"))+
    # ggtitle("MLS")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.x  = element_text(angle=90, vjust=0.5, size=6.5))

svg(filename="FIGURES/lineplots/lineplot_mls_sg.svg", width=6, height=1.5, pointsize=12)
mls_sg
dev.off()

aminoglycosides_sg <- ggplot(df_line4.melt_4.sbs[df_line4.melt_4.sbs$arg == "Aminoglycosides",], aes(Time, relab, group = group, fill = group))+
    geom_area()+
    facet_wrap(~Dataset, ncol = 5)+
    ylab("mean % of reads")+
    xlab("log Time, days")+
    theme_bw()+
    scale_x_log10(breaks=c(sort(unique(df_dp_2.sbs_2$Time))))+
    theme(legend.position = "none")+
    scale_fill_manual(values = c("#4DAF4A", "#984EA3"))+
    # ggtitle("Aminoglycosides")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.x  = element_text(angle=90, vjust=0.5, size=6.5))

svg(filename="FIGURES/lineplots/lineplot_aminoglycosides_sg.svg", width=6, height=1.5, pointsize=12)
aminoglycosides_sg
dev.off()