workdir <- "/home/acari/github/RECAST_project/"

setwd(workdir)

library(stringr)
library(zCompositions)
library(pheatmap)
library(vegan)
library(ggplot2)
library(tidyr)
library(reshape2)
library(dplyr)

# import sample metadata
sample.metadata <- read.csv("DATA/sample.metadata", sep = "\t")

# import KEGG data
ko_gene <- read.csv("DATA/kegg_links/ko_gene", sep = "\t", stringsAsFactors = F)

# import datasets
df.sort <- read.csv("DATA/df.lantibioitcs.sorting.txt", sep = "\t", stringsAsFactors = F, row.names = NULL)
# df.sort$width_cov <- df.sort$width_cov/df.sort$gene_length
# df.sort <- df.sort[df.sort$width_cov > 0.95,][-c(3,4)]

df.sort.ko <- merge(df.sort, ko_gene, by = "gene_name")
df.sort.ko <- df.sort.ko %>% 
    group_by(sample, ko) %>% 
    summarise(hits = sum(hits))
df.sort.ko <- as.data.frame(df.sort.ko)

df.sort_benchmark <- read.csv("DATA/df.lantibioitcs.sorting_benchmark.txt", sep = "\t", stringsAsFactors = F, row.names = NULL)
# df.sort_benchmark$width_cov <- df.sort_benchmark$width_cov/df.sort_benchmark$gene_length
# df.sort_benchmark <- df.sort_benchmark[df.sort_benchmark$width_cov > 0.95,][-c(3,4)]

df.sort_benchmark.ko <- merge(df.sort_benchmark, ko_gene, by = "gene_name")
df.sort_benchmark.ko <- df.sort_benchmark.ko %>% 
    group_by(sample, ko) %>% 
    summarise(hits = sum(hits))
df.sort_benchmark.ko <- as.data.frame(df.sort_benchmark.ko)

# make groups table
groups <- sapply(str_split(df.sort.ko$sample, "\\_", n = 3), function(x) x[3])
pattern <- paste0("_",unique(groups), collapse = "|")
non_sorting_sample_id <- gsub(pattern, "", df.sort.ko$sample)

group.df.fmt <- data.frame(
    Sample = df.sort.ko$sample, 
    non_sorting_sample_id, 
    group = groups, 
    dataset = "FMT"
)

groups.benchmark <- sapply(str_split(df.sort_benchmark.ko$sample, "\\_", n = 4), function(x) x[4])
pattern.benchmark <- paste0(paste0("_donor", "_", unique(groups.benchmark)), collapse = "|")
non_sorting_sample_id_benchmark <- gsub(pattern.benchmark, "", df.sort_benchmark.ko$sample)
non_sorting_sample_id_benchmark <- gsub("_daisy|_bugkiller|_scavenger|_peacemaker|_tigress", "", non_sorting_sample_id_benchmark)

group.df.benchmark <- data.frame(
    Sample = df.sort_benchmark.ko$sample, 
    non_sorting_sample_id = non_sorting_sample_id_benchmark, 
    group = groups.benchmark, 
    dataset = "Control"
)

group.df <- rbind(group.df.fmt, group.df.benchmark)

group.df$group <- as.character(group.df$group)
group.df$group[group.df$group == "s"] <- "settle"
group.df$group[group.df$group == "n"] <- "not_settle"
group.df$group <- gsub("came_", "", group.df$group)

ind_sbs <- as.numeric(rownames(group.df[!group.df$group %in% c("from_donor", "from_both", "from_before", "itself"),]))

# make combine ko table
df.sort.all <- rbind(df.sort.ko, df.sort_benchmark.ko)

write.table(df.sort.all, "OUTPUT/df_lantibiotics.txt", sep = "\t", quote = F, row.names = F)

# make count ko table
df.sort.sp <- spread(df.sort.all[ind_sbs,], ko, hits, fill = 0)
rownames(df.sort.sp) <- df.sort.sp$sample
df.sort.sp <- df.sort.sp[-1]

# Make Bray-Curtis MDS biplot
bray.mds <- metaMDS(df.sort.sp, distance = "bray", k = 2)
bray.mds.points <- as.data.frame(bray.mds$points)

group.df.sbs <- unique(group.df[group.df$Sample %in% rownames(bray.mds.points),])
group.df.sbs$Sample <- as.character(group.df.sbs$Sample)

bray.mds.points <- bray.mds.points[group.df.sbs$Sample,]
bray.mds.points <- cbind(group.df.sbs, bray.mds.points)
bray.mds.points <- bray.mds.points[-1]

bray.mds.points$group <- factor(bray.mds.points$group, levels = c("settle", "not_settle", "stay", "gone"))

bray_biplot <- ggplot(bray.mds.points, aes(MDS1, MDS2, col = group, shape = dataset))+
    geom_point(size = 2)+
    # xlim(c(-0.5,0.5))+
    # ylim(c(-0.5,0.5))+
    scale_color_brewer(palette = "Set1")+
    theme_bw()+
    theme(legend.position = "right")+
    scale_shape_manual(values = c(19,0))

svg(filename="FIGURES/biplot_lantibiotics.svg", width=4.5, height=3.3, pointsize=12)
bray_biplot
dev.off()

# add to combine table group data
df.sort.group <- merge(group.df, df.sort.all, by = 1)[-1]
df.sort.group.sbs <- df.sort.group[df.sort.group$group %in% c("settle", "not_settle", "stay", "gone"),]

df.sort.sum <- df.sort.group.sbs %>% 
    group_by(group, dataset) %>% 
    summarise(mean = mean(hits), sd = sd(hits))
df.sort.sum <- as.data.frame(df.sort.sum)

df.sort.sum$group <- factor(df.sort.sum$group, levels = c("settle", "not_settle", "stay", "gone"))

lanti_barplot <- ggplot(df.sort.sum, aes(group, mean, fill = group))+
    geom_bar(color="black", position = position_dodge(), stat = "identity", width = 0.55)+
    geom_errorbar(aes(ymin = mean, ymax = mean+sd), width=.2,
                  position=position_dodge(.9))+
    facet_wrap(~dataset)+
    theme_classic()+
    scale_fill_brewer(palette = "Set1")+
    theme(legend.position = "bottom")+
    ylim(c(0,12000))+
    ylab("Hits")+
    xlab("Groups")

svg(filename="FIGURES/barplot_lantibiotics.svg", width=4.5, height=3.3, pointsize=12)
lanti_barplot
dev.off()

############################
df.line_2 <- merge(group.df, df.sort.all, by = 1)[-1]
df.line_2 <- merge(df.line_2, sample.metadata, by = 1)
df.line_2 <- df.line_2[-c(3,8:10)]
df.line_2 <- df.line_2[df.line_2$group %in% c("from_donor", "from_both", "from_before", "itself"),]

df.line_2_2 <- df.line_2 %>% 
    group_by(group, Dataset, Time) %>%
    summarise(hits = sum(hits))
df.line_2_2 <- as.data.frame(df.line_2_2)

df.line_2_3 <- df.line_2 %>% 
    group_by(Dataset, Time) %>%
    summarise(sum_hits = sum(hits))
df.line_2_3 <- as.data.frame(df.line_2_3)
colnames(df.line_2_2)
df.line_2_4 <- merge(df.line_2_2, df.line_2_3, by = c("Dataset", "Time"))
df.line_2_4$relab <- 100*(df.line_2_4$hits/df.line_2_4$sum_hits)


df.line_2_4$Dataset <- factor(df.line_2_4$Dataset, 
                                  levels = c("SPB18", "LEE17", "VRIEZE12", "VOIGT15"))
df.line_2_4.sbs <- df.line_2_4[df.line_2_4$Time < 400,]    

df.line_2_4.sbs$group <- factor(df.line_2_4.sbs$group, levels = c("from_donor", "from_both", "from_before", "itself"))
df.line_2_4.sbs <- df.line_2_4.sbs[df.line_2_4.sbs$Time != 45 &  df.line_2_4.sbs$Time != 75,]

lineplot_lanti_line_2 <- ggplot(df.line_2_4.sbs, aes(Time, relab, group = group, fill = group))+
    geom_area()+
    facet_wrap(~Dataset, ncol = 5)+
    ylab("mean % of reads")+
    xlab("log Time, days")+
    theme_bw()+
    scale_x_log10(breaks=c(sort(unique(df.line_2_4.sbs$Time))))+
    theme(legend.position = "none")+
    scale_fill_manual(values = c("#E41A1C", "#FF7F00", "#4DAF4A", "#984EA3"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.x  = element_text(angle=90, vjust=0.5, size=6.5))

svg(filename="FIGURES/lineplot_lantibiotics.svg", width=6, height=1.5, pointsize=12)
lineplot_lanti_line_2
dev.off()

################
df.line_3 <- merge(group.df, df.sort.all, by = 1)[-1]
df.line_3 <- merge(df.line_3, sample.metadata, by = 1)
df.line_3 <- df.line_3[-c(3,8:10)]
df.line_3 <- df.line_3[df.line_3$group %in% c("settle", "not_settle"),]

df.line_3_2 <- df.line_3 %>% 
    group_by(group, Dataset, Time) %>%
    summarise(hits = sum(hits))
df.line_3_2 <- as.data.frame(df.line_3_2)

df.line_3_3 <- df.line_3 %>% 
    group_by(Dataset, Time) %>%
    summarise(sum_hits = sum(hits))
df.line_3_3 <- as.data.frame(df.line_3_3)

df.line_3_4 <- merge(df.line_3_2, df.line_3_3, by = c("Dataset", "Time"))
df.line_3_4$relab <- 100*(df.line_3_4$hits/df.line_3_4$sum_hits)

df.line_3_4$Dataset <- factor(df.line_3_4$Dataset, 
                              levels = c("SPB18", "LEE17", "VRIEZE12", "VOIGT15"))
df.line_3_4.sbs <- df.line_3_4[df.line_3_4$Time < 400,]    

df.line_3_4.sbs$group <- factor(df.line_3_4.sbs$group, levels = c("settle", "not_settle"))
df.line_3_4.sbs <- df.line_3_4.sbs[df.line_3_4.sbs$Time != 45 &  df.line_3_4.sbs$Time != 75,]

lineplot_lanti_line_3 <- ggplot(df.line_3_4.sbs, aes(Time, relab, group = group, fill = group))+
    geom_area()+
    facet_wrap(~Dataset, ncol = 5)+
    ylab("mean % of reads")+
    xlab("log Time, days")+
    theme_bw()+
    scale_x_log10(breaks=c(sort(unique(df.line_3_4.sbs$Time))))+
    theme(legend.position = "none")+
    scale_fill_brewer(palette="Set1")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.x  = element_text(angle=90, vjust=0.5, size=6.5))

svg(filename="FIGURES/lineplot_lantibiotics_sns.svg", width=6, height=1.5, pointsize=12)
lineplot_lanti_line_3
dev.off()

########
df.line_4 <- merge(group.df, df.sort.all, by = 1)[-1]
df.line_4 <- merge(df.line_4, sample.metadata, by = 1)
df.line_4 <- df.line_4[-c(3,8:10)]
df.line_4 <- df.line_4[df.line_4$group %in% c("stay", "gone"),]

df.line_4_2 <- df.line_4 %>% 
    group_by(group, Dataset, Time) %>%
    summarise(hits = sum(hits))
df.line_4_2 <- as.data.frame(df.line_4_2)

df.line_4_3 <- df.line_4 %>% 
    group_by(Dataset, Time) %>%
    summarise(sum_hits = sum(hits))
df.line_4_3 <- as.data.frame(df.line_4_3)

df.line_4_4 <- merge(df.line_4_2, df.line_4_3, by = c("Dataset", "Time"))
df.line_4_4$relab <- 100*(df.line_4_4$hits/df.line_4_4$sum_hits)

df.line_4_4$Dataset <- factor(df.line_4_4$Dataset, 
                              levels = c("SPB18", "LEE17", "VRIEZE12", "VOIGT15"))
df.line_4_4.sbs <- df.line_4_4[df.line_4_4$Time < 400,]    

df.line_4_4.sbs$group <- factor(df.line_4_4.sbs$group, levels = c("stay", "gone"))
df.line_4_4.sbs <- df.line_4_4.sbs[df.line_4_4.sbs$Time != 45 &  df.line_4_4.sbs$Time != 75,]

lineplot_lanti_line_4 <- ggplot(df.line_4_4.sbs, aes(Time, relab, group = group, fill = group))+
    geom_area()+
    facet_wrap(~Dataset, ncol = 5)+
    ylab("mean % of reads")+
    xlab("log Time, days")+
    theme_bw()+
    scale_x_log10(breaks=c(sort(unique(df.line_4_4.sbs$Time))))+
    theme(legend.position = "none")+
    scale_fill_manual(values = c("#4DAF4A", "#984EA3"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.x  = element_text(angle=90, vjust=0.5, size=6.5))

svg(filename="FIGURES/lineplot_lantibiotics_sg.svg", width=6, height=1.5, pointsize=12)
lineplot_lanti_line_4
dev.off()