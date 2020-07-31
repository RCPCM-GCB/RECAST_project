# set work directory
workdir <- "/home/acari/github/RECAST_project/"
setwd(workdir)

# import libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# import data
df_non_sorting <- read.csv("DATA/df_sample_non_sorting_reads_count.txt", sep = "\t", stringsAsFactors = F)

df_sorting <- read.csv("DATA/df_sample_sorting_reads_count.txt", sep = "\t", stringsAsFactors = F)
df_sorting_benchmark <- read.csv("DATA/df_sample_sorting_benchmark_reads_count.txt", sep = "\t", stringsAsFactors = F)

sample.metadata <- read.csv("DATA/sample.metadata", sep = "\t", stringsAsFactors = F)
sample.metadata[sample.metadata$Status == "Allogenic",]
sample.metadata.sbs <- sample.metadata[c(1:4,6:7)]
sample.metadata.sbs <- sample.metadata.sbs[!sample.metadata.sbs$Status %in% c("Autologous", "Donor"),]

# baseline samples idx
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

# make sorting groups table
groups <- sapply(str_split(df_sorting$Sample, "\\_", n = 3), function(x) x[3])
pattern <- paste0("_",unique(groups), collapse = "|")
non_sorting_sample_id <- gsub(pattern, "", df_sorting$Sample)

sample.metadata.sbs_2 <- sample.metadata.sbs[sample.metadata.sbs$Status == "Allogenic",]
sample.metadata.sbs_2 <- sample.metadata.sbs_2[c(1,5,6)]
sample.metadata.sbs_2 <- sample.metadata.sbs_2[sample.metadata.sbs_2$Time != 0,][-3]

sample.metadata.sbs_2$Donor[which(!is.na(str_extract(sample.metadata.sbs_2$Sample, "FAT_006")))] <- "FAT_DON_11-22-pooled"
sample.metadata.sbs_2$Donor[which(!is.na(str_extract(sample.metadata.sbs_2$Sample, "FAT_008")))] <- "FAT_DON_11-22-pooled"
sample.metadata.sbs_2$Donor[which(!is.na(str_extract(sample.metadata.sbs_2$Sample, "FAT_015")))] <- "FAT_DON_11-22-pooled"
sample.metadata.sbs_2$Donor[which(!is.na(str_extract(sample.metadata.sbs_2$Sample, "FAT_012")))] <- "FAT_DON_8-22-0-0"
sample.metadata.sbs_2$Donor[which(!is.na(str_extract(sample.metadata.sbs_2$Sample, "FAT_020")))] <- "FAT_DON_19-22-0-0"
colnames(sample.metadata.sbs_2) <- c("non_sorting_sample_id", "donor")

group.df.fmt <- data.frame(
    Sample = df_sorting$Sample, 
    non_sorting_sample_id, 
    group = groups, 
    dataset = "FMT"
)

group.df.fmt <- merge(group.df.fmt, sample.metadata.sbs_2, by = "non_sorting_sample_id")[c(2,1,3,5,4)]
colnames(group.df.fmt)
groups.benchmark <- sapply(str_split(df_sorting_benchmark$Sample, "\\_", n = 4), function(x) x[4])
pattern.benchmark <- paste0(paste0("_donor", "_", unique(groups.benchmark)), collapse = "|")

non_sorting_sample_id_benchmark <- gsub(pattern.benchmark, "", df_sorting_benchmark$Sample)
non_sorting_sample_id_benchmark <- gsub("_daisy|_bugkiller|_scavenger|_peacemaker|_tigress", "", non_sorting_sample_id_benchmark)
non_sorting_sample_id_benchmark <- gsub("_donor_came_from_donor", "", non_sorting_sample_id_benchmark)

donor.benchmark <- gsub("donor_", "", df_sorting_benchmark$Sample)
donor.benchmark <- sapply(str_split(donor.benchmark, "\\_", n = 4), function(x) x[2])
donor.benchmark <- paste0(donor.benchmark, "-11-0-0")

group.df.benchmark <- data.frame(
    Sample = df_sorting_benchmark$Sample, 
    non_sorting_sample_id = non_sorting_sample_id_benchmark, 
    group = groups.benchmark,
    donor = donor.benchmark,
    dataset = "Control"
)

group.df <- rbind(group.df.fmt, group.df.benchmark)

group.df$group <- as.character(group.df$group)
group.df$group[group.df$group == "s"] <- "settle"
group.df$group[group.df$group == "n"] <- "not_settle"
group.df$group <- gsub("came_", "", group.df$group)

# merge with base-meta
colnames(base.meta)[2] <- "non_sorting_sample_id"
group.df <- merge(group.df, base.meta, by = "non_sorting_sample_id")
group.df <- merge(group.df[-5], sample.metadata.sbs[c(1,3)], by = 1)
group.df <- group.df[c(2,1,3,5,4,6)]
colnames(group.df)[6] <- "dataset"

df_sorting_all <- rbind(df_sorting, df_sorting_benchmark)
df_sorting_all <- merge(group.df, df_sorting_all, by = 1)

# settle-not_settle (donor sample)
df_sns <- df_sorting_all[df_sorting_all$group %in% c("settle", "not_settle"),]
df_sns <- df_sns[c(2,3,5:7)]

donor_samples <- df_non_sorting
colnames(donor_samples)[c(1:2)] <- c("donor", "N_donor_reads")

df_sns_d <- merge(df_sns, donor_samples, by = "donor")
df_sns_d$relab <- 100*(df_sns_d$N_reads/df_sns_d$N_donor_reads)
df_sns_d <- merge(df_sns_d[c(2,1,3,4,7)], sample.metadata.sbs[c(1,2,6)], by = 1)
df_sns_d <- df_sns_d[-1]
df_sns_d <- df_sns_d[df_sns_d$dataset != "VOIGT15",]

df_sns_d$group <- factor(df_sns_d$group, levels = c("settle", "not_settle"))
df_sns_d$Subject <- factor(df_sns_d$Subject, levels = c("V1", "V2", "V3", "R01", "R02", "FAT_006", "FAT_008", "FAT_012", "FAT_015", "FAT_020"))

SPB18 <- df_sns_d[df_sns_d$dataset == "SPB18",]
LEE17 <- df_sns_d[df_sns_d$dataset == "LEE17",]
VRIEZE12 <- df_sns_d[df_sns_d$dataset == "VRIEZE12",]

spb_reads_count_plot <- ggplot(SPB18, aes(Time, relab, group = group, fill = group))+
    geom_area()+
    facet_wrap(~Subject, ncol = 1)+
    ylab("% of reads")+
    xlab("Time, days")+
    theme_bw()+
    xlim(0,300)+
    scale_x_log10(breaks=c(5,30,45,75,300))+
    theme(legend.position = "none")+
    scale_y_continuous(breaks = c(100,0))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_fill_brewer(palette="Set1")+
    ylim(c(0,100))

svg(filename="FIGURES/spb_reads_count_plot.svg", width=2.5, height=3.5, pointsize=12)
spb_reads_count_plot
dev.off()

lee_reads_count_plot <- ggplot(LEE17, aes(Time, relab, group = group, fill = group))+
    geom_area()+
    facet_wrap(~Subject, ncol = 1)+
    ylab("% of reads")+
    xlab("Time, days")+
    theme_bw()+
    xlim(5,300)+
    scale_x_log10(breaks=c(28,56))+
    theme(legend.position = "none")+
    scale_y_continuous(breaks = c(100,0))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_fill_brewer(palette="Set1")+
    ylim(c(0,100))

svg(filename="FIGURES/lee_reads_count_plot.svg", width=2.5, height=2.5, pointsize=12)
lee_reads_count_plot
dev.off()

vrieze_reads_count_plot <- ggplot(VRIEZE12, aes(Time, relab, group = group, fill = group))+
    geom_area()+
    facet_wrap(~Subject, ncol = 1)+
    ylab("% of reads")+
    xlab("Time, days")+
    theme_bw()+
    xlim(5,300)+
    scale_x_log10(breaks=c(0,2,14,42,84))+
    theme(legend.position = "none")+
    scale_y_continuous(breaks = c(100,0))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_fill_brewer(palette="Set1")+
    ylim(c(0,100))

svg(filename="FIGURES/vrieze_reads_count_plot.svg", width=2.5, height=5.6, pointsize=12)
vrieze_reads_count_plot
dev.off()

# stay / gone (baseline sample)
df_sg <- df_sorting_all[df_sorting_all$group %in% c("stay", "gone"),]
df_sg <- df_sg[c(2:4,6:7)]

baseline_samples <- df_non_sorting
colnames(baseline_samples)[c(1:2)] <- c("baseline", "N_baseline_reads")

df_sg_d <- merge(df_sg, baseline_samples, by = "baseline")
df_sg_d$relab <- 100*(df_sg_d$N_reads/df_sg_d$N_baseline_reads)
df_sg_d <- merge(df_sg_d[c(2,1,3,4,7)], sample.metadata.sbs[c(1,2,6)], by = 1)
df_sg_d <- df_sg_d[-1]
df_sg_d <- df_sg_d[df_sg_d$dataset != "VOIGT15",]

df_sg_d$group <- factor(df_sg_d$group, levels = c("stay", "gone"))
df_sg_d$Subject <- factor(df_sg_d$Subject, levels = c("V1", "V2", "V3", "R01", "R02", "FAT_006", "FAT_008", "FAT_012", "FAT_015", "FAT_020"))

SPB18_sg <- df_sg_d[df_sg_d$dataset == "SPB18",]
LEE17_sg <- df_sg_d[df_sg_d$dataset == "LEE17",]
VRIEZE12_sg <- df_sg_d[df_sg_d$dataset == "VRIEZE12",]

spb_reads_count_plot_sg <- ggplot(SPB18_sg, aes(Time, relab, group = group, fill = group))+
    geom_area()+
    facet_wrap(~Subject, ncol = 1)+
    ylab("% of reads")+
    xlab("Time, days")+
    theme_bw()+
    xlim(0,300)+
    scale_x_log10(breaks=c(5,30,45,75,300))+
    theme(legend.position = "none")+
    scale_y_continuous(breaks = c(100,0))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_fill_manual(values = c("#4DAF4A", "#984EA3"))+
    ylim(c(0,100))

svg(filename="FIGURES/spb_reads_count_plot_sg.svg", width=2.5, height=3.5, pointsize=12)
spb_reads_count_plot_sg
dev.off()

lee_reads_count_plot_sg <- ggplot(LEE17_sg, aes(Time, relab, group = group, fill = group))+
    geom_area()+
    facet_wrap(~Subject, ncol = 1)+
    ylab("% of reads")+
    xlab("Time, days")+
    theme_bw()+
    xlim(5,300)+
    scale_x_log10(breaks=c(28,56))+
    theme(legend.position = "none")+
    scale_y_continuous(breaks = c(100,0))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_fill_manual(values = c("#4DAF4A", "#984EA3"))+
    ylim(c(0,100))

svg(filename="FIGURES/lee_reads_count_plot_sg.svg", width=2.5, height=2.5, pointsize=12)
lee_reads_count_plot_sg
dev.off()

vrieze_reads_count_plot_sg <- ggplot(VRIEZE12_sg, aes(Time, relab, group = group, fill = group))+
    geom_area()+
    facet_wrap(~Subject, ncol = 1)+
    ylab("% of reads")+
    xlab("Time, days")+
    theme_bw()+
    xlim(5,300)+
    scale_x_log10(breaks=c(0,2,14,42,84))+
    theme(legend.position = "none")+
    scale_y_continuous(breaks = c(100,0))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_fill_manual(values = c("#4DAF4A", "#984EA3"))+
    ylim(c(0,100))

svg(filename="FIGURES/vrieze_reads_count_plot_sg.svg", width=2.5, height=5.6, pointsize=12)
vrieze_reads_count_plot_sg
dev.off()

# from_donor / from_before / from_both / itself (after FMT sample)
df_4_gr <- df_sorting_all[df_sorting_all$group %in% c("from_donor", "from_before", "from_both", "itself"),]
df_4_gr <- df_4_gr[c(2:3,6,7)]

after_samples <- df_non_sorting
colnames(after_samples)[c(1:2)] <- c("non_sorting_sample_id", "N_after_reads")

df_4_gr_d <- merge(df_4_gr, after_samples, by = "non_sorting_sample_id")
df_4_gr_d$relab <- 100*(df_4_gr_d$N_reads/df_4_gr_d$N_after_reads)
df_4_gr_d <- df_4_gr_d[-c(4,5)]
df_4_gr_d <- merge(df_4_gr_d, sample.metadata.sbs[c(1,2,6)], by = 1)
df_4_gr_d <- df_4_gr_d[df_4_gr_d$dataset != "VOIGT15",]

df_4_gr_d$group <- factor(df_4_gr_d$group, levels = c("from_donor", "from_both", "from_before", "itself"))
df_4_gr_d$Subject <- factor(df_4_gr_d$Subject, levels = c("V1", "V2", "V3", "R01", "R02", "FAT_006", "FAT_008", "FAT_012", "FAT_015", "FAT_020"))

SPB18_4_gr <- df_4_gr_d[df_4_gr_d$dataset == "SPB18",]
LEE17_4_gr <- df_4_gr_d[df_4_gr_d$dataset == "LEE17",]
VRIEZE12_4_gr <- df_4_gr_d[df_4_gr_d$dataset == "VRIEZE12",]

spb_reads_count_plot_4_gr <- ggplot(SPB18_4_gr, aes(Time, relab, group = group, fill = group))+
    geom_area()+
    facet_wrap(~Subject, ncol = 1)+
    ylab("% of reads")+
    xlab("Time, days")+
    theme_bw()+
    xlim(0,300)+
    scale_x_log10(breaks=c(5,30,45,75,300))+
    theme(legend.position = "none")+
    scale_y_continuous(breaks = c(100,0))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_fill_manual(values = c("#E41A1C", "#FF7F00", "#4DAF4A", "#984EA3"))+
    ylim(c(0,100))

svg(filename="FIGURES/spb_reads_count_plot_4_gr.svg", width=2.5, height=3.5, pointsize=12)
spb_reads_count_plot_4_gr
dev.off()

lee_reads_count_plot_4_gr <- ggplot(LEE17_4_gr, aes(Time, relab, group = group, fill = group))+
    geom_area()+
    facet_wrap(~Subject, ncol = 1)+
    ylab("% of reads")+
    xlab("Time, days")+
    theme_bw()+
    xlim(5,300)+
    scale_x_log10(breaks=c(28,56))+
    theme(legend.position = "none")+
    scale_y_continuous(breaks = c(100,0))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_fill_manual(values = c("#E41A1C", "#FF7F00", "#4DAF4A", "#984EA3"))+
    ylim(c(0,100))

svg(filename="FIGURES/lee_reads_count_plot_4_gr.svg", width=2.5, height=2.5, pointsize=12)
lee_reads_count_plot_4_gr
dev.off()

vrieze_reads_count_plot_4_gr <- ggplot(VRIEZE12_4_gr, aes(Time, relab, group = group, fill = group))+
    geom_area()+
    facet_wrap(~Subject, ncol = 1)+
    ylab("% of reads")+
    xlab("Time, days")+
    theme_bw()+
    xlim(5,300)+
    scale_x_log10(breaks=c(0,2,14,42,84))+
    theme(legend.position = "none")+
    scale_y_continuous(breaks = c(100,0))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_fill_manual(values = c("#E41A1C", "#FF7F00", "#4DAF4A", "#984EA3"))+
    ylim(c(0,100))

svg(filename="FIGURES/vrieze_reads_count_plot_4_gr.svg", width=2.5, height=5.6, pointsize=12)
vrieze_reads_count_plot_4_gr
dev.off()