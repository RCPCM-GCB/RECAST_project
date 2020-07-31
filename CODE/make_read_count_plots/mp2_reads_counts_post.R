workdir <- "/home/acari/github/RECAST_project/"

library(ggplot2)
library(tidyr)
library(stringr)
library(reshape2)
library(dplyr)
library(gridExtra)

setwd(workdir)

# import data tables
df_sorting <- read.csv("DATA/df_reads_count_sorting.org", sep = "\t", stringsAsFactors = F)
df_sorting_bench <- read.csv("DATA/df_reads_count_benchmark_sorting.org", sep = "\t", stringsAsFactors = F)
df_sorting.all <- rbind(df_sorting, df_sorting_bench)

df_non_sorting <- read.csv("DATA/df_reads_count_non_sorting.org", sep = "\t", stringsAsFactors = F)
df_non_sorting_bench <- read.csv("DATA/df_reads_count_benchmark_non_sorting.org", sep = "\t", stringsAsFactors = F)

sample.metadata <- read.csv("DATA/sample.metadata", sep = "\t", stringsAsFactors = F)
sample.metadata.sbs <- sample.metadata[c(1:3,7)]

# add groups to sorting data frame 
groups <- sapply(str_split(df_sorting$sample, "\\_", n = 3), function(x) x[3])
pattern <- paste0("_",unique(groups), collapse = "|")
non_sorting_sample_id <- gsub(pattern, "", df_sorting$sample)

group.df.fmt <- data.frame(
    sample = df_sorting$sample, 
    non_sorting_sample_id, 
    group = groups, 
    bechmark = "FMT"
)

groups.benchmark <- sapply(str_split(df_sorting_bench$sample, "\\_", n = 4), function(x) x[4])
pattern.benchmark <- paste0(paste0("_donor", "_", unique(groups.benchmark)), collapse = "|")
non_sorting_sample_id_benchmark <- gsub(pattern.benchmark, "", df_sorting_bench$sample)
non_sorting_sample_id_benchmark <- gsub("_daisy|_bugkiller|_scavenger|_peacemaker|_tigress", "", non_sorting_sample_id_benchmark)

group.df.benchmark <- data.frame(
    sample = df_sorting_bench$sample, 
    non_sorting_sample_id = non_sorting_sample_id_benchmark, 
    group = groups.benchmark, 
    bechmark = "Control"
)

group.df <- rbind(group.df.fmt, group.df.benchmark)

group.df$group <- as.character(group.df$group)
group.df$group[group.df$group == "s"] <- "settle"
group.df$group[group.df$group == "n"] <- "not_settle"
group.df$group <- gsub("came_", "", group.df$group)

####
df_sort_gr <- cbind(group.df[-1], df_sorting.all[c(1,2)])
df_sort_gr2 <- merge(sample.metadata.sbs, df_sort_gr, by = 1)
df_sort_gr3 <- df_sort_gr2[df_sort_gr2$group %in% c("from_donor", "from_both", "from_before", "itself"),]

###
df_dp <- df_sort_gr3[-1] %>% 
    group_by(Subject, Dataset, Time, group, bechmark) %>% 
    summarise(n_reads = sum(n_reads))
df_dp <- as.data.frame(df_dp)

df_dp_all <- df_dp[-4] %>% 
    group_by(Subject, Dataset, Time, bechmark) %>% 
    summarise(all_reads = sum(n_reads))
df_dp_all <- as.data.frame(df_dp_all)

df_dp_2 <- merge(df_dp, df_dp_all, by = c("Subject", "Dataset", "Time", "bechmark"))
df_dp_2$relab <- (df_dp_2$n_reads/df_dp_2$all_reads)*100

# make plot
df_dp_2$Dataset <- factor(df_dp_2$Dataset, 
                               levels = c("SPB18", "LEE17", "VRIEZE12", "VOIGT15"))
df_dp_2.sbs <- df_dp_2[df_dp_2$Time < 400,]

df_dp_2.sbs_2 <- df_dp_2.sbs[-1] %>% 
    group_by(Dataset, Time, bechmark, group) %>% 
    summarise(mean_relab = mean(relab))
df_dp_2.sbs_2 <- as.data.frame(df_dp_2.sbs_2)

df_dp_2.sbs_2$group <- factor(df_dp_2.sbs_2$group, levels = c("from_donor", "from_both", "from_before", "itself"))
df_dp_2.sbs_2 <- df_dp_2.sbs_2[df_dp_2.sbs_2$Time != 45 &  df_dp_2.sbs_2$Time != 75,]

lineplot <- ggplot(df_dp_2.sbs_2, aes(Time, mean_relab, group = group, fill = group))+
    geom_area()+
    # geom_point(size = 2)+
    facet_wrap(~Dataset, ncol = 5)+
    ylab("mean % of reads")+
    xlab("log Time, days")+
    theme_bw()+
    scale_x_log10(breaks=c(sort(unique(df_dp_2.sbs_2$Time))))+
    theme(legend.position = "right")+
    scale_fill_manual(values = c("#E41A1C", "#FF7F00", "#4DAF4A", "#984EA3"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x  = element_text(angle=90, vjust=0.5, size=6.5))

svg(filename="FIGURES/lineplot_4_groups.svg", width=10, height=2.5, pointsize=12)
lineplot
dev.off()