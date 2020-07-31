workdir <- "/home/acari/github/RECAST_project/"

setwd(workdir)

library(stringr)
library(pheatmap)
library(vegan)
library(ggplot2)
library(reshape2)
library(randomForest)

# init sample metadata 
sample.metadata <- read.csv("DATA/sample.metadata", sep = "\t", stringsAsFactors = F)
sample.metadata.sbs <- sample.metadata[c(1:4,6:7)]
sample.metadata.sbs <- sample.metadata.sbs[sample.metadata.sbs$Status == "Allogenic",]

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

# init sorting/non_sortind taxonomic data
df.sorting.org <- read.csv("DATA/df_reads_count_sorting.org", sep = "\t", stringsAsFactors = F)[c(3,1,2)]

# init HMP2012 data frame
df.hmp2012 <- read.csv("DATA/df.hmp2012.org", sep = "\t", row.names = 1, stringsAsFactors = F)
colnames(df.hmp2012) <- gsub("\\.", "\\|", colnames(df.hmp2012))
df.hmp2012 <- data.frame(taxa = colnames(df.hmp2012), mean_relab_hmp2012 = apply(df.hmp2012, 2, mean))
rownames(df.hmp2012) <- 1:nrow(df.hmp2012)

# make group data frame
groups <- sapply(str_split(df.sorting.org$sample, "\\_", n = 3), function(x) x[3])
pattern <- paste0(paste0("_", unique(groups)), collapse = "|")
non_sorting_sample_id <- gsub(pattern, "", df.sorting.org$sample)

group.df <- data.frame(sample = df.sorting.org$sample, non_sorting_sample_id, group = groups)
group.df$group <- as.character(group.df$group)
group.df$group[group.df$group == "s"] <- "settle"
group.df$group[group.df$group == "n"] <- "not_settle"
group.df$group <- gsub("came_", "", group.df$group)

ind_sns <- rownames(group.df[group.df$group %in% c("settle", "not_settle"),])
ind_sg <- rownames(group.df[group.df$group %in% c("stay", "gone"),])

# init non_sortind taxonomic data
df.non_sorting.org <- read.csv("DATA/df.non_sorting.mp2.org", sep = "\t", row.names = 1, stringsAsFactors = F)
colnames(df.non_sorting.org) <- gsub("\\.", "\\|", colnames(df.non_sorting.org))

baseline.id <- sample.metadata.sbs$Sample[sample.metadata.sbs$Time == 0]
df.baseline <- df.non_sorting.org[baseline.id,]
df.baseline.melt <- melt(cbind(rownames(df.baseline), df.baseline))
colnames(df.baseline.melt) <- c("Sample", "taxa", "baseline_abundance")
df.baseline.melt$Sample <- as.character(df.baseline.melt$Sample)

df.baseline.melt.merge <- merge(base.meta, df.baseline.melt, by = 1)[-1]
df.baseline.melt.merge$Sample <- as.character(df.baseline.melt.merge$Sample)
df.baseline.melt.merge$taxa <- as.character(df.baseline.melt.merge$taxa)

donor.id <- unique(sample.metadata.sbs$Donor)
donor.id <- c("FAT_DON_11-22-pooled", "FAT_DON_19-22-0-0", "FAT_DON_8-22-0-0", donor.id[c(4:6)])
df.donor <- df.non_sorting.org[donor.id,]
df.donor.melt <- melt(cbind(rownames(df.donor), df.donor))
colnames(df.donor.melt) <- c("Donor", "taxa", "donor_abundance")
df.donor.melt$Donor <- (df.donor.melt$Donor)
df.donor.melt$taxa <- as.character(df.donor.melt$taxa)

# subset and melt data frame
df.sns <- df.sorting.org[ind_sns,]
df.sns <- cbind(group.df[ind_sns,], df.sns[-1])
df.sns_spread <- spread(df.sns[-1], group, n_reads, fill = 0)

df.sg <- df.sorting.org[ind_sg,]
df.sg <- cbind(group.df[ind_sg,], df.sg[-1])
df.sg_spread <- spread(df.sg[-1], group, n_reads, fill = 0)

# add donor and baseline bacteria relative abundance
df.melt.meta_sns <- merge(sample.metadata.sbs, df.sns_spread, by = 1)
df.melt.meta_sns.baseline <- merge(df.melt.meta_sns, df.baseline.melt.merge, by = c("Sample", "taxa"))
df.donor.melt$Donor <- as.character(df.donor.melt$Donor)
df.donor.melt$Donor <- gsub("-22-pooled|-22-0-0", "", df.donor.melt$Donor)
df.int_sns <- merge(df.melt.meta_sns.baseline, df.donor.melt, by = c("Donor", "taxa"))
df.int_sns <- merge(df.int_sns, df.hmp2012, by = "taxa")

df.melt.meta_sg <- merge(sample.metadata.sbs, df.sg_spread, by = 1)
df.melt.meta_sg.baseline <- merge(df.melt.meta_sg, df.baseline.melt.merge, by = c("Sample", "taxa"))
df.int_sg <- merge(df.melt.meta_sg.baseline, df.donor.melt, by = c("Donor", "taxa"))
df.int_sg <- merge(df.int_sg, df.hmp2012, by = "taxa")

# make taxa data frame
taxa_sns <- as.character(df.int_sns$taxa)
taxa_sns <- gsub("k__|p__|c__|o__|f__|g__|s__", "", taxa_sns)
df.taxa_sns <- NULL

for (i in 1:7){
    df.taxa_sns <- cbind(sapply(str_split(taxa_sns, "\\|"), function(x) x[i]), 
                         df.taxa_sns)
}

df.taxa_sns <- as.data.frame(df.taxa_sns)
colnames(df.taxa_sns) <- c("species", "genus", "family", "order", "class", "phyla", "kingdom")
df.comp_sns <- cbind(df.taxa_sns, df.int_sns)

###
taxa_sg <- as.character(df.int_sg$taxa)
taxa_sg <- gsub("k__|p__|c__|o__|f__|g__|s__", "", taxa_sg)

df.taxa_sg <- NULL
for (i in 1:7){
    df.taxa_sg <- cbind(sapply(str_split(taxa_sg, "\\|"), function(x) x[i]), 
                        df.taxa_sg)
}

df.taxa_sg <- as.data.frame(df.taxa_sg)
colnames(df.taxa_sg) <- colnames(df.taxa_sns)
df.comp_sg <- cbind(df.taxa_sg, df.int_sg)

###

df.comp_sns$settle <- df.comp_sns$settle+1
df.comp_sns$not_settle <- df.comp_sns$not_settle+1
df.comp_sns$Status <- df.comp_sns$settle/(df.comp_sns$settle+df.comp_sns$not_settle)
df.comp_sns <- df.comp_sns[-c(1,2,7,8,10,15,16)]
df.comp_sns <- df.comp_sns[c(6,5,7,9,1:4,8,10:12)]
df.comp_sns$Status[df.comp_sns$Status > 0.5] <- 1
df.comp_sns$Status[df.comp_sns$Status < 0.5] <- 0

df.comp_sg$stay <- df.comp_sg$stay+1
df.comp_sg$gone <- df.comp_sg$gone+1
df.comp_sg$Status <- df.comp_sg$stay/(df.comp_sg$stay+df.comp_sg$gone)
df.comp_sg <- df.comp_sg[-c(1,2,7,8,10,15,16)]
df.comp_sg <- df.comp_sg[c(6,5,7,9,1:4,8,10:12)]
df.comp_sg$Status[df.comp_sg$Status > 0.5] <- 1
df.comp_sg$Status[df.comp_sg$Status < 0.5] <- 0

# write obtained files
write.table(df.comp_sns, "OUTPUT/df_sns", quote = F, row.names = T, sep = "\t")
write.table(df.comp_sg, "OUTPUT/df_sg", quote = F, row.names = T, sep = "\t")