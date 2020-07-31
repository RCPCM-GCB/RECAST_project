workdir <- "/home/acari/github/RECAST_project/"

setwd(workdir)

library(stringr)
library(vegan)
library(ggplot2)
library(reshape2)

# init taxonomic data 
df.org <- read.csv("DATA/df.non_sorting.mp2.org", sep = "\t", row.names = 1, stringsAsFactors = F)
colnames(df.org) <- sapply(str_split(colnames(df.org), "\\."), function(x) tail(x,1))
colnames(df.org) <- gsub("s__", "", colnames(df.org))

# init metadata
sample.metadata <- read.csv("DATA/sample.metadata", sep = "\t", stringsAsFactors = F)

# filtering taxonomic data.frame
num.zeros <- colSums(df.org==0)
save.zeros  <-  (num.zeros<(nrow(df.org)*0.8))
df.org.no_zeros  <-  df.org[,(save.zeros==TRUE),]

# calculate Bray-Curtis dissimilarity 
bray.dist <- metaMDSdist(df.org.no_zeros)
bray.dist <- as.data.frame(as.matrix(bray.dist))
bray.dist <- cbind(rownames(bray.dist), bray.dist)
colnames(bray.dist)[1] <- "sample"
bray.dist$sample <- as.character(bray.dist$sample)

# merge Bray-Curtis dissimilarity data with sample metadata
sample.metadata.sbs <- sample.metadata[c(1:4,6,7)]
sample.metadata.sbs <- sample.metadata.sbs[sample.metadata.sbs$Status == "Allogenic",]

sample.metadata.sbs$Donor[sample.metadata.sbs$Subject == "FAT_006"] <- "FAT_DON_11-22-0-4"
sample.metadata.sbs$Donor[sample.metadata.sbs$Subject == "FAT_008"] <- "FAT_DON_11-22-0-5"
sample.metadata.sbs$Donor[sample.metadata.sbs$Subject == "FAT_015"] <- "FAT_DON_11-22-0-6"
sample.metadata.sbs$Donor[sample.metadata.sbs$Subject == "FAT_012"] <- "FAT_DON_8-22-0-0"
sample.metadata.sbs$Donor[sample.metadata.sbs$Subject == "FAT_020"] <- "FAT_DON_19-22-0-0"

donor.id <- unique(sample.metadata.sbs$Donor)
bray.dist <- merge(bray.dist[c("sample", donor.id)], sample.metadata.sbs, by = 1)

# reshape data without Spb dataset
subject_id <- unique(bray.dist$Subject)
common_lines <- c("Subject", "Dataset", "Status", "Time")

bray.res <- NULL
for (ids in subject_id){
    
    row_id <- which(bray.dist$Subject == ids)
    d.id <- unique(sample.metadata.sbs$Donor[sample.metadata.sbs$Subject == ids])
    
    bray.dist.sbs <- bray.dist[row_id, c(d.id, common_lines)]
    colnames(bray.dist.sbs)[1] <- "bray_dist"
    
    bray.res <- rbind(bray.res, bray.dist.sbs)  
}

write.table(bray.res, "OUTPUT/mp2_from_donor.txt", sep = "\t", quote = F, row.names = F)

bray.res$Dataset <- factor(bray.res$Dataset, 
                           levels = c("SPB18", "LEE17", "VRIEZE12"))
CR <- 0.4201082

# make plot
mp2_lineplots <- ggplot(bray.res, aes(Time, bray_dist, group = Subject, col = Dataset))+
    geom_line(size = 1)+
    geom_point(size = 2)+
    geom_hline(yintercept = CR, col = "red", linetype = 2)+
    facet_wrap(~Dataset, ncol = 5)+
    ylab("Bray-Curtis dissimilarity")+
    xlab("Time, days")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x  = element_text(angle=90, vjust=0.5, size=6.5))+
    scale_x_log10(breaks=c(sort(unique(bray.res$Time))))+
    theme(legend.position = "none")+
    scale_color_brewer(palette="Set1")+
    ylim(c(0,0.8))

svg(filename="FIGURES/mp2_lineplots_donor.svg", width=5, height=2, pointsize=12)
mp2_lineplots
dev.off()

# wilcoxon runk sum test: are there any differences in the distance change over time between datasets?
bray.res.sbs <- bray.res[bray.res$Time >0,]

bray.res.sbs.wilcox <- pairwise.wilcox.test(
    bray.res.sbs$bray_dist, 
    bray.res.sbs$Dataset, 
    p.adjust.method="BH"
)