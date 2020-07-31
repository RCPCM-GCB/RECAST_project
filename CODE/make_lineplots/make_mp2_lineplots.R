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
sample.metadata.sbs <- sample.metadata[c(1:4,7)]
sample.metadata.sbs <- sample.metadata.sbs[sample.metadata.sbs$Status != "Donor",]
basename.id <- sample.metadata.sbs$Sample[sample.metadata.sbs$Time == 0]
bray.dist <- merge(bray.dist[c("sample", basename.id)], sample.metadata.sbs, by = 1)

# reshape data without Spb dataset
subject_id <- unique(bray.dist$Subject)
common_lines <- c("Subject", "Dataset", "Status", "Time")

bray.res <- NULL
for (ids in subject_id){
  
  row_id <- which(!is.na(str_extract(bray.dist$sample, ids)))
  col_id <- which(!is.na(str_extract(colnames(bray.dist), paste(c(ids, common_lines), collapse = "|"))))
  
  bray.dist.sbs <- bray.dist[row_id, col_id]
  colnames(bray.dist.sbs)[1] <- "bray_dist"
  
  bray.res <- rbind(bray.res, bray.dist.sbs)  
}

# reshape data Spb dataset by subject
subject_id_v1 <- c("SPB_T21|SPB_T23|SPB_T25|SPB_T27|SPB_T29|SPB_T30")
subject_id_v2 <- c("SPB_T31|SPB_T33|SPB_T35|SPB_T36")
subject_id_v3 <- c("SPB_T37|SPB_T39|SPB_T41|SPB_T42")

V1_dist  <- bray.dist[which(!is.na(str_extract(bray.dist$sample, subject_id_v1))),
  which(!is.na(str_extract(colnames(bray.dist), paste(c(subject_id_v1, common_lines), collapse = "|"))))]
V2_dist  <- bray.dist[which(!is.na(str_extract(bray.dist$sample, subject_id_v2))),
  which(!is.na(str_extract(colnames(bray.dist), paste(c(subject_id_v2, common_lines), collapse = "|"))))]
V3_dist  <- bray.dist[which(!is.na(str_extract(bray.dist$sample, subject_id_v3))),
  which(!is.na(str_extract(colnames(bray.dist), paste(c(subject_id_v3, common_lines), collapse = "|"))))]

colnames(V1_dist)[1] <- "bray_dist"
colnames(V2_dist)[1] <- "bray_dist"
colnames(V3_dist)[1] <- "bray_dist"

# make full reshape Bray-Curtis dataset
bray.res <- rbind(bray.res, V1_dist, V2_dist, V3_dist)
bray.res <- bray.res[c(2:5, 1)]
bray.res$Dataset[bray.res$Dataset == "VRIEZE12"] <- paste0(bray.res$Dataset[bray.res$Dataset == "VRIEZE12"], "_",
                                                           bray.res$Status[bray.res$Dataset == "VRIEZE12"])
bray.res <- bray.res[-3]

write.table(bray.res, "OUTPUT/mp2_from_baseline.txt", sep = "\t", quote = F, row.names = F)

# make plot
bray.res$Dataset <- factor(bray.res$Dataset, 
                           levels = c("SPB18", "LEE17", "VRIEZE12_Allogenic", 
                                                        "VRIEZE12_Autologous", "VOIGT15"))
bray.res.sbs <- bray.res[bray.res$Time < 400,]
CR <- max(bray.res.sbs$bray_dist[bray.res.sbs$Time != 0 & bray.res.sbs$Dataset == "VOIGT15"])

mp2_lineplots <- ggplot(bray.res.sbs, aes(Time, bray_dist, group = Subject, col = Dataset))+
    geom_line(size = 1)+
    geom_point(size = 2)+
    geom_hline(yintercept = CR, col = "red", linetype = 2)+
    facet_wrap(~Dataset, ncol = 5)+
    ylab("Bray-Curtis dissimilarity")+
    xlab("Time, days")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x  = element_text(angle=90, vjust=0.5, size=6.5))+
    scale_x_log10(breaks=c(sort(unique(bray.res.sbs$Time))))+
    theme(legend.position = "none")+
    scale_color_brewer(palette="Set1")+
    ylim(c(0,0.8))
  
svg(filename="FIGURES/mp2_lineplots.svg", width=8, height=2, pointsize=12)
mp2_lineplots
dev.off()

# wilcoxon runk sum test: are there any differences in the distance change over time between datasets?
bray.res.sbs <- bray.res[bray.res$Time >0,]

bray.res.sbs.wilcox <- pairwise.wilcox.test(
    bray.res.sbs$bray_dist, 
    bray.res.sbs$Dataset, 
    p.adjust.method="BH"
)