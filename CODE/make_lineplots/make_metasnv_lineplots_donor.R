workdir <- "/home/acari/github/RECAST_project/"

setwd(workdir)

library(stringr)
library(vegan)
library(ggplot2)
library(reshape2)
library(dplyr)

# init metadata
sample.metadata <- read.csv("DATA/sample.metadata", sep = "\t", stringsAsFactors = F)
sample.metadata.sbs <- sample.metadata[sample.metadata$Status == "Allogenic",]

sample.metadata.sbs$Donor[sample.metadata.sbs$Subject == "FAT_006"] <- "FAT_DON_11-22-0-5"
sample.metadata.sbs$Donor[sample.metadata.sbs$Subject == "FAT_008"] <- "FAT_DON_11-22-0-5"
sample.metadata.sbs$Donor[sample.metadata.sbs$Subject == "FAT_015"] <- "FAT_DON_11-22-0-5"
sample.metadata.sbs$Donor[sample.metadata.sbs$Subject == "FAT_012"] <- "FAT_DON_8-22-0-0"
sample.metadata.sbs$Donor[sample.metadata.sbs$Subject == "FAT_020"] <- "FAT_DON_19-22-0-0"

# init basename and common value names
donor.id <- unique(sample.metadata.sbs$Donor)[-5]
common_lines <- c("Subject", "Dataset", "Donor", "Time", "mOTU")

# parsing data by dataset
get_mann_dist <- function(path){
    list.f <- list.files(path, pattern = "mann.dist")
    
    df.all <- NULL
    for (k in list.f){
        
        df <- read.csv(paste0(path, k), sep = "\t", row.names = 1)
        colnames(df) <- gsub(".bam", "", colnames(df))
        rownames(df) <- gsub(".bam", "", rownames(df))
        colnames(df) <- gsub("\\.", "\\-", colnames(df))
        
        if (ncol(df[colnames(df) %in% donor.id]) > 0){
            
            df <- df[colnames(df) %in% donor.id]
            df <- cbind(rownames(df), df)
            colnames(df)[1] <- "samples"
            rownames(df) <- 1:nrow(df)
            
            df <- merge(df, sample.metadata.sbs, by = 1)
            
            if (nrow(df) > 0){
                
                colnames(df)[2] <- "dist"
                df$mOTU <- gsub(".filtered.mann.dist", "", k)
                
                df.all <- rbind(df, df.all)
            }
        }
        
    }
    return(df.all)
}

SPB18 <- get_mann_dist(path = "DATA/metasnv/SPB18/")
LEE17 <- get_mann_dist(path = "DATA/metasnv/LEE17/")
VRIEZE12_allogenic_don_11 <- get_mann_dist(path = "DATA/metasnv/VRIEZE12_allogenic_don_11/")
VRIEZE12_allogenic_don_19 <- get_mann_dist(path = "DATA/metasnv/VRIEZE12_allogenic_don_19/")
VRIEZE12_allogenic_don_8 <- get_mann_dist(path = "DATA/metasnv/VRIEZE12_allogenic_don_8/")

metasnv.df <- rbind(SPB18, LEE17, VRIEZE12_allogenic_don_11, VRIEZE12_allogenic_don_19, VRIEZE12_allogenic_don_8)
metasnv.df <- metasnv.df[c("dist", common_lines)]

# reshape data by mean distance
metasnv.df.r <- metasnv.df %>% 
    group_by(Subject, Dataset, Time) %>% 
    summarise(mean = mean(dist), sd = sd(dist))
metasnv.df.r <- as.data.frame(metasnv.df.r)

metasnv.df.r$Dataset <- factor(metasnv.df.r$Dataset, 
                               levels = c("SPB18", "LEE17", "VRIEZE12"))
CR <- 0.1140371
    
metasnv_lineplots <- ggplot(metasnv.df.r, aes(Time, mean, group = Subject, col = Dataset))+
    geom_line(size = 1)+
    geom_point(size = 2)+
    geom_hline(yintercept = CR, col = "red", linetype = 2)+
    facet_wrap(~Dataset, ncol = 5)+
    ylab("Manhattan distance")+
    xlab("Time, days")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x  = element_text(angle=90, vjust=0.5, size=6.5))+
    scale_x_log10(breaks=c(sort(unique(metasnv.df.r$Time))))+
    theme(legend.position = "none")+
    scale_color_brewer(palette="Set1")+
    ylim(c(0,0.6))

svg(filename="FIGURES/metasnv_lineplots_donor.svg", width=5, height=2, pointsize=12)
metasnv_lineplots
dev.off()

# wilcoxon runk sum test: are there any differences in the distance change over time between datasets?
metasnv.df.sbs <- metasnv.df[metasnv.df$Time >0,]

metasnv.df.sbs.wilcox <- pairwise.wilcox.test(
    metasnv.df.sbs$dist, 
    metasnv.df.sbs$Dataset, 
    p.adjust.method="BH"
)