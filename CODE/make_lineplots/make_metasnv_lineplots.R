workdir <- "/home/acari/github/RECAST_project/"

setwd(workdir)

library(stringr)
library(vegan)
library(ggplot2)
library(reshape2)
library(dplyr)

# init metadata
sample.metadata <- read.csv("DATA/sample.metadata", sep = "\t", stringsAsFactors = F)
sample.metadata.sbs <- sample.metadata[sample.metadata$Status != "Donor",]

# init basename and common value names
basename.id <- sample.metadata.sbs$Sample[sample.metadata.sbs$Time == 0]
common_lines <- c("Subject", "Dataset", "Status", "Time", "mOTU")

# make parcer
get_mann_dist <- function(path){
    
    list.f <- list.files(path, pattern = "mann.dist")
    
    df.all <- NULL
    for (k in list.f){
        
        df <- read.csv(paste0(path, k), sep = "\t", row.names = 1)
        colnames(df) <- gsub(".bam", "", colnames(df))
        rownames(df) <- gsub(".bam", "", rownames(df))
        colnames(df) <- gsub("\\.", "\\-", colnames(df))
        
        if (ncol(df[colnames(df) %in% basename.id]) > 0){
            
            df <- df[colnames(df) %in% basename.id]
            df <- cbind(rownames(df), df)
            colnames(df)[1] <- "samples"
            rownames(df) <- 1:nrow(df)
            
            df <- merge(df, sample.metadata.sbs, by = 1)
            
            df.2 <- NULL
            for (i in unique(df$Subject)){
                
                df.sbs <- df[df$Subject == i,]
                
                if (nrow(df.sbs[df.sbs$samples %in% colnames(df.sbs),]) > 0){
                    df.sbs <- df.sbs[colnames(df.sbs) %in% c(as.character(df.sbs$samples), common_lines)]
                    colnames(df.sbs)[1] <- "dist"
                    
                    df.2 <- rbind(df.sbs, df.2)    
                }
            }
            
            df.2$mOTU <- gsub(".filtered.mann.dist", "", k)
            
            df.all <- rbind(df.2, df.all)    
        }
    }
    
    df.all <- df.all[c(6,2:5,1)]
    return(df.all)
}

# parsing data by dataset
SPB18 <- get_mann_dist(path = "DATA/metasnv/SPB18/")
LEE17 <- get_mann_dist(path = "DATA/metasnv/LEE17/")
VRIEZE12_allogenic_don_11 <- get_mann_dist(path = "DATA/metasnv/VRIEZE12_allogenic_don_11/")
VRIEZE12_allogenic_don_19 <- get_mann_dist(path = "DATA/metasnv/VRIEZE12_allogenic_don_19/")
VRIEZE12_allogenic_don_8 <- get_mann_dist(path = "DATA/metasnv/VRIEZE12_allogenic_don_8/")
VRIEZE12_autologous <- get_mann_dist(path = "DATA/metasnv/VRIEZE12_autologous/")
VOIGT15 <- get_mann_dist(path = "DATA/metasnv/VOIGT15/")

# make summary data frame
VRIEZE12_allogenic <- rbind(VRIEZE12_allogenic_don_11, VRIEZE12_allogenic_don_19, VRIEZE12_allogenic_don_8)
VRIEZE12_allogenic$Dataset <- paste0(VRIEZE12_allogenic$Dataset, "_", VRIEZE12_allogenic$Status)
VRIEZE12_autologous$Dataset <- paste0(VRIEZE12_autologous$Dataset, "_", VRIEZE12_autologous$Status)
metasnv.df <- rbind(SPB18, LEE17, VRIEZE12_allogenic, VRIEZE12_autologous, VOIGT15)[-4]

# reshape data by mean distance
metasnv.df.r <- metasnv.df[-1] %>% 
    group_by(Subject, Dataset, Time) %>% 
    summarise(mean = mean(dist), sd = sd(dist))
metasnv.df.r <- as.data.frame(metasnv.df.r)

# make plot
metasnv.df.r$Dataset <- factor(metasnv.df.r$Dataset, 
                           levels = c("SPB18", "LEE17", "VRIEZE12_Allogenic", 
                                      "VRIEZE12_Autologous", "VOIGT15"))
metasnv.df.r.sbs <- metasnv.df.r[metasnv.df.r$Time < 400,]
CR <- max(metasnv.df.r.sbs$mean[metasnv.df.r.sbs$Time != 0 & metasnv.df.r.sbs$Dataset == "VOIGT15"])

metasnv_lineplots <- ggplot(metasnv.df.r.sbs, aes(Time, mean, group = Subject, col = Dataset))+
    geom_line(size = 1)+
    geom_point(size = 2)+
    facet_wrap(~Dataset, ncol = 5)+
    geom_hline(yintercept = CR, col = "red", linetype = 2)+
    ylab("Manhattan distance")+
    xlab("Time, days")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x  = element_text(angle=90, vjust=0.5, size=6.5))+
    scale_x_log10(breaks=c(sort(unique(metasnv.df.r.sbs$Time))))+
    theme(legend.position = "none")+
    scale_color_brewer(palette="Set1")+
    ylim(c(0,0.6))

svg(filename="FIGURES/metasnv_lineplots.svg", width=8, height=2, pointsize=12)
metasnv_lineplots
dev.off()

# wilcoxon runk sum test: are there any differences in the distance change over time between datasets?
metasnv.df.sbs <- metasnv.df[metasnv.df$Time >0,]

metasnv.df.sbs.wilcox <- pairwise.wilcox.test(
    metasnv.df.sbs$dist, 
    metasnv.df.sbs$Dataset, 
    p.adjust.method="BH"
)