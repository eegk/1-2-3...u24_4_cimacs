### Load packages for analysis
#############################################################################################################################################
libs<-c("corrplot","ROTS","Hmisc","tcR","tidyr","limma","NMF","tsne","Rtsne" ,"UpSetR",
        "AnnotationDbi","RColorBrewer","papmap","sva","GO.db","fitdistrplus","ff","plyranges",
        "annotables","Rsamtools","GenomicFeatures","ggarrange","pheatmap","rms",
        "dplyr","DESeq2","ggplot2","ggrepel","edgeR","doParallel","ggradar","colorspace",
        "reshape2","rmarkdown","org.Hs.eg.db","treemapify",
        "variancePartition","ggpubr","factoextra","Mfuzz", "universalmotif","ggbio","GenomicRanges",
        "scales","tibble","RColorBrewer","tidyr","ggplot2","reshape2","circlize",
        "colorspace","Vennerable","enrichR","cowplot","data.table",
        "ROCR","mlbench","caretEnsemble","gbm","rpart","UpSetR","SCENIC","AUCell","RcisTarget","plyr",
        "tidyverse","hrbrthemes","fmsb","colormap","viridis","survminer","survival","ggalluvial")
lapply(libs, require, character.only = TRUE) ; rm(libs)
#if (!requireNamespace("BiocManager", quietly = TRUE)) ; install.packages("BiocManager") ; BiocManager::install("survminer")
#############################################################################################################################################
setwd("/Users/gonzae34/Documents/projects_gnjatic/U24/Icahn School of Medicine at Mount Sinai (MSSM) CIMAC/")
#############################################################################################################################################
### Olink data analysis for U24 Sacha Gnjatic project
### Initial analysis: Mar 16 2021
### Last review: Mar 16 2021
### Edgar E. G. Kozlova
#############################################################################################################################################

### E4412
#############################################################################################################################################
#############################################################################################################################################
### E4412 has been analyzed extensively elsewhere
#############################################################################################################################################


### S1400I
#############################################################################################################################################
#############################################################################################################################################

### Data processing
#############################################################################################################################################
cimac_s1400i_npx <- readxl::read_excel("S1400I/Olink/CIMAC_S1400I_AllPlates_NPX.xlsx")

### inspect
### dim(cimac_s1400i_npx) #234 101
### cimac_s1400i_npx[1:10,1:5]
### cimac_s1400i_npx[1:10,1:5] #7
### cimac_s1400i_npx[231:234,1:5] #232

### Prptein names
cimac_s1400i_protein_info <- as.data.frame( cbind(t(cimac_s1400i_npx[2:5,]),as.data.frame(t(cimac_s1400i_npx[233:234,])) ))
colnames(cimac_s1400i_protein_info) <- cimac_s1400i_protein_info[1,]
cimac_s1400i_protein_info <- cimac_s1400i_protein_info[-1,]

### values
cimac_s1400i_npx_raw <- as.data.frame(cimac_s1400i_npx[7:231,])
cimac_s1400i_npx_raw$`Olink NPX Manager 2.0.1.175`

### duplicated names
test <- cimac_s1400i_npx_raw$CIMAC_S1400I_AllPlates[grep("IPC", cimac_s1400i_npx_raw$CIMAC_S1400I_AllPlates)] 
test <- paste(1:length(test),test,sep="___")
cimac_s1400i_npx_raw$CIMAC_S1400I_AllPlates[grep("IPC", cimac_s1400i_npx_raw$CIMAC_S1400I_AllPlates)]  <- test

### duplicated names
test <- cimac_s1400i_npx_raw$CIMAC_S1400I_AllPlates[grep("NC", cimac_s1400i_npx_raw$CIMAC_S1400I_AllPlates)] 
test <- paste(1:length(test),test,sep="___")
cimac_s1400i_npx_raw$CIMAC_S1400I_AllPlates[grep("NC", cimac_s1400i_npx_raw$CIMAC_S1400I_AllPlates)]  <- test

### labels
rownames(cimac_s1400i_npx_raw) <- cimac_s1400i_npx_raw$CIMAC_S1400I_AllPlates

### clean
cimac_s1400i_npx_raw <- cimac_s1400i_npx_raw[,-1] ; rm(test)

### check
dim(cimac_s1400i_protein_info) # protein info
dim(cimac_s1400i_npx_raw) # data matrix

### labels
colnames(cimac_s1400i_npx_raw) <- cimac_s1400i_protein_info$Assay

### load metadata
cimac_s1400i_samples <- readxl::read_excel("S1400I/Olink/S1400I_Plasma_ID_Completed_01082020_UPLOADED_CIDC_REDACTED.xlsx",sheet = "Samples")
cimac_s1400i_samples <- as.data.frame(cimac_s1400i_samples)
colnames(cimac_s1400i_samples) <- cimac_s1400i_samples[1,]
cimac_s1400i_samples <- cimac_s1400i_samples[-1,]

### inclusion of tables
table(rownames(cimac_s1400i_npx_raw) %in% cimac_s1400i_samples$`Cimac id`)
rownames(cimac_s1400i_npx_raw)[(!rownames(cimac_s1400i_npx_raw) %in% cimac_s1400i_samples$`Cimac id`)]

### separate QC
test <- cimac_s1400i_npx_raw[,97:100]
cimac_s1400i_npx_raw <- cimac_s1400i_npx_raw[,-c(97:100)]

### separate controls
cimac_s1400i_npx_raw_controls <- cimac_s1400i_npx_raw[which((!rownames(cimac_s1400i_npx_raw) %in% cimac_s1400i_samples$`Cimac id`)),]

### separate samples
cimac_s1400i_npx_raw <- cimac_s1400i_npx_raw[which((rownames(cimac_s1400i_npx_raw) %in% cimac_s1400i_samples$`Cimac id`)),]

### check
cimac_s1400i_samples <- cimac_s1400i_samples[ match(rownames(cimac_s1400i_npx_raw), cimac_s1400i_samples$`Cimac id`), ]
identical(rownames(cimac_s1400i_npx_raw) , cimac_s1400i_samples$`Cimac id`)

### compare QCs
test <- test[which((rownames(test) %in% cimac_s1400i_samples$`Cimac id`)),]
identical(rownames(test) , cimac_s1400i_samples$`Cimac id`)

### combine
cimac_s1400i_samples <- cbind(cimac_s1400i_samples,test) 
rm(test)

### convert to numeric
cimac_s1400i_npx_raw <- cimac_s1400i_npx_raw %>% mutate_at(1:ncol(cimac_s1400i_npx_raw), as.numeric)
cimac_s1400i_npx_raw <- t(cimac_s1400i_npx_raw)

### files
dim(cimac_s1400i_samples) # samples data
dim(cimac_s1400i_npx_raw_controls) # control data
dim(cimac_s1400i_npx) # original data
dim(cimac_s1400i_protein_info) # protein info
dim(cimac_s1400i_npx_raw) # data matrix
################################################
save(file="/Users/gonzae34/Documents/projects_gnjatic/U24/rdata_files/cimac_s1400i_raw_data.RData",
     cimac_s1400i_samples,cimac_s1400i_npx_raw_controls,cimac_s1400i_npx,cimac_s1400i_protein_info,cimac_s1400i_npx_raw)
################################################


### 10021
#############################################################################################################################################
#############################################################################################################################################

### Data processing
#############################################################################################################################################
cimac_10021_npx <- readxl::read_excel("10021/Olink/CIMAC_10021_AllSamples_01212020_NPX.xlsx")


### inspect
### dim(cimac_10021_npx) #207 101
### cimac_s1400i_npx[1:10,1:5]
### cimac_s1400i_npx[1:10,1:5] #7
### cimac_s1400i_npx[231:234,1:5] #232

### Prptein names
cimac_10021_protein_info <- as.data.frame( cbind(t(cimac_10021_npx[2:5,]),as.data.frame(t(cimac_10021_npx[206:207,])) ))
colnames(cimac_10021_protein_info) <- cimac_10021_protein_info[1,]
cimac_10021_protein_info <- cimac_10021_protein_info[-1,]

### values
cimac_10021_npx_raw <- as.data.frame(cimac_10021_npx[7:204,])

### duplicated names
test <- cimac_10021_npx_raw$CIMAC_10021_IO[grep("IPC", cimac_10021_npx_raw$CIMAC_10021_IO)] 
test <- paste(1:length(test),test,sep="___")
cimac_10021_npx_raw$CIMAC_10021_IO[grep("IPC", cimac_10021_npx_raw$CIMAC_10021_IO)]  <- test

### duplicated names
test <- cimac_10021_npx_raw$CIMAC_10021_IO[grep("NC", cimac_10021_npx_raw$CIMAC_10021_IO)] 
test <- paste(1:length(test),test,sep="___")
cimac_10021_npx_raw$CIMAC_10021_IO[grep("NC", cimac_10021_npx_raw$CIMAC_10021_IO)]  <- test

### labels
rownames(cimac_10021_npx_raw) <- cimac_10021_npx_raw$CIMAC_10021_IO

### clean
cimac_10021_npx_raw <- cimac_10021_npx_raw[,-1] ; rm(test)

### check
dim(cimac_10021_protein_info) # protein info
dim(cimac_10021_npx_raw) # data matrix

### labels
colnames(cimac_10021_npx_raw) <- cimac_10021_protein_info$Assay


### load metadata
cimac_10021_samples <- readxl::read_excel("10021/Olink/10021_Plasma_IDs_CORRECTIONS_1_22_20.xlsx",sheet = "Samples")
cimac_10021_samples <- as.data.frame(cimac_10021_samples)
colnames(cimac_10021_samples) <- cimac_10021_samples[1,]
cimac_10021_samples <- cimac_10021_samples[-1,]

### inclusion of tables
table(rownames(cimac_10021_npx_raw) %in% cimac_10021_samples$`Cimac id`)
rownames(cimac_10021_npx_raw)[(!rownames(cimac_10021_npx_raw) %in% cimac_10021_samples$`Cimac id`)]

### separate QC
test <- cimac_10021_npx_raw[,97:100]
cimac_10021_npx_raw <- cimac_10021_npx_raw[,-c(97:100)]

### separate controls
cimac_10021_npx_raw_controls <- cimac_10021_npx_raw[which((!rownames(cimac_10021_npx_raw) %in% cimac_10021_samples$ `Cimac id`)),]

### separate samples
cimac_10021_npx_raw <- cimac_10021_npx_raw[which((rownames(cimac_10021_npx_raw) %in% cimac_10021_samples$`Cimac id`)),]

### check
cimac_10021_samples <- cimac_10021_samples[ match(rownames(cimac_10021_npx_raw), cimac_10021_samples$`Cimac id`), ]
identical(rownames(cimac_10021_npx_raw) , cimac_10021_samples$`Cimac id`)

### compare QCs
test <- test[which((rownames(test) %in% cimac_10021_samples$`Cimac id`)),]
identical(rownames(test) , cimac_10021_samples$`Cimac id`)

### combine
cimac_10021_samples <- cbind(cimac_10021_samples,test) 
rm(test)

### convert to numeric
cimac_10021_npx_raw <- cimac_10021_npx_raw %>% mutate_at(1:ncol(cimac_10021_npx_raw), as.numeric)
cimac_10021_npx_raw <- t(cimac_10021_npx_raw)


### files
dim(cimac_10021_samples) # samples data
dim(cimac_10021_npx_raw_controls) # control data
dim(cimac_10021_npx) # original data
dim(cimac_10021_protein_info) # protein info
dim(cimac_10021_npx_raw) # data matrix
################################################
save(file="/Users/gonzae34/Documents/projects_gnjatic/U24/rdata_files/cimac_10021_raw_data.RData",
     cimac_10021_samples,cimac_10021_npx_raw_controls,cimac_10021_npx,cimac_10021_protein_info,cimac_10021_npx_raw)
################################################

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################



### Data composition
#############################################################################################################################################

### s1400i
#################### 
rm(list=setdiff(ls(), ""))
load(file="/Users/gonzae34/Documents/projects_gnjatic/U24/rdata_files/s1400i_analysis.RData")
rm(list=setdiff(ls(), "cimac_s1400i_samples"))

head(cimac_s1400i_samples)

write.csv("/Users/gonzae34/Documents/cimac_s1400i_samples")

my_tmp_df <- table(cimac_s1400i_samples$condition,cimac_s1400i_samples$`Cohort name`)

class(my_tmp_df) <- "matrix"

my_tmp_df <- as.data.frame(my_tmp_df)
my_tmp_df$Timepoint <- rownames(my_tmp_df)
my_tmp_df <- melt(my_tmp_df)

colnames(my_tmp_df)[2]<- "Cohort"

ggplot(data=my_tmp_df) + aes(x=Cohort, y=Timepoint,fill=value) + geom_tile(alpha=0.5) +
        theme_bw() + geom_text(aes(label=value),color='black') + xlab("") + ylab("") +
        #scale_fill_gradient(low = "black", high = "firebrick",space = "Lab", na.value = "grey50", guide = "colourbar",aesthetics = "fill") #+
        scale_fill_viridis(discrete = FALSE) 

### 10021
#################### 
rm(list=setdiff(ls(), ""))
load(file="/Users/gonzae34/Documents/projects_gnjatic/U24/rdata_files/10021_analysis.RData")
rm(list=setdiff(ls(), "cimac_10021_samples"))

head(cimac_10021_samples)

#write.csv("/Users/gonzae34/Documents/cimac_s1400i_samples")

my_tmp_df <- table(cimac_10021_samples$condition,cimac_10021_samples$cohort)

class(my_tmp_df) <- "matrix"

my_tmp_df <- as.data.frame(my_tmp_df)
my_tmp_df$Timepoint <- rownames(my_tmp_df)
my_tmp_df <- melt(my_tmp_df)

colnames(my_tmp_df)[2]<- "Cohort"

ggplot(data=my_tmp_df) + aes(x=Timepoint, y=Cohort,fill=value) + geom_tile(alpha=0.5) +
        theme_bw() + geom_text(aes(label=value),color='black') + xlab("") + ylab("") +
        theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
        #scale_fill_gradient(low = "black", high = "firebrick",space = "Lab", na.value = "grey50", guide = "colourbar",aesthetics = "fill") #+
        scale_fill_viridis(discrete = FALSE) #+

### 10021
#################### 

rm(list=setdiff(ls(), ""))
load("/Users/gonzae34/Documents/projects_gnjatic/U24/U24/E4412_Diefenbach/e4412_analysis.RData")
rm(list=setdiff(ls(), "cimac_e4412_olink_clustergram_meta"))

my_tmp_df <- table(cimac_e4412_olink_clustergram_meta$Cohort,cimac_e4412_olink_clustergram_meta$Collection_Event)

class(my_tmp_df) <- "matrix"
         
my_tmp_df <- as.data.frame(my_tmp_df)
my_tmp_df$Timepoint <- rownames(my_tmp_df)
my_tmp_df <- melt(my_tmp_df)
         
colnames(my_tmp_df)[2]<- "Cohort"
         
ggplot(data=my_tmp_df) + aes(x=Timepoint, y=Cohort,fill=value) + geom_tile(alpha=0.5) +
                 theme_bw() + geom_text(aes(label=value),color='black') + xlab("") + ylab("") +
                 theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
                 #scale_fill_gradient(low = "black", high = "firebrick",space = "Lab", na.value = "grey50", guide = "colourbar",aesthetics = "fill") #+
                 scale_fill_viridis(discrete = FALSE) #+

### Commonly detected markers 
#############################################################################################################################################
rm(list=setdiff(ls(), ""))
#
load(file="/Users/gonzae34/Documents/projects_gnjatic/U24/rdata_files/s1400i_analysis.RData")
rm(list=setdiff(ls(), "olink_de_results_s1400i"))

load(file="/Users/gonzae34/Documents/projects_gnjatic/U24/rdata_files/10021_analysis.RData")
rm(list=setdiff(ls(), c("olink_de_results_s1400i","olink_de_results_10021")))
###lm_olink_all[lm_olink_all$adj.P.Val<0.05,]
olink_de_results_e4412 <- read.csv(file="/Users/gonzae34/Documents/projects_gnjatic/U24/U24/E4412_Diefenbach/olink_e4412_results_table_eegk.csv")


my_list_res <- list( E4412 = unique(sort(olink_de_results_e4412$protein[olink_de_results_e4412$adj.P.Val<0.05])) ,
                     `10021` = gsub("-","\\.",unique(sort(olink_de_results_10021$Marker[olink_de_results_10021$adj.P.Val<0.05])) ) ,
                     S1400I = gsub("-","\\.",unique(sort(olink_de_results_s1400i$Marker[olink_de_results_s1400i$adj.P.Val<0.05])) ) )

library(Vennerable)
library(UpSetR)
library(ComplexHeatmap)

### Venn figure
W <- compute.Venn(Venn(my_list_res))

gp <- VennThemes(W)
gp[["Face"]][["000"]]$fill <-  "white"
gp[["Face"]][["001"]]$fill <-  "white"
gp[["Face"]][["010"]]$fill <-  "white"
gp[["Face"]][["100"]]$fill <-  "white"

plot(W, gp = gp)


####
m = make_comb_mat(my_list_res)
comb_name(m)
set_size(m)
comb_size(m)
extract_comb(m,'110') 
extract_comb(m,'101') # patients exosomes
extract_comb(m,'011') # patients cells
extract_comb(m,'111') # all

rm(list=setdiff(ls(), ""))

#############################################################################################################################################











#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################