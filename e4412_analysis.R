#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
libs<-c("corrplot","ROTS","Hmisc","tcR","tidyr","limma","NMF","tsne","Rtsne" ,"UpSetR",
        "AnnotationDbi","RColorBrewer","papmap","sva","GO.db","fitdistrplus","ff","plyranges",
        "annotables","Rsamtools","GenomicFeatures","ggarrange","pheatmap","rms",
        "dplyr","DESeq2","ggplot2","ggrepel","edgeR","doParallel","ggradar","colorspace",
        "reshape2","rmarkdown","org.Hs.eg.db","treemapify",
        "variancePartition","ggpubr","factoextra","Mfuzz", "universalmotif","ggbio","GenomicRanges",
        "randomForest","caret","mlbench","psych","scales","tibble","RColorBrewer","tidyr","ggplot2","reshape2","circlize",
        "pROC","plotROC","nnet","caTools","MLmetrics","colorspace","Vennerable","enrichR",
        "glmnet","Gviz","randomForestSRC","ggRandomForests","e1071","cowplot","data.table",
        "ROCR","mlbench","caretEnsemble","gbm","rpart","UpSetR","SCENIC","AUCell","RcisTarget","plyr",
        "tidyverse","hrbrthemes","fmsb","colormap","viridis","BSgenome.Hsapiens.UCSC.hg38","ggbeeswarm",
        "umap","Rtsne")
lapply(libs, require, character.only = TRUE)
rm(libs)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dplyr")
setwd("/Users/gonzae34/Documents/projects_gnjatic/U24/U24/E4412_Diefenbach/")
# setwd("D:/Documents/Tolerized_TCells/tolerized_scRNAseq/")
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
### Reading NPX data from excel
cimac_e4412_olink <- as.data.frame(readxl::read_excel("Olink Data/CIMAC E4412_NPX.xlsx"))
cimac_e4412_olink_feature_data <- cimac_e4412_olink[1:5,]
cimac_e4412_olink_feature_data[1:5,1:5]
rownames(cimac_e4412_olink_feature_data) <- cimac_e4412_olink_feature_data$`CIMAC E4412`
cimac_e4412_olink_feature_data <- cimac_e4412_olink_feature_data[,-1]
###
cimac_e4412_olink <- cimac_e4412_olink[-c(1:5),]
ii <- duplicated(cimac_e4412_olink[,1])
cimac_e4412_olink[,1][ii]
cimac_e4412_olink <- cimac_e4412_olink[!ii,]
rownames(cimac_e4412_olink) <- cimac_e4412_olink[,1]
cimac_e4412_olink <- cimac_e4412_olink[,-1]
###
colnames(cimac_e4412_olink) <- as.character(cimac_e4412_olink_feature_data[3,])
###
ii <- rownames(cimac_e4412_olink)
cimac_e4412_olink <- data.frame(lapply(cimac_e4412_olink, function(x) as.numeric(as.character(x))))
rownames(cimac_e4412_olink) <- ii
###
str(cimac_e4412_olink)
###
save(file="cimac_olink_e4412_raw.RData",cimac_e4412_olink,cimac_e4412_olink_feature_data)
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
### read the metadata file from Cytof; copied from astrolabe into a tsv
cimac_e4412_metadata <- read.table(file="metadata_samples.tsv",sep="\t")
### rename the columns properly
colnames(cimac_e4412_metadata) <- c("File_Name","Sample_Name","Participant_ID","Sample_ID","CIMAC_ID","Timepoint","Cohort","Sample_Batch",
                                    "Veritas_batch","Baseline_all_cohorts","Pre_Day_1_Cycle_2_all_cohorts","Restaging_all_cohorts","End_of_treatment_all_cohorts","Arm_A_all_time_points",
                                    "Arm_B_all_time_points","Arm_C_all_time_points","Arm_D_all_time_points","Arm_E_all_time_points","Arm_F_all_time_points","Arm_G_all_time_points",
                                    "Arm_H_all_time_points","Arm_I_all_time_points")
###
### Process the names into a single variable 
###
my_list <- list( a=which(!is.na(cimac_e4412_metadata$Arm_A_all_time_points)),
                 b=which(!is.na(cimac_e4412_metadata$Arm_B_all_time_points)),
                 c=which(!is.na(cimac_e4412_metadata$Arm_C_all_time_points)),
                 d=which(!is.na(cimac_e4412_metadata$Arm_D_all_time_points)),
                 e=which(!is.na(cimac_e4412_metadata$Arm_E_all_time_points)),
                 f=which(!is.na(cimac_e4412_metadata$Arm_F_all_time_points)),
                 g=which(!is.na(cimac_e4412_metadata$Arm_G_all_time_points)),
                 h=which(!is.na(cimac_e4412_metadata$Arm_H_all_time_points)),
                 i=which(!is.na(cimac_e4412_metadata$Arm_I_all_time_points)) )
### Reduce(intersect,my_list)
###
cimac_e4412_metadata <- cimac_e4412_metadata[ ! is.na(cimac_e4412_metadata$CIMAC_ID) , ]
###
temp_metadata <- data.frame(  pos =  c(which(!is.na(cimac_e4412_metadata$Arm_A_all_time_points)),
                                       which(!is.na(cimac_e4412_metadata$Arm_B_all_time_points)),
                                       which(!is.na(cimac_e4412_metadata$Arm_C_all_time_points)),
                                       which(!is.na(cimac_e4412_metadata$Arm_D_all_time_points)),
                                       which(!is.na(cimac_e4412_metadata$Arm_E_all_time_points)),
                                       which(!is.na(cimac_e4412_metadata$Arm_F_all_time_points)),
                                       which(!is.na(cimac_e4412_metadata$Arm_G_all_time_points)),
                                       which(!is.na(cimac_e4412_metadata$Arm_H_all_time_points)),
                                       which(!is.na(cimac_e4412_metadata$Arm_I_all_time_points)))
                              ,
                              names =  c( cimac_e4412_metadata$Arm_A_all_time_points[ !is.na(cimac_e4412_metadata$Arm_A_all_time_points) ],
                                          cimac_e4412_metadata$Arm_B_all_time_points[ !is.na(cimac_e4412_metadata$Arm_B_all_time_points) ],
                                          cimac_e4412_metadata$Arm_C_all_time_points[ !is.na(cimac_e4412_metadata$Arm_C_all_time_points) ],
                                          cimac_e4412_metadata$Arm_D_all_time_points[ !is.na(cimac_e4412_metadata$Arm_D_all_time_points) ],
                                          cimac_e4412_metadata$Arm_E_all_time_points[ !is.na(cimac_e4412_metadata$Arm_E_all_time_points) ],
                                          cimac_e4412_metadata$Arm_F_all_time_points[ !is.na(cimac_e4412_metadata$Arm_F_all_time_points) ],
                                          cimac_e4412_metadata$Arm_G_all_time_points[ !is.na(cimac_e4412_metadata$Arm_G_all_time_points) ],
                                          cimac_e4412_metadata$Arm_H_all_time_points[ !is.na(cimac_e4412_metadata$Arm_H_all_time_points) ],
                                          cimac_e4412_metadata$Arm_I_all_time_points[ !is.na(cimac_e4412_metadata$Arm_I_all_time_points) ])    )
###
temp_metadata <- temp_metadata [ order(temp_metadata$pos,decreasing = FALSE) , ]
cimac_e4412_metadata$time_points <- temp_metadata$names
rm(temp_metadata,my_list,ii)
### clean temporal files
# cimac_e4412_metadata$time_points
# cimac_e4412_metadata$Sample_Batch
# table(rownames(cimac_e4412_olink) %in% cimac_e4412_metadata$CIMAC_ID)
###
# table( rownames(cimac_e4412_olink) %in% "C29ZS2M01.01" )
# grep( "C29ZS2" , cimac_e4412_metadata$CIMAC_ID ,value = TRUE)
# cimac_e4412_metadata$End_of_treatment_all_cohorts[grep( "C29ZS2" , cimac_e4412_metadata$CIMAC_ID)]
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
### read the csv files from clustergrammer
generic_var <- read.csv(file="ELISA CIMAC E4412/2_npx_CIMAC_E4412_Clustergrammer.csv",fill = TRUE)
###
table( colnames(generic_var) %in% rownames(cimac_e4412_olink) )
table( colnames(generic_var) %in% cimac_e4412_metadata$CIMAC_ID )
###
my_proteins <- generic_var[,1]
###
colnames(generic_var) [ ! colnames(generic_var) %in% rownames(cimac_e4412_olink) ]
###
generic_var <- generic_var[ , colnames(generic_var) %in% rownames(cimac_e4412_olink) ]
###
cimac_e4412_olink_clustergram <- generic_var
###
cimac_e4412_olink_clustergram_meta <- generic_var[1:7,]
cimac_e4412_olink_clustergram <- data.frame(lapply(cimac_e4412_olink_clustergram, function(x) as.numeric(as.character(x))))
###
cimac_e4412_olink_clustergram <- cimac_e4412_olink_clustergram[ -c(1:7) , ]
###
rm(generic_var)
###
rownames(cimac_e4412_olink_clustergram) <- my_proteins[-c(1:7)]
###
table(is.na(cimac_e4412_olink_clustergram))
###
table( rownames(cimac_e4412_olink) %in% colnames(cimac_e4412_olink_clustergram) )
cimac_e4412_olink <- cimac_e4412_olink [ rownames(cimac_e4412_olink) %in% colnames(cimac_e4412_olink_clustergram) , ]
cimac_e4412_olink <- t(cimac_e4412_olink)
###
rownames(cimac_e4412_olink_clustergram) <- gsub("Protein: ","",rownames(cimac_e4412_olink_clustergram) )
rownames(cimac_e4412_olink_clustergram) <- gsub(" ",".",rownames(cimac_e4412_olink_clustergram) )
####
cimac_e4412_olink <-cimac_e4412_olink [ ! rownames(cimac_e4412_olink) %in% c("Plate.ID","QC.Warning","QC.Deviation.from.median","QC.Deviation.from.median.1") , ]
###
grep( "Ctr", rownames(cimac_e4412_olink_clustergram) ,value=TRUE )
grep( "Ctr", rownames(cimac_e4412_olink) ,value=TRUE )
###
cimac_e4412_olink_clustergram <- cimac_e4412_olink_clustergram [ -grep( "Ctr", rownames(cimac_e4412_olink_clustergram)) ,]
cimac_e4412_olink <- cimac_e4412_olink [ -grep( "Ctr", rownames(cimac_e4412_olink)) ,]
###
rownames(cimac_e4412_olink_clustergram)
rownames(cimac_e4412_olink)
###
cimac_e4412_olink_clustergram <- cimac_e4412_olink_clustergram[ , order(colnames(cimac_e4412_olink_clustergram)) ]
cimac_e4412_olink <- cimac_e4412_olink[ , order(colnames(cimac_e4412_olink)) ]
###
identical( colnames(cimac_e4412_olink),colnames(cimac_e4412_olink_clustergram) )
###
rownames(cimac_e4412_olink_clustergram) <- gsub("-",".",rownames(cimac_e4412_olink_clustergram) )
rownames(cimac_e4412_olink_clustergram) <- gsub("/",".",rownames(cimac_e4412_olink_clustergram) )
###
### both files for olink data. Raw and some convertion (quantile said Sacha)
###
identical(rownames(cimac_e4412_olink),rownames(cimac_e4412_olink_clustergram))
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
###
### Pre processing
### Olink data for e4412 ; Verify they are not NAs
###
table(is.na(cimac_e4412_olink_clustergram))
table(is.na(cimac_e4412_olink))
###
rownames(cimac_e4412_olink_clustergram_meta) <- c("Arm","Cohort","Collection_Event","Sample_Type","QC","Plate","ParticipantID")
###
cimac_e4412_olink_clustergram_meta <- t(cimac_e4412_olink_clustergram_meta)
###
cimac_e4412_olink_clustergram_meta <- as.data.frame(cimac_e4412_olink_clustergram_meta)
cimac_e4412_olink_clustergram_meta$Arm <- gsub("Arm: ","", cimac_e4412_olink_clustergram_meta$Arm )
cimac_e4412_olink_clustergram_meta$Cohort <- gsub("Cohort: ","", cimac_e4412_olink_clustergram_meta$Cohort )
cimac_e4412_olink_clustergram_meta$Cohort <- gsub("Cohort_","", cimac_e4412_olink_clustergram_meta$Cohort )
cimac_e4412_olink_clustergram_meta$Collection_Event <- gsub("Collection_Event: ", "", cimac_e4412_olink_clustergram_meta$Collection_Event)
cimac_e4412_olink_clustergram_meta$Sample_Type <- gsub("Sample type: ", "", cimac_e4412_olink_clustergram_meta$Sample_Type)
cimac_e4412_olink_clustergram_meta$QC <- gsub("QC: ", "", cimac_e4412_olink_clustergram_meta$QC)
cimac_e4412_olink_clustergram_meta$Plate <- gsub("Plate: ", "", cimac_e4412_olink_clustergram_meta$Plate)
cimac_e4412_olink_clustergram_meta$ParticipantID <- gsub("Participant_ID: ", "", cimac_e4412_olink_clustergram_meta$ParticipantID)
cimac_e4412_olink_clustergram_meta$Sample_ID <- rownames(cimac_e4412_olink_clustergram_meta)
###  
head(cimac_e4412_olink_clustergram_meta)
#############################################################################################################################################
### for( i in 1:ncol(cimac_e4412_olink_clustergram_meta)){ print(table(cimac_e4412_olink_clustergram_meta[,i])) }
### cimac_e4412_olink_clustergram_meta$Sample_Type ### is the clinical trial ID
#############################################################################################################################################

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
###
### VARIANCE ANALYSIS START
###
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

#############################################################################################################################################
### For Mac with 12 Cores
libs<-c("variancePartition","doParallel")
lapply(libs, require, character.only = TRUE)
registerDoParallel(makeCluster(10))
rm(libs)
#############################################################################################################################################
###
### inspect the data
### Raw:
#cimac_e4412_olink_clustergram[1:5,1:5]
### ***Normalized:
#cimac_e4412_olink[1:5,1:5]
###
### inspect the files
### metadata
# cimac_e4412_olink_clustergram_meta
cimac_e4412_olink_clustergram_meta <- cimac_e4412_olink_clustergram_meta [ match( colnames(cimac_e4412_olink), cimac_e4412_olink_clustergram_meta$Sample_ID ) , ]
identical( colnames(cimac_e4412_olink_clustergram), colnames(cimac_e4412_olink) )
identical( cimac_e4412_olink_clustergram_meta$Sample_ID, colnames(cimac_e4412_olink) )
###
table(cimac_e4412_olink_clustergram_meta$Arm,cimac_e4412_olink_clustergram_meta$Cohort)
table(cimac_e4412_olink_clustergram_meta$Arm,cimac_e4412_olink_clustergram_meta$Plate)
table(cimac_e4412_olink_clustergram_meta$Arm,cimac_e4412_olink_clustergram_meta$Collection_Event)
table(cimac_e4412_olink_clustergram_meta$Cohort,cimac_e4412_olink_clustergram_meta$Collection_Event)
cimac_e4412_olink_clustergram_meta$ParticipantID
#############################################################################################################################################
### Form for model
form <- ~ (1|ParticipantID) + (1|Cohort) + (1|Collection_Event) + (1|Plate)
#############################################################################################################################################
cimac_e4412_olink[is.na(cimac_e4412_olink)] <- 0
### run the model
variance_exprs_matrix <- fitExtractVarPartModel(cimac_e4412_olink, form, cimac_e4412_olink_clustergram_meta)
###
#############################################################################################################################################
pdf(file="figures/VP_npx_ee4412_olink.pdf",width = 4, height = 4)
plotVarPart( sortCols( variance_exprs_matrix ) ) + ggtitle( 'NPX Data') + theme(axis.text.x = element_text(angle = 60,hjust = 1))
dev.off()
#############################################################################################################################################
### LEGACY RESULTS
#############################################################################################################################################
#cimac_e4412_olink_clustergram[is.na(cimac_e4412_olink_clustergram)] <- 0
### run the model
#variance_exprs_matrix_raw <- fitExtractVarPartModel(cimac_e4412_olink_clustergram, form, cimac_e4412_olink_clustergram_meta)
###
#############################################################################################################################################
#pdf(file="figures/VP_raw_ee4412_olink.pdf",width = 4, height = 4)
#plotVarPart( sortCols( variance_exprs_matrix_raw ) ) + ggtitle( 'Raw Data') + theme(axis.text.x = element_text(angle = 60,hjust = 1))
#dev.off()
#dev.off()
#############################################################################################################################################

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

summary(colSums(cimac_e4412_olink_clustergram))
summary(log(colSums(cimac_e4412_olink_clustergram)))
summary(colSums(cimac_e4412_olink))

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
### DIMENSIONALITY REDUCTION
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
###
mydata <- prcomp( t(cimac_e4412_olink[,]) , scale=TRUE, center = TRUE ) 
###
fviz_eig(mydata, addlabels = T, barfill = "lightsteelblue3", barcolor = "lightsteelblue3") + theme_bw() + ggtitle("")
###
pca_mts <- as.data.frame(mydata$x)
pca_mts$Sample_ID <- rownames(pca_mts)
###
umap_clustering <- umap(t(cimac_e4412_olink[,]))
pca_mts$UMAP1 <- umap_clustering$layout[,1]
pca_mts$UMAP2 <- umap_clustering$layout[,2]
###
tsne_clustering <- Rtsne(t(cimac_e4412_olink), dims = 2, initial_dims = 30,
                         perplexity = 15, theta = 1, check_duplicates = TRUE,
                         pca = TRUE, partial_pca = FALSE, max_iter = 5000,
                         normalize = TRUE, momentum = 0.5, final_momentum = 0.8, eta = 200, exaggeration_factor = 12, num_threads = 1)
###      
pca_mts$tsne1 <- tsne_clustering$Y[,1]
pca_mts$tsne2 <- tsne_clustering$Y[,2]
###
### over PCA
umap_clustering <- umap(pca_mts[,1:30])
pca_mts$UMAP1_s <- umap_clustering$layout[,1]
pca_mts$UMAP2_s <- umap_clustering$layout[,2]
### over PCA
tsne_clustering <- Rtsne(pca_mts[,1:30], dims = 2, initial_dims = 30,
                         perplexity = 15, theta = 1, check_duplicates = FALSE,
                         pca = FALSE, partial_pca = FALSE, max_iter = 5000,
                         normalize = TRUE, momentum = 0.5, final_momentum = 0.8, eta = 200, 
                         exaggeration_factor = 12, num_threads = 1)
### 
pca_mts$tsne1_s <- tsne_clustering$Y[,1]
pca_mts$tsne2_s <- tsne_clustering$Y[,2]

### merge to add metadata columns
pca_mts <- merge(pca_mts, cimac_e4412_olink_clustergram_meta, by="Sample_ID")

#############################################################################################################################################
pdf(file="figures/reductions_all_data_ee4412_olink.pdf",width = 14, height = 10)
ggarrange( ggplot(pca_mts, aes(PC1, PC2, color = Condition)) + theme_bw() + geom_point(alpha=0.69) ,
           ggplot(pca_mts, aes(UMAP1, UMAP2, color = Condition)) + theme_bw() + geom_point(alpha=0.69) , 
           fviz_eig(mydata, addlabels = T, barfill = "lightsteelblue3", barcolor = "lightsteelblue3") + theme_bw() + ggtitle(""), 
           ggplot(pca_mts, aes(tsne1, tsne2, color = Condition)) + theme_bw() + geom_point(alpha=0.69) , 
           ncol=2,nrow=2 )
dev.off()
#############################################################################################################################################

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
### Linear Model Pipeline
# design <- model.matrix(~ 0 + Cohort + Collection_Event + Plate , data = cimac_e4412_olink_clustergram_meta )
# cimac_e4412_olink_clustergram_meta$Collection_Event %in% "Baseline"
cimac_e4412_olink_clustergram_meta$Condition <- paste(cimac_e4412_olink_clustergram_meta$Collection_Event, cimac_e4412_olink_clustergram_meta$Cohort,sep="_")
###
identical( cimac_e4412_olink_clustergram_meta$Sample_ID, colnames(cimac_e4412_olink) )
###
design <- model.matrix(~ 0 + Condition + Plate , data = cimac_e4412_olink_clustergram_meta )
###
table(duplicated(cimac_e4412_olink_clustergram_meta$ParticipantID))
###
colnames(design) <- c("C1_B","C2_B","C3_B",
                      "C1_OS","C2_OS","C3_OS",
                      "C1_Pre","C2_Pre","C3_Pre",
                      "C1_Res","C2_Res","C3_Res","Plaet")
### trying with normalized
block2 <- as.numeric(as.factor( cimac_e4412_olink_clustergram_meta$ParticipantID ))
### function estimates the correlation between repeated observations
dupcor2 <- duplicateCorrelation(cimac_e4412_olink, design, block=block2)
### model accounting for replicates
vfit <- lmFit(cimac_e4412_olink, design, block=block2, correlation=dupcor2$consensus)
efit <- eBayes(vfit)
### summary(decideTests(efit))
#############################################################################################################################################
contr.matrix <- makeContrasts( "C1vsC2_B" = C1_B - C2_B,
                               "C1vsC3_B" = C1_B - C3_B,
                               "C2vsC3_B" = C2_B - C3_B,
                               
                               "C1vsC2_OS" = C1_OS - C2_OS,
                               "C1vsC3_OS" = C1_OS - C3_OS,
                               "C2vsC3_OS" = C2_OS - C3_OS,
                               
                               "C1vsC2_Pre" = C1_Pre - C2_Pre,
                               "C1vsC3_Pre" = C1_Pre - C3_Pre,
                               "C2vsC3_Pre" = C2_Pre - C3_Pre,
                               
                               "C1vsC2_Res" = C1_Res - C2_Res,
                               "C1vsC3_Res" = C1_Res - C3_Res,
                               "C2vsC3_Res" = C2_Res - C3_Res,
                               
                               "C1vsC1_Base_pre" = C1_B - C1_Pre,
                               "C2vsC2_Base_pre" = C2_B - C2_Pre,
                               "C3vsC3_Base_pre" = C3_B - C3_Pre,
                               
                               "C1vsC1_Base_OS" = C1_B - C1_OS,
                               "C2vsC2_Base_OS" = C2_B - C2_OS,
                               "C3vsC3_Base_OS" = C3_B - C3_OS,
                               
                               "C1vsC1_Base_Res" = C1_B - C1_Res,
                               "C2vsC2_Base_Res" = C2_B - C2_Res,
                               "C3vsC3_Base_Res" = C3_B - C3_Res,
                               
                               "C1vsC1_Pre_OS" = C1_Pre - C1_OS,
                               "C2vsC2_Pre_OS" = C2_Pre - C2_OS,
                               "C3vsC3_Pre_OS" = C3_Pre - C3_OS,
                               
                               "C1vsC1_Res_OS" = C1_Res - C1_OS,
                               "C2vsC2_Res_OS" = C1_Res - C1_OS,
                               "C3vsC3_Res_OS" = C1_Res - C1_OS,
                               
                               "C1vsC1_Pre_Res" = C1_Pre - C1_Res,
                               "C2vsC2_Pre_Res" = C2_Pre - C2_Res,
                               "C3vsC3_Pre_Res" = C3_Pre - C3_Res,
                               ### Average baseline - condition
                               "BaselinevsC1_Pre" = (C1_B+C2_B+C3_B)/3 - C1_Pre,
                               "BaselinevsC2_Pre" = (C1_B+C2_B+C3_B)/3 - C2_Pre,
                               "BaselinevsC3_Pre" = (C1_B+C2_B+C3_B)/3 - C3_Pre,
                               
                               "BaselinevsC1_Res" = (C1_B+C2_B+C3_B)/3 - C1_Res,
                               "BaselinevsC2_Res" = (C1_B+C2_B+C3_B)/3 - C2_Res,
                               "BaselinevsC3_Res" = (C1_B+C2_B+C3_B)/3 - C3_Res,
                               
                               "BaselinevsC1_Res" = (C1_B+C2_B+C3_B)/3 - C1_OS,
                               "BaselinevsC2_Res" = (C1_B+C2_B+C3_B)/3 - C2_OS,
                               "BaselinevsC3_Res" = (C1_B+C2_B+C3_B)/3 - C3_OS,
                               
                               levels = colnames(design) )
###
vfit <- lmFit(cimac_e4412_olink, design, block=block2, correlation=dupcor2$consensus, method = "robust")
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE)
summary(decideTests(efit))
#############################################################################################################################################
### Summarize into a single data.frame object
lm_olink <- list()
for ( i in 1:length( colnames(summary(decideTests(efit))) ) ) {
  lm_olink[[i]] <- topTable(efit,coef= colnames(summary(decideTests(efit)))[i] ,n=Inf)  
  lm_olink[[i]]$Comparison <- colnames(summary(decideTests(efit)))[i]
  lm_olink[[i]]$protein <- rownames(lm_olink[[i]])
}

lm_olink_all <- do.call(rbind,lm_olink)
###
dim(lm_olink_all)
###
table(lm_olink_all$adj.P.Val<0.05) ## 99
###
# lm_olink_all[ lm_olink_all$adj.P.Val<0.05 , ]
###
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
###
### Diff Expr. Summary Figure. 
###
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
stacked_bar_plot <- lm_olink_all 
stacked_bar_plot$logFC_dir <- stacked_bar_plot$logFC > 0 
stacked_bar_plot$logFC_dir[stacked_bar_plot$logFC > 0] <- "Up"
stacked_bar_plot$logFC_dir[stacked_bar_plot$logFC < 0] <- "Down"
stacked_bar_plot$sig <- stacked_bar_plot$adj.P.Val<0.05
#############################################################################################################################################
stacked_bar_plot <- data.frame(table(stacked_bar_plot$Comparison, stacked_bar_plot$logFC_dir, stacked_bar_plot$sig))
colnames(stacked_bar_plot) <- c('Comparison', 'logFC_dir','sig', 'count')
#############################################################################################################################################
regional_summary <- stacked_bar_plot %>%  filter(sig=='TRUE') %>% dplyr::select(-sig)
regional_summary$count[regional_summary$logFC_dir=='Down'] <- 0 - regional_summary$count[regional_summary$logFC_dir=='Down']
#############################################################################################################################################
regional_summary$Comparison <- as.character(regional_summary$Comparison) 
#############################################################################################################################################
pdf(file="figures/linear_model_all_data_ee4412_olink.pdf",width = 4, height = 6)
ggplot(regional_summary) + 
  aes(x=count, y=Comparison, fill=logFC_dir, label=abs(count)) + 
  geom_bar(stat='identity', width=0.5) +
  geom_text(data=regional_summary[regional_summary$count != 0,], show.legend=FALSE, color='black') +
  theme_classic() + theme(axis.text.y = element_text(face='bold'),
                          axis.text.x = element_text(face='bold'),
                          axis.title = element_text(face = "bold"),
                          legend.title= element_text(face='bold'),
                          strip.text=element_text(face='bold'),
                          plot.title=element_text(face='bold')) + 
  scale_fill_manual(values=c('firebrick', 'steelblue'), name='LogFC Direction') +
  labs(x='no. DGPs (FDR < 0.05)', y='Contrast') + theme(legend.position="bottom")
dev.off()
#############################################################################################################################################

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

#############################################################################################################################################
### test to see if PCA further separate the samples
ii <- rownames(cimac_e4412_olink) %in% lm_olink_all$protein[lm_olink_all$adj.P.Val<0.05]
###
mydata <- prcomp( t(cimac_e4412_olink[ii,]) , scale=TRUE, center = TRUE ) 
###
fviz_eig(mydata, addlabels = T, barfill = "lightsteelblue3", barcolor = "lightsteelblue3") + theme_bw() + ggtitle("") + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=rel(4.8)) )
###
pca_mts <- as.data.frame(mydata$x)
pca_mts$Sample_ID <- rownames(pca_mts)
pca_mts <- merge(pca_mts, cimac_e4412_olink_clustergram_meta, by="Sample_ID")
###
# pdf(file="figures/pca_1sample_raw_somalogic_data_matrix_batch_corrected_problem_samples.pdf",width = 8,height = 6)
# Condition, Cohort, Collection_Event
ggplot(pca_mts, aes(PC1, PC2, color=Condition)) + 
  geom_point() +
  theme_bw() + ggtitle('') + 
  xlab( paste( format(round(get_eigenvalue(mydata)$variance.percent[1], 2), nsmall = 2),"%",sep="" ) ) + 
  ylab( paste( format(round(get_eigenvalue(mydata)$variance.percent[2], 2), nsmall = 2),"%",sep="" ) ) #+
#geom_label_repel(data=pca_mts[pca_mts$Subject_ID %in% my_problems$`Patient ID`[9],],aes(label=SampleId))
# dev.off()
#############################################################################################################################################

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
###
### Examples for presentation
###
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
table(lm_olink_all$Comparison[lm_olink_all$adj.P.Val<0.05])

lm_olink_all$Comparison[lm_olink_all$adj.P.Val<0.05]

test <- lm_olink_all[grep("Baseline",lm_olink_all$Comparison),]
test <- test[test$adj.P.Val<0.05,]

#############################################################################################################################################
generic_var <- melt(cimac_e4412_olink[,])
colnames(generic_var)<- c("protein","Sample_ID","value")
generic_var <- merge(generic_var,cimac_e4412_olink_clustergram_meta,by="Sample_ID")
###
### order the variables using a defined order
generic_var$Condition <- factor(generic_var$Condition, levels = c( "Baseline_1_BV+ipi","Pre_D1_C2_1_BV+ipi","Restaging_1_BV+ipi","Off_study_1_BV+ipi",
                                                                   "Baseline_2_BV+nivo","Pre_D1_C2_2_BV+nivo","Restaging_2_BV+nivo","Off_study_2_BV+nivo",
                                                                   "Baseline_3_BV+nivo+ipi","Pre_D1_C2_3_BV+nivo+ipi","Restaging_3_BV+nivo+ipi","Off_study_3_BV+nivo+ipi" ))
#############################################################################################################################################
my_comparisons <- list( c("Baseline_1_BV+ipi","Pre_D1_C2_1_BV+ipi"),
                        c("Baseline_1_BV+ipi","Restaging_1_BV+ipi"),
                        c("Baseline_1_BV+ipi","Off_study_1_BV+ipi"),
                        c("Restaging_1_BV+ipi","Off_study_1_BV+ipi"),
                        
                        c("Baseline_2_BV+nivo","Pre_D1_C2_2_BV+nivo"),
                        c("Baseline_2_BV+nivo","Restaging_2_BV+nivo"),
                        c("Baseline_2_BV+nivo","Off_study_2_BV+nivo"),
                        c("Restaging_2_BV+nivo","Off_study_2_BV+nivo"),
                        
                        c("Baseline_3_BV+nivo+ipi","Pre_D1_C2_3_BV+nivo+ipi"),
                        c("Baseline_3_BV+nivo+ipi","Restaging_3_BV+nivo+ipi"),
                        c("Baseline_3_BV+nivo+ipi","Off_study_3_BV+nivo+ipi"),
                        c("Restaging_3_BV+nivo+ipi","Off_study_3_BV+nivo+ipi"))
#############################################################################################################################################
table( levels(generic_var$Condition) %in% unlist(my_comparisons) )
levels(generic_var$Condition)[! levels(generic_var$Condition) %in% unlist(my_comparisons)]
#############################################################################################################################################

#############################################################################################################################################
### Fix ordering while vertical plot
generic_var$Condition <- fct_rev(generic_var$Condition)

pdf(file="figures/examples_olink_diff_expr.pdf",width = 6, height = 4)

ii <- generic_var$protein %in% "PDCD1"

ggboxplot(data=generic_var[ii,], 
          x="Condition", y="value",fill="Cohort") + theme_minimal() +
  #geom_boxplot(width=0.1) + facet_wrap(~protein,ncol=4) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("PDCD1") + labs(x = "", y = "NPX Values") + theme(legend.position="none")

ii <- generic_var$protein %in% "CCL17"
ggboxplot(data=generic_var[ii,], 
          x="Condition", y="value",fill="Cohort") + theme_minimal() +
  #geom_boxplot(width=0.1) + facet_wrap(~protein,ncol=4) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("CCL17") + labs(x = "", y = "NPX Values") + theme(legend.position="none")

ii <- generic_var$protein %in% "ANGPT2"

ggboxplot(data=generic_var[ii,], 
          x="Condition", y="value",fill="Cohort") + theme_minimal() +
  #geom_boxplot(width=0.1) + facet_wrap(~protein,ncol=4) + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = TRUE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) + 
  coord_flip()+ 
  ggtitle("ANGPT2") + labs(x = "", y = "NPX Values") + theme(legend.position="none")

ii <- generic_var$protein %in% "IL13" & generic_var$value < 5

ggboxplot(data=generic_var[ii,], 
          x="Condition", y="value",fill="Cohort") + theme_minimal() +
  #geom_boxplot(width=0.1) + facet_wrap(~protein,ncol=4) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("IL13") + labs(x = "", y = "NPX Values") + theme(legend.position="none")

ii <- generic_var$protein %in% "CXCL13" #& generic_var$value < 5

ggboxplot(data=generic_var[ii,], 
          x="Condition", y="value",fill="Cohort") + theme_minimal() +
  #geom_boxplot(width=0.1) + facet_wrap(~protein,ncol=4) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("CXCL13") + labs(x = "", y = "NPX Values") + theme(legend.position="none")

ii <- generic_var$protein %in% "IL10" #& generic_var$value < 5

ggboxplot(data=generic_var[ii,], 
          x="Condition", y="value",fill="Cohort") + theme_minimal() +
  #geom_boxplot(width=0.1) + facet_wrap(~protein,ncol=4) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("IL10") + labs(x = "", y = "NPX Values") + theme(legend.position="none")


ii <- generic_var$protein %in% "CXCL5" #& generic_var$value < 5

ggboxplot(data=generic_var[ii,], 
          x="Condition", y="value",fill="Cohort") + theme_minimal() +
  #geom_boxplot(width=0.1) + facet_wrap(~protein,ncol=4) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("CXCL5") + labs(x = "", y = "NPX Values") + theme(legend.position="none")

dev.off()

ii <- generic_var$protein %in% "GZMA"

ggboxplot(data=generic_var[ii,], 
          x="Condition", y="value",fill="Cohort") + theme_minimal() +
  #geom_boxplot(width=0.1) + facet_wrap(~protein,ncol=4) + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = TRUE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) + 
  coord_flip()+ 
  ggtitle("GZMA") + labs(x = "", y = "NPX Values") + theme(legend.position="none")

ii <- generic_var$protein %in% "LAG3"

ggboxplot(data=generic_var[ii,], 
          x="Condition", y="value",fill="Cohort") + theme_minimal() +
  #geom_boxplot(width=0.1) + facet_wrap(~protein,ncol=4) + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = TRUE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) + 
  coord_flip()+ 
  ggtitle("LAG3") + labs(x = "", y = "NPX Values") + theme(legend.position="none")

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################






#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
###
### Many violin plots, for all markers with at least 1 sig in the linear model
###
ii <- rownames(cimac_e4412_olink) %in% lm_olink_all$protein[lm_olink_all$adj.P.Val<0.05]
###
my_plot  <- list()
for( i in 1:length(rownames(cimac_e4412_olink)[ii])) {
  ix <- generic_var$protein %in% rownames(cimac_e4412_olink)[ii][i]
  my_plot[[i]] <- ggviolin(data=generic_var[ix,], 
                           x="Condition", y="value", color="Cohort",fill="Cohort") + theme_minimal() +
    geom_boxplot(width=0.1) + facet_wrap(~protein,ncol=4) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
    theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
    ggtitle("") + labs(x = "", y = "Values") + theme(legend.position="none")
}
###
pdf(file="figures/violinplots_all_data_ee4412_olink.pdf",width = 8, height = 6)
my_plot
dev.off()
###

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
###
### Composed HEATMAP for OLINK
###
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

### TEST START
#############################################################################################################################################
###
ii <- rownames(cimac_e4412_olink) %in% lm_olink_all$protein[ lm_olink_all$adj.P.Val < 0.05 ]
###
annotation_col <- data.frame(
  Collection = cimac_e4412_olink_clustergram_meta$Collection_Event ,
  Cohort = cimac_e4412_olink_clustergram_meta$Cohort )
###
rownames(annotation_col) <- cimac_e4412_olink_clustergram_meta$Sample_ID
###
dim(cimac_e4412_olink_clustergram_meta)
dim(cimac_e4412_metadata)
### table(cimac_e4412_metadata$Sample_ID %in% cimac_e4412_olink_clustergram_meta$Sample_ID)
pheatmap( cimac_e4412_olink[ii,] , 
          main = "",
          color = colorRampPalette(c("steelblue","white","firebrick"))(255),         
          cluster_cols = T, 
          cluster_rows = T, 
          show_colnames = FALSE,
          border_color = FALSE, 
          show_rownames = FALSE, 
          annotation_col = annotation_col,
          clustering_distance_rows ="euclidean", 
          scale="row")
###
#############################################################################################################################################
### TEST END
###
### START HEATMAP FILES
###
generic_var <- lm_olink_all[lm_olink_all$adj.P.Val<0.05 , ]
###
annotation_col <- data.frame(
  Collection = cimac_e4412_olink_clustergram_meta$Collection_Event ,
  Cohort = cimac_e4412_olink_clustergram_meta$Cohort )
###
rownames(annotation_col) <- cimac_e4412_olink_clustergram_meta$Sample_ID
###
ii <- rownames(cimac_e4412_olink) %in% generic_var$protein
###
annotation_row <- as.matrix(table(generic_var$protein,generic_var$Comparison))
class(annotation_row) <- "matrix"
annotation_row <- as.data.frame(annotation_row)
###
annotation_colors <- list( Cohort = c( `1_BV+ipi`= brewer.pal(7, "Dark2")[1], 
                                       `2_BV+nivo` = brewer.pal(7, "Dark2")[2], 
                                       `3_BV+nivo+ipi` = brewer.pal(7, "Dark2")[3] ) , 
                           Collection = c( Baseline = "grey30",
                                           Pre_D1_C2 = "grey50",
                                           Restaging = "grey70",
                                           Off_study = "grey90"),
                           C1vsC1_Base_OS = c( `1`="white", `0` = "black"),
                           C1vsC1_Base_pre = c( `1`="white", `0` = "black"),
                           C1vsC1_Base_Res = c( `1`="white", `0` = "black"), 
                           C1vsC1_Pre_OS = c( `1`="white", `0` = "black"),  
                           C1vsC2_B = c( `1`="white", `0` = "black"), 
                           C1vsC2_OS = c( `1`="white", `0` = "black"), 
                           C1vsC2_Pre = c( `1`="white", `0` = "black"), 
                           C1vsC2_Res = c( `1`="white", `0` = "black"),     
                           C1vsC3_B = c( `1`="white", `0` = "black"), 
                           C1vsC3_OS = c( `1`="white", `0` = "black"), 
                           C1vsC3_Pre = c( `1`="white", `0` = "black"), 
                           C1vsC3_Res = c( `1`="white", `0` = "black"), 
                           C2vsC2_Base_OS = c( `1`="white", `0` = "black"), 
                           C2vsC2_Base_pre = c( `1`="white", `0` = "black"), 
                           C2vsC2_Base_Res = c( `1`="white", `0` = "black"), 
                           C2vsC3_B = c( `1`="white", `0` = "black"), 
                           C3vsC3_Base_OS = c( `1`="white", `0` = "black"), 
                           C3vsC3_Base_pre = c( `1`="white", `0` = "black"), 
                           C3vsC3_Base_Res = c( `1`="white", `0` = "black"), 
                           C3vsC3_Pre_OS = c( `1`="white", `0` = "black"),
                           BaselinevsC1_Pre = c( `1`="white", `0` = "black"),
                           BaselinevsC1_Res = c( `1`="white", `0` = "black"),
                           BaselinevsC2_Pre = c( `1`="white", `0` = "black"),
                           BaselinevsC2_Res = c( `1`="white", `0` = "black"),
                           BaselinevsC3_Pre = c( `1`="white", `0` = "black"),
                           BaselinevsC3_Res = c( `1`="white", `0` = "black") )
###
quantile_breaks <- function(xs, n = 20) { breaks <- quantile(xs, probs = seq(0, 1, length.out = n)) ; breaks[!duplicated(breaks)] }
mat_breaks1 <- quantile_breaks( as.matrix(cimac_e4412_olink) , n = 20)
mat_breaks <- seq(min(cimac_e4412_olink), max(cimac_e4412_olink), length.out = 20)
mat_breaks2 <- seq(min(cimac_e4412_olink), max(cimac_e4412_olink), length.out = 10)
mat_breaks3 <- seq(min(cimac_e4412_olink), max(cimac_e4412_olink), length.out = 15)
###
cimac_e4412_olink_clustergram_meta$sorter <- cimac_e4412_olink_clustergram_meta$Collection_Event
cimac_e4412_olink_clustergram_meta$sorter[cimac_e4412_olink_clustergram_meta$sorter %in% "Baseline"] <- "a"
cimac_e4412_olink_clustergram_meta$sorter[cimac_e4412_olink_clustergram_meta$sorter %in% "Pre_D1_C2"] <- "b"
cimac_e4412_olink_clustergram_meta$sorter[cimac_e4412_olink_clustergram_meta$sorter %in% "Restaging"] <- "c"
cimac_e4412_olink_clustergram_meta$sorter[cimac_e4412_olink_clustergram_meta$sorter %in% "Off_study"] <- "d"
###
ix <- order(paste(cimac_e4412_olink_clustergram_meta$Cohort,cimac_e4412_olink_clustergram_meta$sorter))
###
cumsum(table(as.character(paste(annotation_col$Collection,annotation_col$Cohort))[ix]) [ match(unique(as.character(paste(annotation_col$Collection,annotation_col$Cohort))[ix]) ,  names(table(as.character(paste(annotation_col$Collection,annotation_col$Cohort))[ix])) ) ])
iy <- as.numeric(cumsum(table(as.character(paste(annotation_col$Collection,annotation_col$Cohort))[ix]) [ match(unique(as.character(paste(annotation_col$Collection,annotation_col$Cohort))[ix]) ,  names(table(as.character(paste(annotation_col$Collection,annotation_col$Cohort))[ix])) ) ]))
###
### NOTE: change mat_breaks, between 20 and 10 for highlighting colors
###
pdf(file="figures/heatmap_sig_conditions_NPX_15.pdf",width = 14, height = 6)
pheatmap( cimac_e4412_olink[ii,ix] , 
          #color = colorRampPalette(c("darkpurple","yellow"))(20),         
          color = viridis(15),
          cluster_cols = F, 
          angle_col = 90,
          cluster_rows = T, 
          show_colnames = FALSE,
          border_color = "black",
          treeheight_row = 0 ,
          show_rownames = TRUE, 
          #gaps_col = iy,
          gaps_col = c(67,115),
          annotation_col = annotation_col,
          annotation_row = annotation_row, 
          annotation_colors = annotation_colors, 
          clustering_distance_rows ="euclidean", 
          breaks = mat_breaks3,
          scale="none")
dev.off()
###
pdf(file="figures/heatmap_sig_conditions_zscore.pdf",width = 14, height = 6)
pheatmap( cimac_e4412_olink[ii,ix] , 
          #color = colorRampPalette(c("darkpurple","yellow"))(20),         
          color = viridis(20),
          cluster_cols = F, 
          angle_col = 90,
          cluster_rows = T, 
          show_colnames = FALSE,
          border_color = "black",
          treeheight_row = 0 ,
          show_rownames = TRUE, 
          #gaps_col = iy,
          gaps_col = c(67,115),
          annotation_col = annotation_col,
          annotation_row = annotation_row, 
          annotation_colors = annotation_colors, 
          clustering_distance_rows ="euclidean", 
          #breaks = mat_breaks,
          scale="row")
dev.off()
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
####
### saved at the end of the script
### write.csv(file="olink_e4412_results_table_eegk.csv",lm_olink_all[ lm_olink_all$adj.P.Val<0.05,])
####
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
### CYTOF
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
### orloj package dependencies: "flowCore", "FlowSOM", "edger","patchwork
### devtools::install_github("astrolabediagnostics/orloj")
library(orloj)
#############################################################################################################################################
### load data
e4412_cytof_experiment <- orloj::loadExperiment("/Users/gonzae34/Documents/projects_gnjatic/U24/U24/E4412_Diefenbach/uncompressed_files_astrolabel/")
### view summary
orloj::experimentSummary(e4412_cytof_experiment)
###
print(e4412_cytof_experiment$samples)
### single sample summary
# sample <- orloj::loadSample(e4412_cytof_experiment, sample_name = "sampleNAME")
# orloj::sampleSummary(sample)
# subset cell counts and channel intensity medians
# sampleCellSubsetCounts(sample)
e4412_sample_features <- orloj::getExperimentSampleFeatures(e4412_cytof_experiment)
###
cell_frequency_per_sample <- experimentCellSubsetCounts(e4412_cytof_experiment) %>% dplyr::arrange(CellSubset)
### subset channel intensity statistics:
channel_intesity_statistics <- experimentCellSubsetChannelStatistics(e4412_cytof_experiment) %>% dplyr::arrange(CellSubset)
### DA Analysis
e4412_daa <- orloj::differentialAbundanceAnalysis(e4412_cytof_experiment)
###
#############################################################################################################################################
### START wrangling
#############################################################################################################################################
table(cimac_e4412_metadata$Participant_ID %in% cell_frequency_per_sample$SampleId)
generic_var <- cell_frequency_per_sample[,c(1:3)]
generic_var <- unique(generic_var)
###
table(generic_var$Filename %in% cimac_e4412_metadata$File_Name)
table(cimac_e4412_metadata$File_Name %in% generic_var$Filename)
###
colnames(generic_var) <- c("astrolabe_ID","Name","File_Name")
###
cimac_e4412_metadata <- merge(cimac_e4412_metadata, generic_var, by="File_Name")
###
dim(cimac_e4412_metadata)
rm(generic_var)
###
generic_var <-  cimac_e4412_olink_clustergram_meta[, colnames(cimac_e4412_olink_clustergram_meta) %in% c("Arm","Cohort") ]
generic_var <- unique(generic_var)
###
colnames(cimac_e4412_metadata) [ colnames(cimac_e4412_metadata) %in% "Cohort" ] <- "Arm"
###
cimac_e4412_metadata <- merge(cimac_e4412_metadata, generic_var, by="Arm")
###
cimac_e4412_metadata$condition <- paste(cimac_e4412_metadata$Timepoint,cimac_e4412_metadata$Cohort, sep="_")
###
#############################################################################################################################################
### END wrangling
#############################################################################################################################################

#############################################################################################################################################
### DA re-analysis
#############################################################################################################################################
#############################################################################################################################################
# head(cell_frequency_per_sample)
### unmelt
e4412_frequency_matrix <- dcast(data = cell_frequency_per_sample,formula = SampleId~CellSubset, value.var = "N")
###
# head(e4412_frequency_matrix)
###
e4412_frequency_matrix <- e4412_frequency_matrix [ e4412_frequency_matrix$SampleId %in%  cimac_e4412_metadata$astrolabe_ID , ]
rownames(e4412_frequency_matrix) <- e4412_frequency_matrix[,1]
e4412_frequency_matrix <- e4412_frequency_matrix[,-1]
#############################################################################################################################################
###
generic_var <- t(e4412_frequency_matrix)### sort
generic_var <- generic_var[ , match(cimac_e4412_metadata$astrolabe_ID , colnames(generic_var)) ]
### check
identical( colnames(generic_var) ,cimac_e4412_metadata$astrolabe_ID ) # check TRUE
### convert to a log scale using DEG/voom
d.cpm.x <- DGEList(counts=generic_var) #Create a Abundance object
d.cpm.x <- calcNormFactors(d.cpm.x, method = "TMM") #method TMM
## d.cpm.x$samples$norm.factors
###
cimac_e4412_metadata$sorter1 <- as.character(cimac_e4412_metadata$Timepoint)
cimac_e4412_metadata$sorter1[cimac_e4412_metadata$sorter1 %in% "Baseline"] <- "A"
cimac_e4412_metadata$sorter1[cimac_e4412_metadata$sorter1 %in% "Pre_Day_1_Cycle_2"] <- "B"
cimac_e4412_metadata$sorter1[cimac_e4412_metadata$sorter1 %in% "Restaging"] <- "C"
cimac_e4412_metadata$sorter1[cimac_e4412_metadata$sorter1 %in% "End_of_treatment"] <- "D"
cimac_e4412_metadata$sorter1[cimac_e4412_metadata$sorter1 %in% "Off_study"] <- "E"
###
cimac_e4412_metadata$sorter2 <- as.character(cimac_e4412_metadata$Cohort)
cimac_e4412_metadata$sorter2[cimac_e4412_metadata$sorter2 %in% "3_BV+nivo+ipi"] <- "3"
cimac_e4412_metadata$sorter2[cimac_e4412_metadata$sorter2 %in% "2_BV+nivo"] <- "2"
cimac_e4412_metadata$sorter2[cimac_e4412_metadata$sorter2 %in% "1_BV+ipi"] <- "1"
###
cimac_e4412_metadata$sorter <- paste(cimac_e4412_metadata$sorter2,cimac_e4412_metadata$sorter1,sep="_")
###
table(cimac_e4412_metadata$sorter)
###
### Design and Voom
###
design <- model.matrix(~ 0 + sorter , data=cimac_e4412_metadata)
###
colnames(design) <- gsub("sorter","s",colnames(design))
###
d.cpm.x = estimateCommonDisp(d.cpm.x, verbose=TRUE)
v_abundance <- voom(d.cpm.x, design, plot=TRUE)
###
### Linear Model Applied
###
vfit <- lmFit(v_abundance, design)
###
efit <- eBayes(vfit)
###
summary(decideTests(efit)) 
###
contr.matrix <- makeContrasts( "1A_vs_1B" = s1_A - s1_B,
                               "1A_vs_1C" = s1_A - s1_C,
                               "1A_vs_1D" = s1_A - s1_D,
                               "1B_vs_1C" = s1_B - s1_C,
                               "1B_vs_1D" = s1_B - s1_D,
                               "1D_vs_1C" = s1_D - s1_C,
                               
                               "2A_vs_2B" = s2_A - s2_B,
                               "2A_vs_2C" = s2_A - s2_C,
                               "2A_vs_2D" = s2_A - s2_D,
                               "2B_vs_2C" = s2_B - s2_C,
                               "2B_vs_2D" = s2_B - s2_D,
                               "2D_vs_2C" = s2_D - s2_C,
                               
                               "3A_vs_3B" = s3_A - s3_B,
                               "3A_vs_3C" = s3_A - s3_C,
                               "3A_vs_3D" = s3_A - s3_D,
                               "3B_vs_3C" = s3_B - s3_C,
                               "3B_vs_3D" = s3_B - s3_D,
                               "3D_vs_3C" = s3_D - s3_C,
                               
                               "1A_vs_2A" = s1_A - s2_A,
                               "2A_vs_3A" = s2_A - s3_A,
                               "1A_vs_3A" = s1_A - s3_A,
                               
                               "1B_vs_2B" = s1_B - s2_B,
                               "2B_vs_3B" = s2_B - s3_B,
                               "1B_vs_3B" = s1_B - s3_B,
                               
                               "1C_vs_2C" = s1_C - s2_C,
                               "2C_vs_3C" = s2_C - s3_C,
                               "1C_vs_3C" = s1_C - s3_C,
                               
                               "1D_vs_2D" = s1_D - s2_D,
                               "2D_vs_3D" = s2_D - s3_D,
                               "1D_vs_3D" = s1_D - s3_D,
                               
                               "Baseline_vs_1B" = (s1_A+s2_A+s3_A)/3 - s1_B,
                               "Baseline_vs_1C" = (s1_A+s2_A+s3_A)/3 - s1_C,
                               "Baseline_vs_1D" = (s1_A+s2_A+s3_A)/3 - s1_D,
                               "Baseline_vs_2B" = (s1_A+s2_A+s3_A)/3 - s2_B,
                               "Baseline_vs_2C" = (s1_A+s2_A+s3_A)/3 - s2_C,
                               "Baseline_vs_2D" = (s1_A+s2_A+s3_A)/3 - s2_D,
                               "Baseline_vs_3B" = (s1_A+s2_A+s3_A)/3 - s3_B,
                               "Baseline_vs_3C" = (s1_A+s2_A+s3_A)/3 - s3_C,
                               "Baseline_vs_3D" = (s1_A+s2_A+s3_A)/3 - s3_D,
                               
                               levels = colnames(design) )
###
### my_disp <- estimateDisp(generic_var, design, trend="none")
block2 <- as.numeric(as.factor( cimac_e4412_metadata$Participant_ID ))
### function estimates the correlation between repeated observations
dupcor2 <- duplicateCorrelation(v_abundance, design, block=block2)
### model accounting for replicates
vfit <- lmFit(v_abundance, design, block=block2, correlation=dupcor2$consensus)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix) ; efit <- eBayes(vfit)
summary(decideTests(efit))
#############################################################################################################################################
lm_frequencies <- list()
for ( i in 1:length( colnames(summary(decideTests(efit))) ) ) {
  lm_frequencies[[i]] <- topTable(efit,coef= colnames(summary(decideTests(efit)))[i] ,n=Inf)  
  lm_frequencies[[i]]$Comparison <- colnames(summary(decideTests(efit)))[i]
  lm_frequencies[[i]]$Celltype <- rownames(lm_frequencies[[i]])
}
###
lm_frequencies_all <- do.call(rbind,lm_frequencies)
#############################################################################################################################################
### STRONG Correction for baseline/cohort
modcombat = model.matrix(~1, data=cimac_e4412_metadata)
sva_x_matrix = ComBat(dat=v_abundance$E, batch=cimac_e4412_metadata$Cohort, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
### Contrast matrix START
contr.matrix <- makeContrasts( "1A_vs_1B" = s1_A - s1_B,
                               "1A_vs_1C" = s1_A - s1_C,
                               "1A_vs_1D" = s1_A - s1_D,
                               "1B_vs_1C" = s1_B - s1_C,
                               "1B_vs_1D" = s1_B - s1_D,
                               "1D_vs_1C" = s1_D - s1_C,
                               
                               "2A_vs_2B" = s2_A - s2_B,
                               "2A_vs_2C" = s2_A - s2_C,
                               "2A_vs_2D" = s2_A - s2_D,
                               "2B_vs_2C" = s2_B - s2_C,
                               "2B_vs_2D" = s2_B - s2_D,
                               "2D_vs_2C" = s2_D - s2_C,
                               
                               "3A_vs_3B" = s3_A - s3_B,
                               "3A_vs_3C" = s3_A - s3_C,
                               "3A_vs_3D" = s3_A - s3_D,
                               "3B_vs_3C" = s3_B - s3_C,
                               "3B_vs_3D" = s3_B - s3_D,
                               "3D_vs_3C" = s3_D - s3_C,
                               
                               "1A_vs_2A" = s1_A - s2_A,
                               "2A_vs_3A" = s2_A - s3_A,
                               "1A_vs_3A" = s1_A - s3_A,
                               
                               "1B_vs_2B" = s1_B - s2_B,
                               "2B_vs_3B" = s2_B - s3_B,
                               "1B_vs_3B" = s1_B - s3_B,
                               
                               "1C_vs_2C" = s1_C - s2_C,
                               "2C_vs_3C" = s2_C - s3_C,
                               "1C_vs_3C" = s1_C - s3_C,
                               
                               "1D_vs_2D" = s1_D - s2_D,
                               "2D_vs_3D" = s2_D - s3_D,
                               "1D_vs_3D" = s1_D - s3_D,
                               
                               "Baseline_vs_1B" = (s1_A+s2_A+s3_A)/3 - s1_B,
                               "Baseline_vs_1C" = (s1_A+s2_A+s3_A)/3 - s1_C,
                               "Baseline_vs_1D" = (s1_A+s2_A+s3_A)/3 - s1_D,
                               "Baseline_vs_2B" = (s1_A+s2_A+s3_A)/3 - s2_B,
                               "Baseline_vs_2C" = (s1_A+s2_A+s3_A)/3 - s2_C,
                               "Baseline_vs_2D" = (s1_A+s2_A+s3_A)/3 - s2_D,
                               "Baseline_vs_3B" = (s1_A+s2_A+s3_A)/3 - s3_B,
                               "Baseline_vs_3C" = (s1_A+s2_A+s3_A)/3 - s3_C,
                               "Baseline_vs_3D" = (s1_A+s2_A+s3_A)/3 - s3_D,
                               
                               levels = colnames(design) )
### re-running model 
vfit <- lmFit(sva_x_matrix, design) #,block=block2, correlation=dupcor2$consensus, method = "robust")
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE) ;summary(decideTests(efit))
#############################################################################################################################################
lm_frequencies_cohort <- list()
for ( i in 1:length( colnames(summary(decideTests(efit))) ) ) {
  lm_frequencies_cohort[[i]] <- topTable(efit,coef= colnames(summary(decideTests(efit)))[i] ,n=Inf)  
  lm_frequencies_cohort[[i]]$Comparison <- colnames(summary(decideTests(efit)))[i]
  lm_frequencies_cohort[[i]]$Celltype <- rownames(lm_frequencies_cohort[[i]])
}
###
lm_frequencies_cohort_all <- do.call(rbind,lm_frequencies_cohort)
#############################################################################################################################################
###
rm(lm_frequencies,lm_frequencies_cohort)
dim(lm_frequencies_cohort_all)
dim(lm_frequencies_all)
###
table(lm_frequencies_all$adj.P.Val<0.05) ###
lm_frequencies_all[ lm_frequencies_all$adj.P.Val<0.05 , ]
###
head(lm_frequencies_all)
###
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
### Diff Expr. Summary.
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
stacked_bar_plot <- lm_frequencies_all 
stacked_bar_plot$logFC_dir <- stacked_bar_plot$logFC > 0 
stacked_bar_plot$logFC_dir[stacked_bar_plot$logFC > 0] <- "Up"
stacked_bar_plot$logFC_dir[stacked_bar_plot$logFC < 0] <- "Down"
stacked_bar_plot$sig <- stacked_bar_plot$adj.P.Val<0.05
#############################################################################################################################################
stacked_bar_plot <- data.frame(table(stacked_bar_plot$Comparison, stacked_bar_plot$logFC_dir, stacked_bar_plot$sig))
colnames(stacked_bar_plot) <- c('Comparison', 'logFC_dir','sig', 'count')
#############################################################################################################################################
regional_summary <- stacked_bar_plot %>%  filter(sig=='TRUE') %>% dplyr::select(-sig)
regional_summary$count[regional_summary$logFC_dir=='Down'] <- 0 - regional_summary$count[regional_summary$logFC_dir=='Down']
#############################################################################################################################################
regional_summary$Comparison <- as.character(regional_summary$Comparison) 
#############################################################################################################################################
pdf(file="figures/linear_model_all_data_ee4412_cytof_frequencies.pdf",width = 4, height = 6)
ggplot(regional_summary) + 
  aes(x=count, y=Comparison, fill=logFC_dir, label=abs(count)) + 
  geom_bar(stat='identity', width=0.5) +
  geom_text(data=regional_summary[regional_summary$count != 0,], show.legend=FALSE, color='black') +
  theme_classic() + theme(axis.text.y = element_text(face='bold'),
                          axis.text.x = element_text(face='bold'),
                          axis.title = element_text(face = "bold"),
                          legend.title= element_text(face='bold'),
                          strip.text=element_text(face='bold'),
                          plot.title=element_text(face='bold')) + 
  scale_fill_manual(values=c('firebrick', 'steelblue'), name='LogFC Direction') +
  labs(x='no. DA CellTypes (FDR < 0.05)', y='Contrast') + theme(legend.position="bottom")
dev.off()
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
###
### DAA HEATMAP
###
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
table(lm_frequencies_all[ lm_frequencies_all$adj.P.Val<0.05 , ]$Celltype)
table(e4412_daa$Cohort$FDR<0.05)
### ????????????????????????????
generic_var <- table(lm_frequencies_all[ lm_frequencies_all$adj.P.Val<0.05 , ]$Celltype,lm_frequencies_all[ lm_frequencies_all$adj.P.Val<0.05 , ]$Comparison)
class(generic_var) <- "matrix"
generic_var <- as.data.frame(generic_var)
###
#############################################################################################################################################
### plot
#############################################################################################################################################
###
annotation_row <- data.frame( Cohort = cimac_e4412_metadata$Cohort, Timepoint = cimac_e4412_metadata$Timepoint )
rownames(annotation_row) <- cimac_e4412_metadata$astrolabe_ID
###
e4412_frequency_matrix <- e4412_frequency_matrix[ match(rownames(annotation_row), rownames(e4412_frequency_matrix) ), ]
###
identical(rownames(e4412_frequency_matrix),rownames(annotation_row))
identical(rownames(e4412_frequency_matrix),cimac_e4412_metadata$astrolabe_ID)
###
annotation_col <- data.frame( AveFreq = rowMedians( t(e4412_frequency_matrix) )*100 )
rownames(annotation_col) <- rownames(t(e4412_frequency_matrix))
###
a <- annotation_col ; b <- generic_var
a$id <- rownames(a) ; b$id <- rownames(b)
annotation_col <- merge(a,b,by="id",all=TRUE) ; rm(a,b)
###
rownames(annotation_col) <- annotation_col$id
annotation_col[is.na(annotation_col)] <- 0
annotation_col <- annotation_col[,-1]
###
annotation_colors <- list(
  Cohort = c( `1_BV+ipi`= "black",`2_BV+nivo`= "grey80", `3_BV+nivo+ipi`= "grey60") ,
  Timepoint = c( Baseline= "black",End_of_treatment="white",Other="red",Off_study="blue",Pre_Day_1_Cycle_2="orange",Restaging="grey50"),
  `1A_vs_1B` = c(`1`="white",`0`="black"), `1A_vs_1C` = c(`1`="white",`0`="black"), `1A_vs_1D` = c(`1`="white",`0`="black"), 
  `1A_vs_2A` = c(`1`="white",`0`="black"), `1A_vs_3A` = c(`1`="white",`0`="black"), `1C_vs_2C` = c(`1`="white",`0`="black"), 
  `1D_vs_2D` = c(`1`="white",`0`="black"), `1D_vs_3D` = c(`1`="white",`0`="black"), `2A_vs_2B` = c(`1`="white",`0`="black"),
  `2B_vs_2D` = c(`1`="white",`0`="black"), `2D_vs_3D` = c(`1`="white",`0`="black"), `3A_vs_3B` = c(`1`="white",`0`="black"), 
  `3B_vs_3D` = c(`1`="white",`0`="black"), `2A_vs_2D` = c(`1`="white",`0`="black"), `Baseline_vs_1B` = c(`1`="white",`0`="black"),
  `Baseline_vs_1C` = c(`1`="white",`0`="black"), `Baseline_vs_1D` = c(`1`="white",`0`="black"), `Baseline_vs_3B` = c(`1`="white",`0`="black")
) 
### Ave Freq
annotation_col$AveFreq <- round(annotation_col$AveFreq/sum(annotation_col$AveFreq)*100,2)
###
cumsum(table(annotation_row$Cohort[order(cimac_e4412_metadata$sorter)]))
###
for_heatmap <- v_abundance$E[,order(cimac_e4412_metadata$sorter)]
table(annotation_row$Timepoint %in% c("Other","Off_study"))
### boring samples to remove
ii <- ! colnames(for_heatmap) %in% rownames(annotation_row)[annotation_row$Timepoint %in% c("Other","Off_study")]
for_heatmap <- for_heatmap[,ii]
###
identical(rownames(annotation_col),rownames(for_heatmap))
ix <- ! rownames(for_heatmap) %in% c("Root_unassigned")
ii <- order(annotation_col$AveFreq[ix],decreasing = TRUE)
for_heatmap<-for_heatmap[ix,]
for_heatmap<-for_heatmap[ii,]
# pheatmap(e4412_frequency_matrix[ order(paste(annotation_row$Cohort,annotation_row$Timepoint)),], 
pdf(file="figures/cytof_cell_frequencies_per_sample.pdf",width = 12, height = 6)
pheatmap(for_heatmap,
         scale="none", 
         clustering_distance_rows ="euclidean", 
         show_rownames = TRUE, 
         show_colnames = FALSE, 
         #gaps_col = c(65,113),
         gaps_col = c(62,109),
         annotation_col = annotation_row, 
         annotation_row = annotation_col,
         annotation_colors = annotation_colors, 
         # color = colorRampPalette(c("blue","white","orangered"))(255), 
         color = viridis(255), 
         angle_col = 90,
         border_color = "grey70",
         cluster_cols = F, 
         cluster_rows = F)
dev.off()
#############################################################################################################################################

#############################################################################################################################################

### channel_intesity_statistics
### cell_frequency_per_sample
### e4412_daa
### e4412_daa$Cohort[e4412_daa$Cohort$FDR<0.05,]

#############################################################################################################################################
my_comparisons <- list( c("Baseline","Pre_Day_1_Cycle_2"),
                        c("Baseline","Restaging"),
                        c("Pre_Day_1_Cycle_2","Restaging"),
                        c("Baseline","End_of_treatment"),
                        c("Pre_Day_1_Cycle_2","End_of_treatment"),
                        c("Restaging","End_of_treatment") )
#############################################################################################################################################

for ( i in 1:length(levels(as.factor(cell_frequency_per_sample$CellSubset))) ) {
  # levels(as.factor(cell_frequency_per_sample$CellSubset))
  ii <- cell_frequency_per_sample$CellSubset %in% levels(as.factor(cell_frequency_per_sample$CellSubset))[i] & cell_frequency_per_sample$SampleId %in% cimac_e4412_metadata$astrolabe_ID
  generic_var <- cell_frequency_per_sample[ii,]       
  
  colnames(generic_var)[colnames(generic_var) %in% "SampleId"] <- "astrolabe_ID"
  generic_var <- merge(generic_var,cimac_e4412_metadata,by="astrolabe_ID")
  
  my_cytof_plots[[i]] <- ggviolin(data=generic_var, x="time_points", y="Freq", color="Cohort",fill="Cohort") + 
    theme_linedraw() + facet_wrap(~Cohort,ncol=3) + 
    geom_boxplot(width=0.33) + 
    geom_quasirandom(varwidth = TRUE,alpha=0.7) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
    theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
    ggtitle( paste(levels(as.factor(cell_frequency_per_sample$CellSubset))[i]) ) + labs(x = "", y = "Frequency") + theme(legend.position="none")
}

#############################################################################################################################################

pdf(file="figures/violinplots_all_data_ee4412_cytof.pdf",width = 8, height = 6)
# my_cytof_plots[[8]] ### issue outputing the list as pdf, hence loop bypasses the issue
for(i in 1:length(my_cytof_plots)){print((my_cytof_plots[[i]]))}
dev.off()

#############################################################################################################################################

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
###
### dim red
###
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

mydata <- prcomp(e4412_frequency_matrix , scale=TRUE, center = TRUE ) 
###
scree_pca <- fviz_eig(mydata, addlabels = T, barfill = "lightsteelblue3", barcolor = "lightsteelblue3") + theme_bw() + ggtitle("") #+ 
#theme(plot.title = element_text(hjust = 0.5), text = element_text(size=rel(4.8)) )
###
pca_mts <- as.data.frame(mydata$x)
pca_mts$astrolabe_ID <- rownames(pca_mts)
###
library(umap)
umap_clustering <- umap((e4412_frequency_matrix[,]))
pca_mts$UMAP1 <- umap_clustering$layout[,1]
pca_mts$UMAP2 <- umap_clustering$layout[,2]
###
library(Rtsne)
tsne_clustering <- Rtsne((e4412_frequency_matrix), dims = 2, initial_dims = 30,
                         perplexity = 15, theta = 1, check_duplicates = TRUE,
                         pca = TRUE, partial_pca = FALSE, max_iter = 5000,
                         normalize = TRUE, momentum = 0.5, final_momentum = 0.8, eta = 200, exaggeration_factor = 12, num_threads = 1)
###      
pca_mts$tsne1 <- tsne_clustering$Y[,1]
pca_mts$tsne2 <- tsne_clustering$Y[,2]
###
pca_mts <- merge(pca_mts, cimac_e4412_metadata, by="astrolabe_ID")
###
pdf(file="figures/reductions_all_data_ee4412_cytoff.pdf",width = 14, height = 10)
ggarrange( ggplot(pca_mts, aes(PC1, PC2, color = time_points)) + theme_bw() + geom_point(alpha=0.5) ,
           ggplot(pca_mts, aes(UMAP1, UMAP2, color = time_points)) + theme_bw() + geom_point(alpha=0.5) , 
           scree_pca, 
           ggplot(pca_mts, aes(tsne1, tsne2, color = time_points)) + theme_bw() + geom_point(alpha=0.5) , 
           ncol=2,nrow=2 )
dev.off()
###

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
### 
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

### channel_intesity_statistics
### cell_frequency_per_sample
### e4412_daa
###
ii <- channel_intesity_statistics$CellSubset %in% channel_intesity_statistics$CellSubset[1] & channel_intesity_statistics$SampleId %in% cimac_e4412_metadata$astrolabe_ID
###
generic_var <- channel_intesity_statistics[ii,]
###
generic_var <- dcast(data = generic_var,formula = SampleId~ChannelName, value.var = "Median")
###
rownames(generic_var) <- generic_var[,1]
generic_var <- generic_var[,-1]
###
annotation_row <- data.frame( Cohort = cimac_e4412_metadata$Cohort, Timepoint = cimac_e4412_metadata$Timepoint )
rownames(annotation_row) <- cimac_e4412_metadata$astrolabe_ID
###
annotation_colors <- list(
  Cohort = c( `1_BV+ipi`= "black",`2_BV+nivo`= "grey80", `3_BV+nivo+ipi`= "grey60") ,
  Timepoint = c( Baseline= "black",End_of_treatment="white",Other="red",Off_study="blue",Pre_Day_1_Cycle_2="orange",Restaging="grey50")
) 
###
identical(rownames(generic_var),rownames(annotation_row))
###
generic_var [ is.na(generic_var) ] <- 0
###

#############################################################################################################################################
pdf(file="figures/cytof_marker_frequencies_per_sample.pdf",width = 8, height = 4)
### Markers with at least some expression
pheatmap(generic_var[,colMeans(generic_var) > 0.1],
         scale="none", 
         main=channel_intesity_statistics$CellSubset[1],
         clustering_distance_rows ="euclidean", 
         show_rownames = FALSE, 
         show_colnames = TRUE, 
         annotation_row = annotation_row, 
         annotation_colors = annotation_colors, 
         color = colorRampPalette(c("blue","white","orangered"))(255), 
         angle_col = 90,
         border_color = "grey70",
         cluster_cols = T, 
         cluster_rows = T)
dev.off()
#############################################################################################################################################

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

# levels(as.factor(channel_intesity_statistics$CellSubset))
channel_intesity_statistics$sorter <- paste(channel_intesity_statistics$SampleId,channel_intesity_statistics$CellSubset,sep="_")
###
generic_var <- dcast(data = channel_intesity_statistics,formula = sorter~ChannelName, value.var = "Median")
###
dim(generic_var)
###
rownames(generic_var) <- generic_var[,1]
generic_var <- generic_var[,-1]
###
generic_var[1:5,1:5]
###
a <- data.frame(  astrolabe_ID = tidyr::separate(data.frame( rownames(generic_var) ), 1, sep="_", c("a","b"))$a , 
                  CellSubset = tidyr::separate(data.frame( rownames(generic_var) ), 1, sep="_", c("a","b"))$b )
###
generic_var <- generic_var[a$astrolabe_ID %in% cimac_e4412_metadata$astrolabe_ID, ]
a <- a[a$astrolabe_ID %in% cimac_e4412_metadata$astrolabe_ID,]
###
generic_var <- generic_var[order(a$astrolabe_ID),]
a <- a[order(a$astrolabe_ID),]
###
annotation_row <- a ; rm(a)
annotation_row <- merge(annotation_row,cimac_e4412_metadata,by="astrolabe_ID")
###
rownames(annotation_row) <- paste(annotation_row$astrolabe_ID,annotation_row$CellSubset,sep="_")
rownames(annotation_row) <- rownames(generic_var)
###
annotation_row <- annotation_row[, colnames(annotation_row) %in% c("Cohort","Timepoint") ]

###
annotation_colors <- list(
  Cohort = c( `1_BV+ipi`= "black",`2_BV+nivo`= "grey80", `3_BV+nivo+ipi`= "grey60") ,
  Timepoint = c( Baseline= "black",End_of_treatment="white",Other="red",Off_study="blue",Pre_Day_1_Cycle_2="orange",Restaging="grey50")
) 
###
identical(rownames(generic_var),rownames(annotation_row))
###
generic_var [ is.na(generic_var) ] <- 0
###
#############################################################################################################################################
pdf(file="figures/cytof_markers_mean_frequencies_celltype_sample.pdf",width = 8, height = 4)
### Markers with at least some expression
pheatmap(generic_var,
         scale="none", 
         # main=channel_intesity_statistics$CellSubset[1],
         clustering_distance_rows ="euclidean", 
         show_rownames = FALSE, 
         show_colnames = TRUE, 
         annotation_row = annotation_row, 
         annotation_colors = annotation_colors, 
         #color = colorRampPalette(c("steelblue","white","firebrick"))(255), 
         color = viridis(255),
         angle_col = 90,
         border_color = "grey70",
         cluster_cols = T, 
         cluster_rows = T)
dev.off()
#############################################################################################################################################

#############################################################################################################################################
### test to identify if sample_batch had an effect on the samples
#############################################################################################################################################

### PCA first
#rownames(generic_var)
#annotation_row
###
a <- data.frame(  astrolabe_ID = tidyr::separate(data.frame( rownames(generic_var) ), 1, sep="_", c("a","b"))$a , 
                  CellSubset = tidyr::separate(data.frame( rownames(generic_var) ), 1, sep="_", c("a","b"))$b )
generic_var <- generic_var[a$astrolabe_ID %in% cimac_e4412_metadata$astrolabe_ID, ]
a <- a[a$astrolabe_ID %in% cimac_e4412_metadata$astrolabe_ID,]
generic_var <- generic_var[order(a$astrolabe_ID),]
a <- a[order(a$astrolabe_ID),]
annotation_row <- a ; rm(a)
annotation_row <- merge(annotation_row,cimac_e4412_metadata,by="astrolabe_ID")
rownames(annotation_row) <- paste(annotation_row$astrolabe_ID,annotation_row$CellSubset,sep="_")
rownames(annotation_row) <- rownames(generic_var)
annotation_row$sorter <- rownames(generic_var)
###

mydata <- prcomp( generic_var , scale=TRUE, center = TRUE ) 
###
scree_pca <- fviz_eig(mydata, addlabels = T, barfill = "lightsteelblue3", barcolor = "lightsteelblue3") + theme_bw() + ggtitle("") #+ 
#theme(plot.title = element_text(hjust = 0.5), text = element_text(size=rel(4.8)) )
###
pca_mts <- as.data.frame(mydata$x)
pca_mts$sorter <- rownames(pca_mts)
###
library(umap)
umap_clustering <- umap(generic_var)
pca_mts$UMAP1 <- umap_clustering$layout[,1]
pca_mts$UMAP2 <- umap_clustering$layout[,2]
###
library(Rtsne)
tsne_clustering <- Rtsne(generic_var, dims = 2, initial_dims = 30,
                         perplexity = 15, theta = 1, check_duplicates = FALSE,
                         pca = TRUE, partial_pca = FALSE, max_iter = 5000,
                         normalize = TRUE, momentum = 0.5, final_momentum = 0.8, eta = 200, exaggeration_factor = 12, num_threads = 1)
###      
pca_mts$tsne1 <- tsne_clustering$Y[,1]
pca_mts$tsne2 <- tsne_clustering$Y[,2]
###
head(as.character(annotation_row$sorter))
head(as.character(pca_mts$sorter))
table(as.character(annotation_row$sorter) %in% as.character(pca_mts$sorter))
pca_mts <- merge(pca_mts, annotation_row, by="sorter")
table(pca_mts$Sample_Batch)
table(pca_mts$Veritas_batch)
###

pdf(file="figures/reductions_all_data_cytof_cell_subtype_specificity.pdf",width = 24, height = 20)
ggarrange( ggplot(pca_mts, aes(PC1, PC2, color = CellSubset)) + theme_bw() + geom_point(alpha=0.5) ,
           ggplot(pca_mts, aes(UMAP1, UMAP2, color = CellSubset)) + theme_bw() + geom_point(alpha=0.5) , 
           scree_pca, 
           ggplot(pca_mts, aes(tsne1, tsne2, color = CellSubset)) + theme_bw() + geom_point(alpha=0.5) , 
           ncol=2,nrow=2 )
dev.off()

#############################################################################################################################################

pdf(file="figures/reduction_tsne_all_data_cytof_cell_subtype_specificity.pdf",width = 12, height = 6)
# ggplot(pca_mts, aes(UMAP1, UMAP2, color = CellSubset)) + theme_bw() + geom_point(alpha=0.5)
ggplot(pca_mts, aes(tsne1, tsne2, color = CellSubset)) + theme_bw() + geom_point(alpha=0.5)
#ii<-colnames(pca_mts) %in% c("UMAP1","UMAP2")
#length(table(pca_mts$CellSubset))
#autoplot(fanny(pca_mts[,ii], 29), frame = TRUE) ### Takes Forever
dev.off()

#############################################################################################################################################

### variance_partition second

### For all subsets
form <- ~ Participant_ID + Sample_Batch + CellSubset + Cohort + time_points
variance_exprs_matrix <- fitExtractVarPartModel(t(generic_var), form, pca_mts)
pdf(file="figures/reductions_all_data_cytof_variance_all_celltypes.pdf",width = 6, height = 4)
plotVarPart( sortCols( variance_exprs_matrix ) ) + ggtitle( 'Cytof Data') + theme(axis.text.x = element_text(angle = 60,hjust = 1))
dev.off()
### 
ii <- pca_mts$CellSubset %in% pca_mts$CellSubset[1]
form <- ~ Participant_ID + Sample_Batch + Cohort + time_points
variance_exprs_matrix <- fitExtractVarPartModel(t(generic_var[ii,]), form, pca_mts[ii,])
pdf(file="figures/reductions_all_data_cytof_variance_all_celltypes_B_cells_sample.pdf",width = 6, height = 4)
plotVarPart( sortCols( variance_exprs_matrix ) ) + ggtitle( 'Cytof Data') + theme(axis.text.x = element_text(angle = 60,hjust = 1))
dev.off()
#### 

#############################################################################################################################################


#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
### Marker Analysis per Cell Type
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################


#############################################################################################################################################
### Long loop to calculate the DEs per celltype [ START ]
#############################################################################################################################################

### add condition to metadata
cimac_e4412_metadata$Condition <- paste(cimac_e4412_metadata$Cohort,cimac_e4412_metadata$time_points,sep="_")
astrolabe_markers <- c("TIGIT", "PD_L1", "PD_1", "KLRG1", "ICOS", "CXCR5", "CXCR3", "CD73", "CD69", "CD39", "CD25", "CCR7", "CCR6", "CCR4", "4_1BB")
# levels(as.factor(channel_intesity_statistics$CellSubset))
###
library(umap)
### list of plots and object to store the results
lm_cytof_markers <- list()
lm_cytof_markers_cohort <- list()
my_plot <- list()
my_plot_cohort <- list()
### START OF LOOP
for (i in 1:length(levels(as.factor(channel_intesity_statistics$CellSubset)))) {
  ### subset 
  ii <- channel_intesity_statistics$CellSubset %in% levels(as.factor(channel_intesity_statistics$CellSubset))[i]
  ### cast a matrix from long table to rectangular matrix
  # generic_var <- dcast(data = channel_intesity_statistics[ii,],formula = SampleId~ChannelName, value.var = "Median") #
  generic_var <- dcast(data = channel_intesity_statistics[ii,],formula = SampleId~ChannelName, value.var = "Q0_95")
  ### prepare and filter matrix
  rownames(generic_var) <- generic_var[,1] ; generic_var <- generic_var[,-1]
  ### filter 50 control samples from astrolabe
  ii <- rownames(generic_var) %in% cimac_e4412_metadata$astrolabe_ID
  generic_var <- generic_var[ii,]
  ### order the matrix to match the metadata (easier to calculate the stats)
  generic_var <- generic_var[ match(cimac_e4412_metadata$astrolabe_ID,rownames(generic_var)) ,]
  ### verify they are correct // verified manually, could be inplemented a checkpoint if conditional.
  identical(rownames(generic_var),cimac_e4412_metadata$astrolabe_ID)
  ### Filter samples that have 1 sample per condition.
  kk <- c(grep("Off_study",cimac_e4412_metadata$Condition),grep("Other",cimac_e4412_metadata$Condition))
  generic_var <- generic_var[-kk,]
  ### subset metatada
  cimac_e4412_metadata_v2 <- cimac_e4412_metadata[ -kk , ]
  ### verify they are correct // verified manually, could be inplemented a checkpoint if conditional.
  identical(rownames(generic_var),cimac_e4412_metadata_v2$astrolabe_ID)
  ### design matrix fot the statistics
  design <- model.matrix( ~ 0 + Condition , data = cimac_e4412_metadata_v2 )
  ### Manually Relabel the files on design
  # colnames(design)
  # a = bseline, b = preday1cycle2, c=restaging, d= end of treatment
  colnames(design) <- c("C1_A","C1_D","C1_B","C1_C",
                        "C2_A","C2_D","C2_B","C2_C",
                        "C3_A","C3_D","C3_B","C3_C")
  ### model patient specificity into the mixed linear model
  block2 <- as.numeric(as.factor( cimac_e4412_metadata_v2$Participant_ID ))
  ### function estimates the correlation between repeated observations
  dupcor2 <- duplicateCorrelation(t(generic_var), design, block=block2)
  ### model accounting for replicates
  vfit <- lmFit(t(generic_var), design, block=block2, correlation=dupcor2$consensus)
  ### stats & summary
  efit <- eBayes(vfit) ; summary(decideTests(efit))
  ### Contrast matrix START
  contr.matrix <- makeContrasts( "C1vsC1_Base_pre" = C1_A - C1_B,
                                 "C2vsC2_Base_pre" = C2_A - C2_B,
                                 "C3vsC3_Base_pre" = C3_A - C3_B,
                                 
                                 "C1vsC1_Base_res" = C1_A - C1_C,
                                 "C2vsC2_Base_res" = C2_A - C2_C,
                                 "C3vsC3_Base_res" = C3_A - C3_C,
                                 
                                 "C1vsC1_Base_end" = C1_A - C1_D,
                                 "C2vsC2_Base_end" = C2_A - C2_D,
                                 "C3vsC3_Base_end" = C3_A - C3_D,
                                 
                                 "C1vsC1_Pre_res" = C1_B - C1_C,
                                 "C2vsC2_Pre_res" = C2_B - C2_C,
                                 "C3vsC3_Pre_res" = C3_B - C3_C,
                                 
                                 "C1vsC1_pre_end" = C1_B - C1_D,
                                 "C2vsC2_pre_end" = C2_B - C2_D,
                                 "C3vsC3_pre_end" = C3_B - C3_D,
                                 
                                 "C1vsC1_Res_end" = C1_C - C1_D,
                                 "C2vsC2_Res_end" = C2_C - C2_D,
                                 "C3vsC3_Res_end" = C3_C - C3_D,
                                 
                                 levels = colnames(design) )
  ### Contrast Matrix END
  
  ### Linear model and contrast fit
  vfit <- lmFit(t(generic_var), design, block=block2, correlation=dupcor2$consensus, method = "robust")
  vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
  efit <- eBayes(vfit, robust = TRUE) ;summary(decideTests(efit))
  ### nested loop for storing the efit results per comparison "per celltype"
  temp_variable <- list()
  for ( iz in 1:length( colnames(summary(decideTests(efit))) ) ) {
    temp_variable[[iz]] <- topTable(efit,coef= colnames(summary(decideTests(efit)))[iz] ,n=Inf)  
    temp_variable[[iz]]$Comparison <- colnames(summary(decideTests(efit)))[iz]
    temp_variable[[iz]]$Marker <- rownames(temp_variable[[iz]])
  } ; lm_cytof_markers[[i]] <- do.call(rbind,temp_variable)
  
  ### STRONG Correction for baseline/cohort
  modcombat = model.matrix(~1, data=cimac_e4412_metadata_v2)
  sva_x_matrix = ComBat(dat=t(generic_var), batch=cimac_e4412_metadata_v2$Cohort, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
  ### Contrast matrix START
  contr.matrix <- makeContrasts( "C1vsC2_Base" = C1_A - C2_A,
                                 "C1vsC3_Base" = C1_A - C3_A,
                                 "C2vsC3_Base" = C2_A - C3_A,
                                 
                                 "C1vsC2_Pre" = C1_B - C2_B,
                                 "C1vsC3_Pre" = C1_B - C3_B,
                                 "C2vsC3_Pre" = C2_B - C3_B,
                                 
                                 "C1vsC2_Res" = C1_C - C2_C,
                                 "C1vsC3_Res" = C1_C - C3_C,
                                 "C2vsC3_Res" = C2_C - C3_C,
                                 
                                 "C1vsC2_End" = C1_D - C2_D,
                                 "C1vsC3_End" = C1_D - C3_D,
                                 "C2vsC3_End" = C2_D - C3_D,
                                 levels = colnames(design) )
  ### re-running model 
  vfit <- lmFit(sva_x_matrix, design, block=block2, correlation=dupcor2$consensus, method = "robust")
  vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
  efit <- eBayes(vfit, robust = TRUE) ;summary(decideTests(efit))
  ### nested loop for storing the efit results per comparison "per celltype"
  temp_variable <- list()
  for ( iz in 1:length( colnames(summary(decideTests(efit))) ) ) {
    temp_variable[[iz]] <- topTable(efit,coef= colnames(summary(decideTests(efit)))[iz] ,n=Inf)  
    temp_variable[[iz]]$Comparison <- colnames(summary(decideTests(efit)))[iz]
    temp_variable[[iz]]$Marker <- rownames(temp_variable[[iz]])
  } ; lm_cytof_markers_cohort[[i]] <- do.call(rbind,temp_variable)
  
  
  ### section for making the heatmap
  annotation_row <- cimac_e4412_metadata[, colnames(cimac_e4412_metadata) %in% c("Cohort","Timepoint") ]
  rownames(annotation_row) <- cimac_e4412_metadata$astrolabe_ID
  ###
  # identical(rownames(generic_var),rownames(annotation_row)) ### this is false; not necessary, heatmap is labeled so it autosorts.
  ###
  generic_var [ is.na(generic_var) ] <- 0
  ###
  ### ix <- colnames(generic_var) %in% astrolabe_markers
  #############################################################################################################################################
  print(levels(as.factor(channel_intesity_statistics$CellSubset))[i])
  #############################################################################################################################################
  if(length(which(lm_cytof_markers[[i]]$adj.P.Val<0.05)) == 0) { print(paste("ERROR",levels(as.factor(channel_intesity_statistics$CellSubset))[i])) } else{
    ### construction of matrices for heatmap
    temp_variable <- table( lm_cytof_markers[[i]][ lm_cytof_markers[[i]]$adj.P.Val<0.05 , ]$Marker, lm_cytof_markers[[i]][ lm_cytof_markers[[i]]$adj.P.Val<0.05 , ]$Comparison )
    class(temp_variable) <- "matrix"
    temp_variable <- as.data.frame(temp_variable)
    ###
    annotation_col <- data.frame( Markers = colnames(generic_var) )
    rownames(annotation_col) <- annotation_col$Markers
    ### 
    a <- annotation_col ; b <- temp_variable
    a$id <- rownames(a) ; b$id <- rownames(b)
    annotation_col <- merge(a,b,by="id",all=TRUE) ; rm(a,b)
    ### annotation column construction
    rownames(annotation_col) <- annotation_col$id
    annotation_col[is.na(annotation_col)] <- 0
    annotation_col <- annotation_col[,-1]
    annotation_col <- annotation_col[,-1]
    ### colors for heatmap
    annotation_colors <- list(
      Cohort = c( `1_BV+ipi`= "black",`2_BV+nivo`= "grey80", `3_BV+nivo+ipi`= "grey60") ,
      Timepoint = c( Baseline= "black",End_of_treatment="white",Other="red",Off_study="blue",Pre_Day_1_Cycle_2="orange",Restaging="grey50"),
      `C1vsC1_Base_end` = c(`1`="white",`0`="black"), `C1vsC1_Base_pre` = c(`1`="white",`0`="black"), `C1vsC1_Base_res` = c(`1`="white",`0`="black"), 
      `C1vsC1_pre_end` = c(`1`="white",`0`="black"), `C1vsC1_Pre_res` = c(`1`="white",`0`="black"), `C1vsC1_Res_end` = c(`1`="white",`0`="black"), 
      `C1vsC2_A` = c(`1`="white",`0`="black"), `C1vsC2_B` = c(`1`="white",`0`="black"), `C1vsC2_C` = c(`1`="white",`0`="black"),
      `C1vsC2_D` = c(`1`="white",`0`="black"), 
      `C1vsC3_A` = c(`1`="white",`0`="black"), `C1vsC3_B` = c(`1`="white",`0`="black"), `C1vsC3_C` = c(`1`="white",`0`="black"),
      `C1vsC3_D` = c(`1`="white",`0`="black"), 
      `C2vsC2_Base_end` = c(`1`="white",`0`="black"), `C2vsC2_Base_pre` = c(`1`="white",`0`="black"), `C2vsC2_Base_res` = c(`1`="white",`0`="black"),
      `C2vsC2_pre_end` = c(`1`="white",`0`="black"), `C2vsC2_Pre_res` = c(`1`="white",`0`="black"), `C2vsC2_Res_end` = c(`1`="white",`0`="black"),
      `C2vsC3_A` = c(`1`="white",`0`="black"), `C2vsC3_B` = c(`1`="white",`0`="black"), `C2vsC3_C` = c(`1`="white",`0`="black"),`C2vsC3_D` = c(`1`="white",`0`="black"),
      `C3vsC3_Base_end` = c(`1`="white",`0`="black"), `C3vsC3_Base_pre` = c(`1`="white",`0`="black"), `C3vsC3_Base_res` = c(`1`="white",`0`="black"),`C3vsC3_pre_end` = c(`1`="white",`0`="black"),
      `C3vsC3_Pre_res` = c(`1`="white",`0`="black"),`C3vsC3_Res_end` = c(`1`="white",`0`="black")
    )
    ### for defining the number of clusters to sort in heatmap; using pca and kmeans.
    # set.seed(12345)
    ### PCA
    mydata <- prcomp( generic_var , scale=FALSE, center = TRUE ) 
    pca_mts <- as.data.frame(mydata$x) ; pca_mts$Sample_ID <- rownames(pca_mts)
    pca_mts$kmean <- as.character(kmeans(as.data.frame(mydata$x),3)$cluster)
    pca_mts$astrolabe_ID <- rownames(pca_mts)
    ### UMAP
    umap_clustering <- umap((generic_var[,]))
    pca_mts$UMAP1 <- umap_clustering$layout[,1]
    pca_mts$UMAP2 <- umap_clustering$layout[,2]
    pca_mts <- merge(pca_mts,cimac_e4412_metadata_v2,by="astrolabe_ID")
    # ggplot(pca_mts, aes(x=PC1,y=PC2,color=kmean)) + geom_point()
    # ggplot(pca_mts, aes(x=PC1,y=PC2,color=Cohort)) + geom_point()
    # ggplot(pca_mts, aes(x=UMAP1,y=UMAP2,color=Condition)) + geom_point()
    ### REMOVE OUTLIERS
    to_remove <- pca_mts$astrolabe_ID[ pca_mts$kmean %in% names(sort(table(pca_mts$kmean))[1])]
    ### annotation row annotation construction
    annotation_row <- annotation_row[ rownames(annotation_row) %in% rownames(generic_var) , ]
    xx <- order(paste(annotation_row$Cohort,annotation_row$Timepoint))
    ### check
    identical(rownames(generic_var),rownames(annotation_row))
    ### filter matrix
    generic_var <- generic_var[xx,]
    ### 
    if (length(to_remove) > 5){ 
      iy <- rownames(generic_var) %in% rownames(generic_var) 
    } else { iy <- ! rownames(generic_var) %in% to_remove }
    # summary(generic_var)
    # summary(generic_var [ , colnames(generic_var) %in% astrolabe_markers ])
    # pheatmap(t(generic_var[,ix]),
    ### check
    ### identical(rownames(generic_var)[iy],rownames(annotation_row)[iy])
    ### filter matrix
    generic_var <- generic_var[iy,]
    generic_var <- generic_var [ rowSums(generic_var) > 0 , colSums(generic_var) > 0 ]
    ### identical(rownames(generic_var),rownames(annotation_row))
    ###
    #print(paste(levels(as.factor(channel_intesity_statistics$CellSubset))[i], table(iy)))
    
    #############################################################################################################################################
    setwd("/Users/gonzae34/Documents/projects_gnjatic/U24/U24/E4412_Diefenbach/figures/celltype_figures")
    #############################################################################################################################################
    
    #############################################################################################################################################
    ### Results within each cohort
    stacked_bar_plot <- lm_cytof_markers[[i]] 
    stacked_bar_plot$logFC_dir <- stacked_bar_plot$logFC > 0 
    stacked_bar_plot$logFC_dir[stacked_bar_plot$logFC > 0] <- "Up"
    stacked_bar_plot$logFC_dir[stacked_bar_plot$logFC < 0] <- "Down"
    stacked_bar_plot <- stacked_bar_plot [ ! stacked_bar_plot$logFC == 0 , ]
    stacked_bar_plot$sig <- stacked_bar_plot$adj.P.Val<0.05
    ###
    stacked_bar_plot <- data.frame(table(stacked_bar_plot$Comparison, stacked_bar_plot$logFC_dir, stacked_bar_plot$sig))
    colnames(stacked_bar_plot) <- c('Comparison', 'logFC_dir','sig', 'count')
    ###
    regional_summary <- stacked_bar_plot %>%  filter(sig=='TRUE') %>% dplyr::select(-sig)
    regional_summary$count[regional_summary$logFC_dir=='Down'] <- 0 - regional_summary$count[regional_summary$logFC_dir=='Down']
    ###
    regional_summary$Comparison <- as.character(regional_summary$Comparison) 
    #############################################################################################################################################
    
    #############################################################################################################################################
    ### Results between cohorts
    stacked_bar_plot_cohort <- lm_cytof_markers_cohort[[i]] 
    stacked_bar_plot_cohort$logFC_dir <- stacked_bar_plot_cohort$logFC > 0 
    stacked_bar_plot_cohort$logFC_dir[stacked_bar_plot_cohort$logFC > 0] <- "Up"
    stacked_bar_plot_cohort$logFC_dir[stacked_bar_plot_cohort$logFC < 0] <- "Down"
    stacked_bar_plot_cohort <- stacked_bar_plot_cohort [ ! stacked_bar_plot_cohort$logFC == 0 , ]
    stacked_bar_plot_cohort$sig <- stacked_bar_plot_cohort$adj.P.Val<0.05
    ###
    stacked_bar_plot_cohort <- data.frame(table(stacked_bar_plot_cohort$Comparison, stacked_bar_plot_cohort$logFC_dir, stacked_bar_plot_cohort$sig))
    colnames(stacked_bar_plot_cohort) <- c('Comparison', 'logFC_dir','sig', 'count')
    ###
    regional_summary_cohort <- stacked_bar_plot_cohort %>%  filter(sig=='TRUE') %>% dplyr::select(-sig)
    regional_summary_cohort$count[regional_summary_cohort$logFC_dir=='Down'] <- 0 - regional_summary_cohort$count[regional_summary_cohort$logFC_dir=='Down']
    ###
    regional_summary_cohort$Comparison <- as.character(regional_summary_cohort$Comparison) 
    #############################################################################################################################################
    
    #############################################################################################################################################
    my_plot[[i]] <- ggplot(regional_summary) + 
      aes(x=count, y=Comparison, fill=logFC_dir, label=abs(count)) + 
      geom_bar(stat='identity', width=0.5) +
      geom_text(data=regional_summary[regional_summary$count != 0,], show.legend=FALSE, color='black') +
      theme_classic() + theme(axis.text.y = element_text(face='bold'),
                              axis.text.x = element_text(face='bold'),
                              axis.title = element_text(face = "bold"),
                              legend.title= element_text(face='bold'),
                              strip.text=element_text(face='bold'),
                              plot.title=element_text(face='bold')) + 
      scale_fill_manual(values=c('firebrick','steelblue'), name='LogFC Direction') +
      labs(x='no. Markers (FDR < 0.05)', y='Contrast') + theme(legend.position="bottom")
    #############################################################################################################################################
    
    #############################################################################################################################################
    my_plot_cohort[[i]] <- ggplot(regional_summary_cohort) + 
      aes(x=count, y=Comparison, fill=logFC_dir, label=abs(count)) + 
      geom_bar(stat='identity', width=0.5) +
      geom_text(data=regional_summary_cohort[regional_summary_cohort$count != 0,], show.legend=FALSE, color='black') +
      theme_classic() + theme(axis.text.y = element_text(face='bold'),
                              axis.text.x = element_text(face='bold'),
                              axis.title = element_text(face = "bold"),
                              legend.title= element_text(face='bold'),
                              strip.text=element_text(face='bold'),
                              plot.title=element_text(face='bold')) + 
      scale_fill_manual(values=c('firebrick','steelblue'), name='LogFC Direction') +
      labs(x='no. Markers (FDR < 0.05)', y='Contrast') + theme(legend.position="bottom")
    #############################################################################################################################################
    pdf(file=paste(levels(as.factor(channel_intesity_statistics$CellSubset))[i],"heatmap","pdf",sep="."), width = 14, height = 8)
    pheatmap(t(generic_var[,]),
             scale="row", 
             main=levels(as.factor(channel_intesity_statistics$CellSubset))[i],
             clustering_distance_rows ="euclidean", 
             show_rownames = TRUE, 
             show_colnames = FALSE, 
             annotation_col = annotation_row[,], 
             annotation_row = annotation_col, 
             annotation_colors = annotation_colors, 
             color = colorRampPalette(c("blue","white","orangered"))(255), 
             angle_col = 90,
             border_color = "grey70",
             cluster_cols = F, 
             cluster_rows = T)
    dev.off()
    #############################################################################################################################################
    
    #############################################################################################################################################
    pdf(file=paste(levels(as.factor(channel_intesity_statistics$CellSubset))[i],"cohort","heatmap","pdf",sep="."), width = 14, height = 8)
    pheatmap(sva_x_matrix,
             scale="row", 
             main=levels(as.factor(channel_intesity_statistics$CellSubset))[i],
             clustering_distance_rows ="euclidean", 
             show_rownames = TRUE, 
             show_colnames = FALSE, 
             annotation_col = annotation_row[,], 
             annotation_row = annotation_col, 
             annotation_colors = annotation_colors, 
             color = colorRampPalette(c("blue","white","orangered"))(255), 
             angle_col = 90,
             border_color = "grey70",
             cluster_cols = F, 
             cluster_rows = F)
    dev.off()
    #############################################################################################################################################
  }
  
}

#############################################################################################################################################
### Long loop to calculate the DEs per celltype [ END ]
#############################################################################################################################################

names(lm_cytof_markers) <- levels(as.factor(channel_intesity_statistics$CellSubset))
names(lm_cytof_markers_cohort) <- levels(as.factor(channel_intesity_statistics$CellSubset))
dev.off()

#############################################################################################################################################
### complete figures
for (i in 1:length(levels(as.factor(channel_intesity_statistics$CellSubset)))) {
  pdf(file=paste(levels(as.factor(channel_intesity_statistics$CellSubset))[i],"DE_comparisons","pdf",sep="."),width = 4, height = 6)
  print(my_plot[[i]])
  dev.off() }
###
for (i in 1:length(levels(as.factor(channel_intesity_statistics$CellSubset)))) {
  pdf(file=paste(levels(as.factor(channel_intesity_statistics$CellSubset))[i],"cohort","DE_comparisons","pdf",sep="."),width = 4, height = 6)
  print(my_plot_cohort[[i]])
  dev.off() }
#############################################################################################################################################
for (i in 1:length( lm_cytof_markers )){ lm_cytof_markers[[i]]$CellSubset <- names(lm_cytof_markers)[i] }
lm_cytof_markers_all <- do.call(rbind,lm_cytof_markers)
###
for (i in 1:length( lm_cytof_markers_cohort )){ lm_cytof_markers_cohort[[i]]$CellSubset <- names(lm_cytof_markers_cohort)[i] }
lm_cytof_markers_cohort_all <- do.call(rbind,lm_cytof_markers_cohort)
###
lm_cytof_markers_all <- lm_cytof_markers_all[lm_cytof_markers_all$adj.P.Val<0.05,]
lm_cytof_markers_cohort_all <- lm_cytof_markers_cohort_all[lm_cytof_markers_cohort_all$adj.P.Val<0.05,]
lm_olink_all <- lm_olink_all[lm_olink_all$adj.P.Val<0.05,]
lm_frequencies_all <- lm_frequencies_all[lm_frequencies_all$adj.P.Val<0.05,]
###
rownames(lm_cytof_markers_all) <- NULL
rownames(lm_cytof_markers_cohort_all) <- NULL
rownames(lm_olink_all) <- NULL
rownames(lm_frequencies_all) <- NULL
rm(my_plot,my_plot_cohort)
#############################################################################################################################################
setwd("/Users/gonzae34/Documents/projects_gnjatic/U24/U24/E4412_Diefenbach")
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
### Summary
### Various Figures
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

#############################################################################################################################################
### Summary figure for O-LINK ; LogFC & FDR ; within cohort 
#############################################################################################################################################
generic_var <- lm_olink_all[lm_olink_all$adj.P.Val<0.05,]
# levels(as.factor(generic_var$Comparison))
generic_var$Comparison <- as.character(generic_var$Comparison)
###
ix <- generic_var$Comparison %in%  c("C3vsC3_Base_Res","C2vsC2_Base_pre","C2vsC2_Base_Res","C2vsC2_Base_OS","C1vsC1_OS_Pre",
                                     "C3vsC3_Base_pre","C1vsC1_Base_OS","C3vsC3_Pre_OS" ,"C3vsC3_Base_OS" ,"C1vsC1_Base_pre","C1vsC1_Base_Res")
###
generic_var <- generic_var[ix,]
###
generic_var <- dcast(data = generic_var,formula = Comparison~protein, value.var = "logFC")
generic_var[is.na(generic_var)]<-0
rownames(generic_var) <- generic_var$Comparison
###
# ord <- hclust( dist(scale((generic_var[,-1])), method = "euclidean"), method = "ward.D" )$order
# ord
out <- pheatmap(generic_var[,-1])
dev.off()
###
generic_var <- lm_olink_all[lm_olink_all$adj.P.Val<0.05,]
# levels(as.factor(generic_var$Comparison))
generic_var$Comparison <- as.character(generic_var$Comparison)
generic_var <- generic_var[ix,]
generic_var$protein <- factor( generic_var$protein, levels = out$tree_col$labels[out$tree_col$order] )
generic_var$Comparison <- factor( generic_var$Comparison, levels = out$tree_row$labels[out$tree_row$order] )
generic_var$nLogFDR <- -log10(generic_var$adj.P.Val)
###
pdf(file="figures/olink_summary_of_results_logFC_within.pdf",width = 5,height = 4)
ggplot(generic_var, aes(Comparison, protein, color= logFC, size=nLogFDR )) + 
  geom_point(shape=15)+
  scale_colour_gradient2(low='firebrick', mid="white",high='steelblue', name='LogFC') +
  theme_bw() +
  rotate_x_text(angle=45) + 
  labs(x ='', y='', title='Olink: Within Cohorts') ### +scale_fill_viridis(name='LogFC')
dev.off()
###
pdf(file="figures/olink_summary_of_results_logFC_within_sacha_inv.pdf",width = 5,height = 4)
ggplot(generic_var, aes(Comparison, protein, color= -logFC, size=nLogFDR )) + 
  geom_point(shape=15)+
  scale_colour_gradient2(low='firebrick', mid="white",high='steelblue', name='LogFC') +
  theme_bw() +
  rotate_x_text(angle=45) + 
  labs(x ='', y='', title='Olink: Within Cohorts') ### +scale_fill_viridis(name='LogFC')
dev.off()
#############################################################################################################################################
### Summary figure for O-LINK ; FDR ; between cohort 
#############################################################################################################################################
generic_var <- lm_olink_all[lm_olink_all$adj.P.Val<0.05,]
# levels(as.factor(generic_var$Comparison))
generic_var$Comparison <- as.character(generic_var$Comparison)
###
ix <- generic_var$Comparison %in% c("C1vsC2_B","C2vsC3_B","C1vsC3_Res","C1vsC3_OS","C1vsC3_Pre","C1vsC2_OS","C1vsC2_Pre","C1vsC2_Res",
                                    "BaselinevsC1_Pre","BaselinevsC1_Res","BaselinevsC2_Pre","BaselinevsC2_Res","BaselinevsC3_Pre","BaselinevsC3_Res")
generic_var <- generic_var[ix,]
###
generic_var <- dcast(data = generic_var,formula = Comparison~protein, value.var = "logFC")
generic_var[is.na(generic_var)]<-0
rownames(generic_var) <- generic_var$Comparison
###
# ord <- hclust( dist(scale((generic_var[,-1])), method = "euclidean"), method = "ward.D" )$order
# ord
out <- pheatmap(generic_var[,-1])
dev.off()
###
generic_var <- lm_olink_all[lm_olink_all$adj.P.Val<0.05,]
# levels(as.factor(generic_var$Comparison))
generic_var$Comparison <- as.character(generic_var$Comparison)
generic_var <- generic_var[ix,]
generic_var$protein <- factor( generic_var$protein, levels = out$tree_col$labels[out$tree_col$order] )
generic_var$Comparison <- factor( generic_var$Comparison, levels = out$tree_row$labels[out$tree_row$order] )
###
generic_var$nLogFDR <- -log10(generic_var$adj.P.Val)
###
pdf(file="figures/olink_summary_of_results_logFC_between.pdf",width = 6,height = 4.5)
ggplot(generic_var, aes(Comparison, protein, color= logFC, size=nLogFDR )) + 
  geom_point(shape=15)+
  scale_colour_gradient2(low='firebrick', mid="white",high='steelblue', name='LogFC') +
  theme_bw() +
  rotate_x_text(angle=45) + 
  labs(x ='', y='', title='Olink: Between Cohorts') ### +scale_fill_viridis(name='LogFC')
dev.off()
###
pdf(file="figures/olink_summary_of_results_logFC_between_sacha_inv.pdf",width = 6,height = 4.5)
ggplot(generic_var, aes(Comparison, protein, color= -logFC, size=nLogFDR )) + 
  geom_point(shape=15)+
  scale_colour_gradient2(low='firebrick', mid="white",high='steelblue', name='LogFC') +
  theme_bw() +
  rotate_x_text(angle=45) + 
  labs(x ='', y='', title='Olink: Between Cohorts') ### +scale_fill_viridis(name='LogFC')
dev.off()
###
#############################################################################################################################################

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

#############################################################################################################################################
### Results for number of markers, CYTOF within cohort
#############################################################################################################################################
generic_var <- lm_cytof_markers_all [ lm_cytof_markers_all$adj.P.Val < 0.05 & abs(lm_cytof_markers_all$logFC) > 0.25, ]
##
generic_var$Marker <- gsub("_","",generic_var$Marker)
generic_var$CellSubset <- gsub("_"," ",generic_var$CellSubset)
## generic_var
generic_var$sorter <- paste(generic_var$Marker,generic_var$CellSubset,sep="_")
##
my_levels <- levels(as.factor(generic_var$sorter))
##
my_summary_matrix<-list()
for( i in 1:length(my_levels)){ ii <- which(generic_var$sorter %in% my_levels[i])
my_summary_matrix[[i]] <- colMeans(generic_var[ii,1:5]) 
my_summary_matrix[[i]]$sorter <- my_levels[i]
my_summary_matrix[[i]]$ncomp <- length(ii)
}
##
my_summary_matrix <- do.call(rbind,my_summary_matrix)
my_summary_matrix <- as.data.frame(my_summary_matrix)
##
my_summary_matrix$logFC <- unlist(my_summary_matrix$logFC)
my_summary_matrix$AveExpr <- unlist(my_summary_matrix$AveExpr)
my_summary_matrix$t <- unlist(my_summary_matrix$t)
my_summary_matrix$P.Value <- unlist(my_summary_matrix$P.Value)
my_summary_matrix$adj.P.Val <- unlist(my_summary_matrix$adj.P.Val)
my_summary_matrix$sorter <- unlist(my_summary_matrix$sorter)
my_summary_matrix$ncomp <- unlist(my_summary_matrix$ncomp)
##
my_summary_matrix$Marker <- tidyr::separate(data.frame(as.character(my_summary_matrix$sorter)), 1, sep="_", c("a","b","c"))$a
my_summary_matrix$Celltype <- tidyr::separate(data.frame(as.character(my_summary_matrix$sorter)), 1, sep="_", c("a","b","c"))$b
###
my_summary_matrix$nLogFDR <- -log10(my_summary_matrix$adj.P.Val)
my_summary_matrix$nLogFDR [my_summary_matrix$nLogFDR > quantile(my_summary_matrix$nLogFDR)[4]] <- 7.5
###
astrolabe_markers_v2 <- gsub("_","",astrolabe_markers) 
my_summary_matrix <- my_summary_matrix[ my_summary_matrix$Marker %in% astrolabe_markers_v2 , ]
###
pdf(file="figures/cytof_summary_of_results_within.pdf",width = 5.5,height = 5)
ggplot(my_summary_matrix, aes(Celltype, Marker ,color=nLogFDR,size=ncomp )) + 
  geom_point(shape=15) + 
  theme_bw() +
  rotate_x_text(angle=65) +
  labs(x ='', y='', title='DE Markers per Celltype, Within cohort') +
  scale_color_viridis(name='nLogFDR')  
dev.off()
###  
#############################################################################################################################################
#############################################################################################################################################
### Results for number of markers, CYTOF between cohort
#############################################################################################################################################
generic_var <- lm_cytof_markers_cohort_all [ lm_cytof_markers_cohort_all$adj.P.Val < 0.05 & abs(lm_cytof_markers_cohort_all$logFC) > 0.25, ]
##
generic_var$Marker <- gsub("_","",generic_var$Marker)
generic_var$CellSubset <- gsub("_"," ",generic_var$CellSubset)
## generic_var
generic_var$sorter <- paste(generic_var$Marker,generic_var$CellSubset,sep="_")
##
my_levels <- levels(as.factor(generic_var$sorter))
##
my_summary_matrix<-list()
for( i in 1:length(my_levels)){ ii <- which(generic_var$sorter %in% my_levels[i])
my_summary_matrix[[i]] <- colMeans(generic_var[ii,1:5]) 
my_summary_matrix[[i]]$sorter <- my_levels[i]
my_summary_matrix[[i]]$ncomp <- length(ii)
}
##
my_summary_matrix <- do.call(rbind,my_summary_matrix)
my_summary_matrix <- as.data.frame(my_summary_matrix)
##
my_summary_matrix$logFC <- unlist(my_summary_matrix$logFC)
my_summary_matrix$AveExpr <- unlist(my_summary_matrix$AveExpr)
my_summary_matrix$t <- unlist(my_summary_matrix$t)
my_summary_matrix$P.Value <- unlist(my_summary_matrix$P.Value)
my_summary_matrix$adj.P.Val <- unlist(my_summary_matrix$adj.P.Val)
my_summary_matrix$sorter <- unlist(my_summary_matrix$sorter)
my_summary_matrix$ncomp <- unlist(my_summary_matrix$ncomp)
##
my_summary_matrix$Marker <- tidyr::separate(data.frame(as.character(my_summary_matrix$sorter)), 1, sep="_", c("a","b","c"))$a
my_summary_matrix$Celltype <- tidyr::separate(data.frame(as.character(my_summary_matrix$sorter)), 1, sep="_", c("a","b","c"))$b
###
my_summary_matrix$nLogFDR <- -log10(my_summary_matrix$adj.P.Val)
my_summary_matrix$nLogFDR [my_summary_matrix$nLogFDR > quantile(my_summary_matrix$nLogFDR)[4]] <- 7.5
###
astrolabe_markers_v2 <- gsub("_","",astrolabe_markers) 
my_summary_matrix <- my_summary_matrix[ my_summary_matrix$Marker %in% astrolabe_markers_v2 , ]
###
pdf(file="figures/cytof_summary_of_results_between.pdf",width = 5.5,height = 5)
ggplot(my_summary_matrix, aes(Celltype, Marker ,color=nLogFDR,size=ncomp )) + 
  geom_point(shape=15) + 
  theme_bw() +
  rotate_x_text(angle=65) +
  labs(x ='', y='', title='DE Markers per Celltype, Between cohort') +
  scale_color_viridis(name='nLogFDR')  
dev.off()
###  
#############################################################################################################################################


#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
### loop for figures
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

# levels(as.factor(lm_cytof_markers_all$CellSubset))
my_cells <- levels(as.factor(lm_cytof_markers_all$CellSubset))[c(1,3,4,
                                                                 5,6,7,8,9,10,11,
                                                                 12,13,14,15,
                                                                 18,19,
                                                                 21,
                                                                 22,23,24,
                                                                 25,26,27)]

#############################################################################################################################################
#############################################################################################################################################
setwd("/Users/gonzae34/Documents/projects_gnjatic/U24/U24/E4412_Diefenbach/figures/celltype_cytof_marker_summaries/")

a<-list()
b<-list()
c<-list()
d<-list()
#############################################################################################################################################
for ( i in 1:length(my_cells)) {
  #############################################################################################################################################
  
  #############################################################################################################################################
  ii <- lm_cytof_markers_all$CellSubset %in% my_cells[i]
  ###
  generic_var <- dcast(data = lm_cytof_markers_all[ii,],formula = Comparison~Marker, value.var = "logFC")
  generic_var[is.na(generic_var)]<-0
  rownames(generic_var) <- generic_var$Comparison
  ###
  out <- pheatmap(generic_var[,-1])
  dev.off()
  ###
  generic_var <- lm_cytof_markers_all[ii,]
  generic_var$Marker <- factor( generic_var$Marker, levels = out$tree_col$labels[out$tree_col$order] )
  generic_var$Comparison <- factor( generic_var$Comparison, levels = out$tree_row$labels[out$tree_row$order] )
  ###
  a[[i]] <- ggplot(generic_var, aes(Comparison, Marker, fill= logFC, )) + 
    geom_tile(color='white') + theme_bw() +
    rotate_x_text(angle=45) +
    labs(x ='', y='', title=paste(my_cells[i],"Proteins with FDR<0.05",sep=" ") ) +
    scale_fill_viridis(name='LogFC')
  #############################################################################################################################################
  
  
  #############################################################################################################################################
  ii <- lm_cytof_markers_all$CellSubset %in% my_cells[i]
  ###
  generic_var <- dcast(data = lm_cytof_markers_all[ii,],formula = Comparison~Marker, value.var = "adj.P.Val")
  generic_var[is.na(generic_var)]<-0
  rownames(generic_var) <- generic_var$Comparison
  ###
  out <- pheatmap(generic_var[,-1])
  dev.off()
  ###
  generic_var <- lm_cytof_markers_all[ii,]
  generic_var$Marker <- factor( generic_var$Marker, levels = out$tree_col$labels[out$tree_col$order] )
  generic_var$Comparison <- factor( generic_var$Comparison, levels = out$tree_row$labels[out$tree_row$order] )
  ###
  b[[i]] <- ggplot(generic_var, aes(Comparison, Marker, fill= -log10(adj.P.Val+1e-200), )) + 
    geom_tile(color='white') + theme_bw() +
    rotate_x_text(angle=45) +
    labs(x ='', y='', title=paste(my_cells[i],"Proteins with FDR<0.05",sep=" ") ) +
    scale_fill_viridis(name='-Log10(FDR)')
  #############################################################################################################################################
  
  
  #############################################################################################################################################
  ii <- lm_cytof_markers_cohort_all$CellSubset %in% my_cells[i]
  ###
  generic_var <- dcast(data = lm_cytof_markers_cohort_all[ii,],formula = Comparison~Marker, value.var = "logFC")
  generic_var[is.na(generic_var)]<-0
  rownames(generic_var) <- generic_var$Comparison
  ###
  out <- pheatmap(generic_var[,-1])
  dev.off()
  ###
  generic_var <- lm_cytof_markers_cohort_all[ii,]
  generic_var$Marker <- factor( generic_var$Marker, levels = out$tree_col$labels[out$tree_col$order] )
  generic_var$Comparison <- factor( generic_var$Comparison, levels = out$tree_row$labels[out$tree_row$order] )
  ###
  c[[i]] <- ggplot(generic_var, aes(Comparison, Marker, fill= logFC, )) + 
    geom_tile(color='white') + theme_bw() +
    rotate_x_text(angle=45) +
    labs(x ='', y='', title=paste(my_cells[i],"Proteins with FDR<0.05",sep=" ") ) +
    scale_fill_viridis(name='LogFC')
  #############################################################################################################################################
  
  
  #############################################################################################################################################
  ii <- lm_cytof_markers_cohort_all$CellSubset %in% my_cells[i]
  ###
  generic_var <- dcast(data = lm_cytof_markers_cohort_all[ii,],formula = Comparison~Marker, value.var = "adj.P.Val")
  generic_var[is.na(generic_var)]<-0
  rownames(generic_var) <- generic_var$Comparison
  ###
  out <- pheatmap(generic_var[,-1])
  dev.off()
  ###
  generic_var <- lm_cytof_markers_cohort_all[ii,]
  generic_var$Marker <- factor( generic_var$Marker, levels = out$tree_col$labels[out$tree_col$order] )
  generic_var$Comparison <- factor( generic_var$Comparison, levels = out$tree_row$labels[out$tree_row$order] )
  ###
  d[[i]] <- ggplot(generic_var, aes(Comparison, Marker, fill= -log10(adj.P.Val+1e-200), )) + 
    geom_tile(color='white') + theme_bw() +
    rotate_x_text(angle=45) +
    labs(x ='', y='', title=paste(my_cells[i],"Proteins with FDR<0.05",sep=" ") ) +
    scale_fill_viridis(name='-Log10(FDR)')
  #############################################################################################################################################
  
}

#############################################################################################################################################
for( i in 1:length(my_cells) ){
  pdf(file=paste("Cytof_summary",my_cells[i],"within_cohort_logFC","pdf",sep="."), width = 8, height = 6)
  print(a[[i]])
  dev.off()
  
  pdf(file=paste("Cytof_summary",my_cells[i],"within_cohort_FDR","pdf",sep="."), width = 8, height = 6)
  print(b[[i]])
  dev.off()
  
  pdf(file=paste("Cytof_summary",my_cells[i],"between_cohort_logFC","pdf",sep="."), width = 8, height = 6)
  print(c[[i]])
  dev.off()
  
  pdf(file=paste("Cytof_summary",my_cells[i],"between_cohort_FDR","pdf",sep="."), width = 8, height = 6)
  print(d[[i]])
  dev.off() }
#############################################################################################################################################
#############################################################################################################################################
### end of loop
#############################################################################################################################################

### LEGACY ###
#############################################################################################################################################
### Results for number of markers, cytof within cohort ### LEGACY ###
#############################################################################################################################################
generic_var <- lm_cytof_markers_all[lm_cytof_markers_all$Marker %in% astrolabe_markers,]
generic_var <- table(generic_var$CellSubset,generic_var$Comparison)
class(generic_var) <- "matrix"
generic_var <- as.data.frame(generic_var)
generic_var$CellSubset <- rownames(generic_var)
generic_var <- generic_var[ ! generic_var$CellSubset %in% "Root_unassigned" , ]
generic_var <- generic_var[ -grep("unassigned",generic_var$CellSubset) , ]
###
out <- pheatmap(generic_var[,-ncol(generic_var)])
dev.off()
###
generic_var <- melt(generic_var)
generic_var$CellSubset <- factor( generic_var$CellSubset, levels = out$tree_row$labels[out$tree_row$order]  )
generic_var$variable <- factor( generic_var$variable, levels = out$tree_col$labels[out$tree_col$order] )
###
#pdf(file="figures/cytof_summary_of_results_within_astrolabe_markers.pdf",width = 8,height = 6)
ggplot(generic_var, aes(variable, CellSubset, fill= value )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof Summary, DE Marker per Celltype, Within cohort') +
  scale_fill_viridis(name='# of DE Markers') +
  theme(legend.position = "bottom")
#dev.off()
###  
#############################################################################################################################################
#############################################################################################################################################
### Results for number of markers, cytof between cohort ### LEGACY ###
#############################################################################################################################################
generic_var <- lm_cytof_markers_cohort_all[lm_cytof_markers_cohort_all$Marker %in% astrolabe_markers,]
generic_var <- table(generic_var$CellSubset,generic_var$Comparison)
class(generic_var) <- "matrix"
generic_var <- as.data.frame(generic_var)
generic_var$CellSubset <- rownames(generic_var)
generic_var <- generic_var[ ! generic_var$CellSubset %in% "Root_unassigned" , ]
generic_var <- generic_var[ -grep("unassigned",generic_var$CellSubset) , ]
###
out <- pheatmap(generic_var[,-ncol(generic_var)])
dev.off()
###
generic_var <- melt(generic_var)
generic_var$CellSubset <- factor( generic_var$CellSubset, levels = out$tree_row$labels[out$tree_row$order]  )
generic_var$variable <- factor( generic_var$variable, levels = out$tree_col$labels[out$tree_col$order] )
###
#pdf(file="figures/cytof_summary_of_results_between_astrolabe_markers.pdf",width = 8,height = 6)
ggplot(generic_var, aes(variable, CellSubset, fill= value )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof Summary, DE Marker per Celltype, Between cohorts') +
  scale_fill_viridis(name='# of DE Markers') +
  theme(legend.position = "bottom")
#dev.off()
###  
#############################################################################################################################################


my_cells
#############################################################################################################################################
### Monocyte (CD14+ CD16+)
#############################################################################################################################################
table(lm_cytof_markers_all$CellSubset %in% "Monocyte (CD14+ CD16+)")
ii <- lm_cytof_markers_all$CellSubset %in% "Monocyte (CD14+ CD16+)"
###
generic_var <- dcast(data = lm_cytof_markers_all[ii,],formula = Comparison~Marker, value.var = "logFC")
generic_var[is.na(generic_var)]<-0
rownames(generic_var) <- generic_var$Comparison
generic_var <- generic_var [ , colnames(generic_var) %in% astrolabe_markers ]
###
out <- pheatmap(generic_var)
dev.off()
###
generic_var <- lm_cytof_markers_all[ii,]
generic_var <- generic_var [ generic_var$Marker %in% astrolabe_markers , ]
###
generic_var$Marker <- factor( generic_var$Marker, levels = out$tree_col$labels[out$tree_col$order] )
generic_var$Comparison <- factor( generic_var$Comparison, levels = out$tree_row$labels[out$tree_row$order] )
###
pdf(file="figures/cytof_summary_of_results_monocytes_cd14p_cd16p_logfc.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= logFC, )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, Monocyte (CD14+ CD16+)\nProteins with FDR<0.05') +
  scale_fill_viridis(name='LogFC')
dev.off()
###
pdf(file="figures/cytof_summary_of_results_monocytes_cd14p_cd16p_fdr.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= -log10(adj.P.Val), )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, Monocyte (CD14+ CD16+)\nProteins with FDR<0.05') +
  scale_fill_viridis(name='-Log10(FDR)')
dev.off()
#############################################################################################################################################
### Monocyte (CD14+ CD16+)
#############################################################################################################################################
table(lm_cytof_markers_cohort_all$CellSubset %in% "Monocyte (CD14+ CD16+)")
ii <- lm_cytof_markers_cohort_all$CellSubset %in% "Monocyte (CD14+ CD16+)"
###
generic_var <- dcast(data = lm_cytof_markers_cohort_all[ii,],formula = Comparison~Marker, value.var = "logFC")
generic_var[is.na(generic_var)]<-0
rownames(generic_var) <- generic_var$Comparison
generic_var <- generic_var [ , colnames(generic_var) %in% astrolabe_markers ]
###
out <- pheatmap(generic_var)
dev.off()
###
generic_var <- lm_cytof_markers_cohort_all[ii,]
generic_var <- generic_var [ generic_var$Marker %in% astrolabe_markers , ]
###
generic_var$Marker <- factor( generic_var$Marker, levels = out$tree_col$labels[out$tree_col$order] )
generic_var$Comparison <- factor( generic_var$Comparison, levels = out$tree_row$labels[out$tree_row$order] )
###
pdf(file="figures/cytof_summary_of_results_monocytes_cd14p_cd16p_logfc_cohort.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= logFC, )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, Monocyte (CD14+ CD16+)\nProteins with FDR<0.05') +
  scale_fill_viridis(name='LogFC')
dev.off()
###
pdf(file="figures/cytof_summary_of_results_monocytes_cd14p_cd16p_fdr_cohort.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= -log10(adj.P.Val), )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, Monocyte (CD14+ CD16+)\nProteins with FDR<0.05') +
  scale_fill_viridis(name='-Log10(FDR)')
dev.off()
#############################################################################################################################################



#############################################################################################################################################
### B Cell (CD27-)
#############################################################################################################################################
table(lm_cytof_markers_all$CellSubset %in% "B Cell (CD27-)")
ii <- lm_cytof_markers_all$CellSubset %in% "B Cell (CD27-)"
###
generic_var <- dcast(data = lm_cytof_markers_all[ii,],formula = Comparison~Marker, value.var = "logFC")
generic_var[is.na(generic_var)]<-0
rownames(generic_var) <- generic_var$Comparison
generic_var <- generic_var [ , colnames(generic_var) %in% astrolabe_markers ]
###
out <- pheatmap(generic_var)
dev.off()
###
generic_var <- lm_cytof_markers_all[ii,]
generic_var <- generic_var [ generic_var$Marker %in% astrolabe_markers , ]
###
generic_var$Marker <- factor( generic_var$Marker, levels = out$tree_col$labels[out$tree_col$order] )
generic_var$Comparison <- factor( generic_var$Comparison, levels = out$tree_row$labels[out$tree_row$order] )
###
pdf(file="figures/cytof_summary_of_results_bcell_cd27n_logfc.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= logFC, )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, B Cell (CD27-)\nProteins with FDR<0.05') +
  scale_fill_viridis(name='LogFC')
dev.off()
###
pdf(file="figures/cytof_summary_of_results_bcell_cd27n_fdr.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= -log10(adj.P.Val), )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, B Cell (CD27-)\nProteins with FDR<0.05') +
  scale_fill_viridis(name='-Log10(FDR)')
dev.off()
#############################################################################################################################################
### B Cell (CD27-)
#############################################################################################################################################
table(lm_cytof_markers_cohort_all$CellSubset %in% "B Cell (CD27-)")
ii <- lm_cytof_markers_cohort_all$CellSubset %in% "B Cell (CD27-)"
###
generic_var <- dcast(data = lm_cytof_markers_cohort_all[ii,],formula = Comparison~Marker, value.var = "logFC")
generic_var[is.na(generic_var)]<-0
rownames(generic_var) <- generic_var$Comparison
generic_var <- generic_var [ , colnames(generic_var) %in% astrolabe_markers ]
###
out <- pheatmap(generic_var)
dev.off()
###
generic_var <- lm_cytof_markers_cohort_all[ii,]
generic_var <- generic_var [ generic_var$Marker %in% astrolabe_markers , ]
###
generic_var$Marker <- factor( generic_var$Marker, levels = out$tree_col$labels[out$tree_col$order] )
generic_var$Comparison <- factor( generic_var$Comparison, levels = out$tree_row$labels[out$tree_row$order] )
###
pdf(file="figures/cytof_summary_of_results_bcell_cd27n_logfc_cohort.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= logFC, )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, B Cell (CD27-)\nProteins with FDR<0.05') +
  scale_fill_viridis(name='LogFC')
dev.off()
###
pdf(file="figures/cytof_summary_of_results_bcell_cd27n_fdr_cohort.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= -log10(adj.P.Val), )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, B Cell (CD27-)\nProteins with FDR<0.05') +
  scale_fill_viridis(name='-Log10(FDR)')
dev.off()
#############################################################################################################################################
#############################################################################################################################################
### 
#############################################################################################################################################
ii <- channel_intesity_statistics$CellSubset %in% "B Cell (CD27-)"
###
generic_var <- dcast(data = channel_intesity_statistics[ii,],formula = SampleId~ChannelName, value.var = "Q0_95")
###
rownames(generic_var) <- generic_var$SampleId
generic_var <- generic_var [ , colnames(generic_var) %in% unique(sort(lm_cytof_markers_all$Marker[lm_cytof_markers_all$CellSubset %in% "B Cell (CD27-)"])) ]
###
generic_var$astrolabe_ID <- rownames(generic_var)
generic_var <- (melt(generic_var))
colnames(generic_var)[colnames(generic_var) %in% "variable"] <- "Marker"
generic_var <- merge(generic_var,cimac_e4412_metadata,by="astrolabe_ID")
###
ii <- c(grep("Baseline",generic_var$Condition),grep("Restaging",generic_var$Condition),
        grep("Pre_Day",generic_var$Condition),grep("End_of_",generic_var$Condition))
##
generic_var <- generic_var[ii,]
###
table( names(table(generic_var$condition)) %in% c( "Baseline_1_BV+ipi","Pre_Day_1_Cycle_2_1_BV+ipi","Restaging_1_BV+ipi","End_of_treatment_1_BV+ipi",
                                                   "Baseline_2_BV+nivo","Pre_Day_1_Cycle_2_2_BV+nivo ","Restaging_2_BV+nivo","End_of_treatment_2_BV+nivo ",
                                                   "Baseline_3_BV+nivo+ipi","Pre_Day_1_Cycle_2_3_BV+nivo+ipi","Restaging_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi" ) )

names(table(generic_var$condition))[ ! names(table(generic_var$condition)) %in% c( "Baseline_1_BV+ipi","Pre_Day_1_Cycle_2_1_BV+ipi","Restaging_1_BV+ipi","End_of_treatment_1_BV+ipi",
                                                                                   "Baseline_2_BV+nivo","Pre_Day_1_Cycle_2_2_BV+nivo","Restaging_2_BV+nivo","End_of_treatment_2_BV+nivo",
                                                                                   "Baseline_3_BV+nivo+ipi","Pre_Day_1_Cycle_2_3_BV+nivo+ipi","Restaging_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi" )]

generic_var$condition <- factor(generic_var$condition, levels = c( "Baseline_1_BV+ipi","Pre_Day_1_Cycle_2_1_BV+ipi","Restaging_1_BV+ipi","End_of_treatment_1_BV+ipi",
                                                                   "Baseline_2_BV+nivo","Pre_Day_1_Cycle_2_2_BV+nivo","Restaging_2_BV+nivo","End_of_treatment_2_BV+nivo",
                                                                   "Baseline_3_BV+nivo+ipi","Pre_Day_1_Cycle_2_3_BV+nivo+ipi","Restaging_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi" ))
###
generic_var$condition <- fct_rev(generic_var$condition)
###
my_comparisons <- list( #c("Baseline_1_BV+ipi","Pre_Day_1_Cycle_2_1_BV+ipi"),
  #c("Baseline_1_BV+ipi","Restaging_1_BV+ipi"),
  #c("Baseline_1_BV+ipi","End_of_treatment_1_BV+ipi"),
  #c("Restaging_1_BV+ipi","End_of_treatment_1_BV+ipi"),
  c("End_of_treatment_2_BV+nivo","End_of_treatment_1_BV+ipi"),
  #c("End_of_treatment_2_BV+nivo","End_of_treatment_3_BV+nivo+ipi"),
  c("End_of_treatment_1_BV+ipi","End_of_treatment_3_BV+nivo+ipi"),
  #c("Restaging_1_BV+ipi","Restaging_2_BV+nivo"),
  #c("Restaging_3_BV+nivo+ipi","Restaging_2_BV+nivo"),
  c("Restaging_1_BV+ipi","Restaging_3_BV+nivo+ipi")
  #c("Baseline_2_BV+nivo","Pre_Day_1_Cycle_2_2_BV+nivo"),
  #c("Baseline_2_BV+nivo","Restaging_2_BV+nivo"),
  #c("Baseline_2_BV+nivo","End_of_treatment_2_BV+nivo"),
  #c("Restaging_2_BV+nivo","End_of_treatment_2_BV+nivo"),
  #c("Baseline_3_BV+nivo+ipi","Pre_Day_1_Cycle_2_3_BV+nivo+ipi"),
  #c("Baseline_3_BV+nivo+ipi","Restaging_3_BV+nivo+ipi"),
  #c("Baseline_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi"),
  #c("Restaging_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi")
)
#############################################################################################################################################

pdf(file="figures/examples_b_cell_cd27n_cytof_markers_diff_expr.pdf",width = 7, height = 4)

ii <- generic_var$Marker %in% "TIGIT" & generic_var$value < 2
test <- generic_var[ii,] ; test <- test[!is.na(test$condition),]

ggboxplot(data=test, 
          x="condition", y="value",fill="Cohort") + theme_minimal() +
  # geom_boxplot(width=0.1) + facet_wrap(~protein,ncol=4) +
  # stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("TIGIT") + labs(x = "", y = "Marker Exp. Values") + theme(legend.position="none")

ii <- generic_var$Marker %in% "CCR4"  & generic_var$value < 2
test <- generic_var[ii,] ; test <- test[!is.na(test$condition),]

ggboxplot(data=test,
          x="condition", y="value",fill="Cohort") + theme_minimal() +
  #geom_boxplot(width=0.1) + facet_wrap(~protein,ncol=4) +
  #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("CCR4") + labs(x = "", y = "Marker Exp. Values") + theme(legend.position="none")

ii <- generic_var$Marker %in% "CD25" & generic_var$value < 2
test <- generic_var[ii,] ; test <- test[!is.na(test$condition),]

ggboxplot(data=test,
          x="condition", y="value",fill="Cohort") + theme_minimal() +
  #geom_boxplot(width=0.1) + facet_wrap(~protein,ncol=4) +
  #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("CD25") + labs(x = "", y = "Marker Exp. Values") + theme(legend.position="none")

dev.off()

#############################################################################################################################################



#############################################################################################################################################
### B Cell (Plasmablast)
#############################################################################################################################################
table(lm_cytof_markers_all$CellSubset %in% "B Cell (Plasmablast)")
ii <- lm_cytof_markers_all$CellSubset %in% "B Cell (Plasmablast)"
###
generic_var <- dcast(data = lm_cytof_markers_all[ii,],formula = Comparison~Marker, value.var = "logFC")
generic_var[is.na(generic_var)]<-0
rownames(generic_var) <- generic_var$Comparison
generic_var <- generic_var [ , colnames(generic_var) %in% astrolabe_markers ]
###
out <- pheatmap(generic_var)
dev.off()
###
generic_var <- lm_cytof_markers_all[ii,]
generic_var <- generic_var [ generic_var$Marker %in% astrolabe_markers , ]
###
generic_var$Marker <- factor( generic_var$Marker, levels = out$tree_col$labels[out$tree_col$order] )
generic_var$Comparison <- factor( generic_var$Comparison, levels = out$tree_row$labels[out$tree_row$order] )
###
pdf(file="figures/cytof_summary_of_results_bcell_plasmablast_logfc.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= logFC, )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, B Cell (Plasmablast)\nProteins with FDR<0.05') +
  scale_fill_viridis(name='LogFC')
dev.off()
###
pdf(file="figures/cytof_summary_of_results_bcell_plasmablast_fdr.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= -log10(adj.P.Val), )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, B Cell (Plasmablast)\nProteins with FDR<0.05') +
  scale_fill_viridis(name='-Log10(FDR)')
dev.off()
#############################################################################################################################################
### B Cell (Plasmablast)
#############################################################################################################################################
table(lm_cytof_markers_cohort_all$CellSubset %in% "B Cell (Plasmablast)")
ii <- lm_cytof_markers_cohort_all$CellSubset %in% "B Cell (Plasmablast)"
###
generic_var <- dcast(data = lm_cytof_markers_cohort_all[ii,],formula = Comparison~Marker, value.var = "logFC")
generic_var[is.na(generic_var)]<-0
rownames(generic_var) <- generic_var$Comparison
generic_var <- generic_var [ , colnames(generic_var) %in% astrolabe_markers ]
###
out <- pheatmap(generic_var)
dev.off()
###
generic_var <- lm_cytof_markers_cohort_all[ii,]
generic_var <- generic_var [ generic_var$Marker %in% astrolabe_markers , ]
###
generic_var$Marker <- factor( generic_var$Marker, levels = out$tree_col$labels[out$tree_col$order] )
generic_var$Comparison <- factor( generic_var$Comparison, levels = out$tree_row$labels[out$tree_row$order] )
###
pdf(file="figures/cytof_summary_of_results_bcell_plasmablast_logfc_cohort.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= logFC, )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, B Cell (Plasmablast)\nProteins with FDR<0.05') +
  scale_fill_viridis(name='LogFC')
dev.off()
###
pdf(file="figures/cytof_summary_of_results_bcell_plasmablast_fdr_cohort.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= -log10(adj.P.Val), )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, B Cell (Plasmablast)\nProteins with FDR<0.05') +
  scale_fill_viridis(name='-Log10(FDR)')
dev.off()
#############################################################################################################################################
#############################################################################################################################################
### 
#############################################################################################################################################
ii <- channel_intesity_statistics$CellSubset %in% "B Cell (Plasmablast)"
###
generic_var <- dcast(data = channel_intesity_statistics[ii,],formula = SampleId~ChannelName, value.var = "Q0_95")
###
rownames(generic_var) <- generic_var$SampleId
generic_var <- generic_var [ , colnames(generic_var) %in% unique(sort(lm_cytof_markers_all$Marker[lm_cytof_markers_all$CellSubset %in% "B Cell (Plasmablast)"])) ]
###
generic_var$astrolabe_ID <- rownames(generic_var)
generic_var <- (melt(generic_var))
colnames(generic_var)[colnames(generic_var) %in% "variable"] <- "Marker"
generic_var <- merge(generic_var,cimac_e4412_metadata,by="astrolabe_ID")
###
ii <- c(grep("Baseline",generic_var$Condition),grep("Restaging",generic_var$Condition),
        grep("Pre_Day",generic_var$Condition),grep("End_of_",generic_var$Condition))
##
generic_var <- generic_var[ii,]
###
table( names(table(generic_var$condition)) %in% c( "Baseline_1_BV+ipi","Pre_Day_1_Cycle_2_1_BV+ipi","Restaging_1_BV+ipi","End_of_treatment_1_BV+ipi",
                                                   "Baseline_2_BV+nivo","Pre_Day_1_Cycle_2_2_BV+nivo ","Restaging_2_BV+nivo","End_of_treatment_2_BV+nivo ",
                                                   "Baseline_3_BV+nivo+ipi","Pre_Day_1_Cycle_2_3_BV+nivo+ipi","Restaging_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi" ) )

names(table(generic_var$condition))[ ! names(table(generic_var$condition)) %in% c( "Baseline_1_BV+ipi","Pre_Day_1_Cycle_2_1_BV+ipi","Restaging_1_BV+ipi","End_of_treatment_1_BV+ipi",
                                                                                   "Baseline_2_BV+nivo","Pre_Day_1_Cycle_2_2_BV+nivo","Restaging_2_BV+nivo","End_of_treatment_2_BV+nivo",
                                                                                   "Baseline_3_BV+nivo+ipi","Pre_Day_1_Cycle_2_3_BV+nivo+ipi","Restaging_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi" )]

generic_var$condition <- factor(generic_var$condition, levels = c( "Baseline_1_BV+ipi","Pre_Day_1_Cycle_2_1_BV+ipi","Restaging_1_BV+ipi","End_of_treatment_1_BV+ipi",
                                                                   "Baseline_2_BV+nivo","Pre_Day_1_Cycle_2_2_BV+nivo","Restaging_2_BV+nivo","End_of_treatment_2_BV+nivo",
                                                                   "Baseline_3_BV+nivo+ipi","Pre_Day_1_Cycle_2_3_BV+nivo+ipi","Restaging_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi" ))
###
generic_var$condition <- fct_rev(generic_var$condition)
###
my_comparisons <- list( #c("Baseline_1_BV+ipi","Pre_Day_1_Cycle_2_1_BV+ipi"),
  #c("Baseline_1_BV+ipi","Restaging_1_BV+ipi"),
  #c("Baseline_1_BV+ipi","End_of_treatment_1_BV+ipi"),
  #c("Restaging_1_BV+ipi","End_of_treatment_1_BV+ipi"),
  c("End_of_treatment_2_BV+nivo","End_of_treatment_1_BV+ipi"),
  #c("End_of_treatment_2_BV+nivo","End_of_treatment_3_BV+nivo+ipi"),
  c("End_of_treatment_1_BV+ipi","End_of_treatment_3_BV+nivo+ipi"),
  #c("Restaging_1_BV+ipi","Restaging_2_BV+nivo"),
  #c("Restaging_3_BV+nivo+ipi","Restaging_2_BV+nivo"),
  c("Restaging_1_BV+ipi","Restaging_3_BV+nivo+ipi")
  #c("Baseline_2_BV+nivo","Pre_Day_1_Cycle_2_2_BV+nivo"),
  #c("Baseline_2_BV+nivo","Restaging_2_BV+nivo"),
  #c("Baseline_2_BV+nivo","End_of_treatment_2_BV+nivo"),
  #c("Restaging_2_BV+nivo","End_of_treatment_2_BV+nivo"),
  #c("Baseline_3_BV+nivo+ipi","Pre_Day_1_Cycle_2_3_BV+nivo+ipi"),
  #c("Baseline_3_BV+nivo+ipi","Restaging_3_BV+nivo+ipi"),
  #c("Baseline_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi"),
  #c("Restaging_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi")
)
#############################################################################################################################################

pdf(file="figures/examples_b_cell_plasmablast_cytof_markers_diff_expr.pdf",width = 7, height = 4)

ii <- generic_var$Marker %in% "TIGIT" & generic_var$value < 1
test <- generic_var[ii,] ; test <- test[!is.na(test$condition),]

ggboxplot(data=test, 
          x="condition", y="value",fill="Cohort") + theme_minimal() +
  # geom_boxplot(width=0.1) + facet_wrap(~protein,ncol=4) +
  # stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("TIGIT") + labs(x = "", y = "Marker Exp. Values") + theme(legend.position="none")

ii <- generic_var$Marker %in% "4_1BB"  & generic_var$value < 2
test <- generic_var[ii,] ; test <- test[!is.na(test$condition),]

ggboxplot(data=test,
          x="condition", y="value",fill="Cohort") + theme_minimal() +
  #geom_boxplot(width=0.1) + facet_wrap(~protein,ncol=4) +
  #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("4_1BB") + labs(x = "", y = "Marker Exp. Values") + theme(legend.position="none")

ii <- generic_var$Marker %in% "ICOS" & generic_var$value < 1.5
test <- generic_var[ii,] ; test <- test[!is.na(test$condition),]

ggboxplot(data=test,
          x="condition", y="value",fill="Cohort") + theme_minimal() +
  #geom_boxplot(width=0.1) + facet_wrap(~protein,ncol=4) +
  #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("ICOS") + labs(x = "", y = "Marker Exp. Values") + theme(legend.position="none")

ii <- generic_var$Marker %in% "CXCR5" #& generic_var$value < 2
test <- generic_var[ii,] ; test <- test[!is.na(test$condition),]

ggboxplot(data=test,
          x="condition", y="value",fill="Cohort") + theme_minimal() +
  #geom_boxplot(width=0.1) + facet_wrap(~protein,ncol=4) +
  #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("CXCR5") + labs(x = "", y = "Marker Exp. Values") + theme(legend.position="none")

dev.off()

#############################################################################################################################################

my_cells
#############################################################################################################################################
### Granulocyte (Basophil)
#############################################################################################################################################
table(lm_cytof_markers_all$CellSubset %in% "Granulocyte (Basophil)")
ii <- lm_cytof_markers_all$CellSubset %in% "Granulocyte (Basophil)"
###
generic_var <- dcast(data = lm_cytof_markers_all[ii,],formula = Comparison~Marker, value.var = "logFC")
generic_var[is.na(generic_var)]<-0
rownames(generic_var) <- generic_var$Comparison
generic_var <- generic_var [ , colnames(generic_var) %in% astrolabe_markers ]
###
out <- pheatmap(generic_var)
dev.off()
###
generic_var <- lm_cytof_markers_all[ii,]
generic_var <- generic_var [ generic_var$Marker %in% astrolabe_markers , ]
###
generic_var$Marker <- factor( generic_var$Marker, levels = out$tree_col$labels[out$tree_col$order] )
generic_var$Comparison <- factor( generic_var$Comparison, levels = out$tree_row$labels[out$tree_row$order] )
###
pdf(file="figures/cytof_summary_of_results_Basophil_logfc.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= logFC, )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, Granulocyte (Basophil)\nProteins with FDR<0.05') +
  scale_fill_viridis(name='LogFC')
dev.off()
###
pdf(file="figures/cytof_summary_of_results_Basophil_fdr.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= -log10(adj.P.Val), )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, Granulocyte (Basophil)\nProteins with FDR<0.05') +
  scale_fill_viridis(name='-Log10(FDR)')
dev.off()
#############################################################################################################################################
### Granulocyte (Basophil)
#############################################################################################################################################
table(lm_cytof_markers_cohort_all$CellSubset %in% "Granulocyte (Basophil)")
ii <- lm_cytof_markers_cohort_all$CellSubset %in% "Granulocyte (Basophil)"
###
generic_var <- dcast(data = lm_cytof_markers_cohort_all[ii,],formula = Comparison~Marker, value.var = "logFC")
generic_var[is.na(generic_var)]<-0
rownames(generic_var) <- generic_var$Comparison
generic_var <- generic_var [ , colnames(generic_var) %in% astrolabe_markers ]
###
out <- pheatmap(generic_var)
dev.off()
###
generic_var <- lm_cytof_markers_cohort_all[ii,]
generic_var <- generic_var [ generic_var$Marker %in% astrolabe_markers , ]
###
generic_var$Marker <- factor( generic_var$Marker, levels = out$tree_col$labels[out$tree_col$order] )
generic_var$Comparison <- factor( generic_var$Comparison, levels = out$tree_row$labels[out$tree_row$order] )
###
pdf(file="figures/cytof_summary_of_results_Basophil_logfc_cohort.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= logFC, )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, Granulocyte (Basophil)\nProteins with FDR<0.05') +
  scale_fill_viridis(name='LogFC')
dev.off()
###
pdf(file="figures/cytof_summary_of_results_Basophil_fdr_cohort.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= -log10(adj.P.Val), )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, Granulocyte (Basophil)\nProteins with FDR<0.05') +
  scale_fill_viridis(name='-Log10(FDR)')
dev.off()
#############################################################################################################################################


my_cells
#############################################################################################################################################
### "NKT Cell"
#############################################################################################################################################
table(lm_cytof_markers_all$CellSubset %in% "NKT Cell")
ii <- lm_cytof_markers_all$CellSubset %in% "NKT Cell"
###
generic_var <- dcast(data = lm_cytof_markers_all[ii,],formula = Comparison~Marker, value.var = "logFC")
generic_var[is.na(generic_var)]<-0
rownames(generic_var) <- generic_var$Comparison
generic_var <- generic_var [ , colnames(generic_var) %in% astrolabe_markers ]
###
#out <- pheatmap(generic_var)
#dev.off()
###
generic_var <- lm_cytof_markers_all[ii,]
generic_var <- generic_var [ generic_var$Marker %in% astrolabe_markers , ]
###
#generic_var$Marker <- factor( generic_var$Marker, levels = out$tree_col$labels[out$tree_col$order] )
#generic_var$Comparison <- factor( generic_var$Comparison, levels = out$tree_row$labels[out$tree_row$order] )
###
pdf(file="figures/cytof_summary_of_results_nkt_cell_logfc.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= logFC, )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, "NKT Cell"\nProteins with FDR<0.05') +
  scale_fill_viridis(name='LogFC')
dev.off()
###
pdf(file="figures/cytof_summary_of_results_nkt_cell_fdr.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= -log10(adj.P.Val), )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, "NKT Cell"\nProteins with FDR<0.05') +
  scale_fill_viridis(name='-Log10(FDR)')
dev.off()
#############################################################################################################################################
### "NKT Cell"
#############################################################################################################################################
table(lm_cytof_markers_cohort_all$CellSubset %in% "NKT Cell")
ii <- lm_cytof_markers_cohort_all$CellSubset %in% "NKT Cell"
###
generic_var <- dcast(data = lm_cytof_markers_cohort_all[ii,],formula = Comparison~Marker, value.var = "logFC")
generic_var[is.na(generic_var)]<-0
rownames(generic_var) <- generic_var$Comparison
generic_var <- generic_var [ , colnames(generic_var) %in% astrolabe_markers ]
###
out <- pheatmap(generic_var)
dev.off()
###
generic_var <- lm_cytof_markers_cohort_all[ii,]
generic_var <- generic_var [ generic_var$Marker %in% astrolabe_markers , ]
###
generic_var$Marker <- factor( generic_var$Marker, levels = out$tree_col$labels[out$tree_col$order] )
generic_var$Comparison <- factor( generic_var$Comparison, levels = out$tree_row$labels[out$tree_row$order] )
###
pdf(file="figures/cytof_summary_of_results_nkt_cell_logfc_cohort.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= logFC, )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, "NKT Cell"\nProteins with FDR<0.05') +
  scale_fill_viridis(name='LogFC')
dev.off()
###
pdf(file="figures/cytof_summary_of_results_nkt_cell_fdr_cohort.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= -log10(adj.P.Val), )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, "NKT Cell"\nProteins with FDR<0.05') +
  scale_fill_viridis(name='-Log10(FDR)')
dev.off()
#############################################################################################################################################




my_cells
#############################################################################################################################################
### "CD4+ T Cell (Treg)"
#############################################################################################################################################
table(lm_cytof_markers_all$CellSubset %in% "CD4+ T Cell (Treg)")
ii <- lm_cytof_markers_all$CellSubset %in% "CD4+ T Cell (Treg)"
###
generic_var <- dcast(data = lm_cytof_markers_all[ii,],formula = Comparison~Marker, value.var = "logFC")
generic_var[is.na(generic_var)]<-0
rownames(generic_var) <- generic_var$Comparison
generic_var <- generic_var [ , colnames(generic_var) %in% astrolabe_markers ]
###
out <- pheatmap(generic_var)
dev.off()
###
generic_var <- lm_cytof_markers_all[ii,]
generic_var <- generic_var [ generic_var$Marker %in% astrolabe_markers , ]
###
generic_var$Marker <- factor( generic_var$Marker, levels = out$tree_col$labels[out$tree_col$order] )
generic_var$Comparison <- factor( generic_var$Comparison, levels = out$tree_row$labels[out$tree_row$order] )
###
pdf(file="figures/cytof_summary_of_results_cd4p_t_reg_logfc.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= logFC, )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, "CD4+ T Cell (Treg)"\nProteins with FDR<0.05') +
  scale_fill_viridis(name='LogFC')
dev.off()
###
pdf(file="figures/cytof_summary_of_results_cd4p_t_reg_fdr.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= -log10(adj.P.Val), )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, "CD4+ T Cell (Treg)"\nProteins with FDR<0.05') +
  scale_fill_viridis(name='-Log10(FDR)')
dev.off()
#############################################################################################################################################
### "CD4+ T Cell (Treg)"
#############################################################################################################################################
table(lm_cytof_markers_cohort_all$CellSubset %in% "CD4+ T Cell (Treg)")
ii <- lm_cytof_markers_cohort_all$CellSubset %in% "CD4+ T Cell (Treg)"
###
generic_var <- dcast(data = lm_cytof_markers_cohort_all[ii,],formula = Comparison~Marker, value.var = "logFC")
generic_var[is.na(generic_var)]<-0
rownames(generic_var) <- generic_var$Comparison
generic_var <- generic_var [ , colnames(generic_var) %in% astrolabe_markers ]
###
out <- pheatmap(generic_var)
dev.off()
###
generic_var <- lm_cytof_markers_cohort_all[ii,]
generic_var <- generic_var [ generic_var$Marker %in% astrolabe_markers , ]
###
generic_var$Marker <- factor( generic_var$Marker, levels = out$tree_col$labels[out$tree_col$order] )
generic_var$Comparison <- factor( generic_var$Comparison, levels = out$tree_row$labels[out$tree_row$order] )
###
pdf(file="figures/cytof_summary_of_results_cd4p_t_reg_logfc_cohort.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= logFC, )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, "CD4+ T Cell (Treg)"\nProteins with FDR<0.05') +
  scale_fill_viridis(name='LogFC')
dev.off()
###
pdf(file="figures/cytof_summary_of_results_cd4p_t_reg_fdr_cohort.pdf",width = 5,height = 3)
ggplot(generic_var, aes(Comparison, Marker, fill= -log10(adj.P.Val), )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, "CD4+ T Cell (Treg)"\nProteins with FDR<0.05') +
  scale_fill_viridis(name='-Log10(FDR)')
dev.off()
#############################################################################################################################################





#############################################################################################################################################
### Monocyte (CD14- CD16+)
#############################################################################################################################################
table(lm_cytof_markers_all$CellSubset %in% "Monocyte (CD14- CD16+)")
ii <- lm_cytof_markers_all$CellSubset %in% "Monocyte (CD14- CD16+)"
###
generic_var <- dcast(data = lm_cytof_markers_all[ii,],formula = Comparison~Marker, value.var = "logFC")
generic_var[is.na(generic_var)]<-0
rownames(generic_var) <- generic_var$Comparison
generic_var <- generic_var [ , colnames(generic_var) %in% astrolabe_markers ]
###
out <- pheatmap(generic_var)
dev.off()
###
generic_var <- lm_cytof_markers_all[ii,]
generic_var <- generic_var [ generic_var$Marker %in% astrolabe_markers , ]
###
generic_var$Marker <- factor( generic_var$Marker, levels = out$tree_col$labels[out$tree_col$order] )
generic_var$Comparison <- factor( generic_var$Comparison, levels = out$tree_row$labels[out$tree_row$order] )
###
pdf(file="figures/cytof_summary_of_results_monocytes_cd14n_cd16p_logfc.pdf",width = 8,height = 6)
ggplot(generic_var, aes(Comparison, Marker, fill= logFC, )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, Monocyte (CD14- CD16+), Proteins with FDR<0.05') +
  scale_fill_viridis(name='LogFC')
dev.off()
###
pdf(file="figures/cytof_summary_of_results_monocytes_cd14n_cd16p_fdr.pdf",width = 8,height = 6)
ggplot(generic_var, aes(Comparison, Marker, fill= -log10(adj.P.Val), )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, Monocyte (CD14- CD16+), Proteins with FDR<0.05') +
  scale_fill_viridis(name='-Log10(FDR)')
dev.off()
#############################################################################################################################################




#############################################################################################################################################
### Monocyte (CD14+ CD16-)
#############################################################################################################################################
table(lm_cytof_markers_all$CellSubset %in% "Monocyte (CD14+ CD16-)")
ii <- lm_cytof_markers_all$CellSubset %in% "Monocyte (CD14+ CD16-)"
###
generic_var <- lm_cytof_markers_all[ii,]
generic_var <- generic_var [ generic_var$Marker %in% astrolabe_markers , ]
###
pdf(file="figures/cytof_summary_of_results_monocytes_cd14p_cd16n_logfc.pdf",width = 8,height = 6)
ggplot(generic_var, aes(Comparison, Marker, fill= logFC, )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, Monocyte (CD14+ CD16-), Proteins with FDR<0.05') +
  scale_fill_viridis(name='LogFC')
dev.off()
###
pdf(file="figures/cytof_summary_of_results_monocytes_cd14p_cd16n_fdr.pdf",width = 8,height = 6)
ggplot(generic_var, aes(Comparison, Marker, fill= -log10(adj.P.Val), )) + 
  geom_tile(color='white') + theme_bw() +
  rotate_x_text(angle=45) +
  labs(x ='', y='', title='Cytof, Monocyte (CD14+ CD16-), Proteins with FDR<0.05') +
  scale_fill_viridis(name='-Log10(FDR)')
dev.off()
#############################################################################################################################################
#############################################################################################################################################
### 
#############################################################################################################################################
# head(channel_intesity_statistics)
# summary(channel_intesity_statistics$Max)
ii <- channel_intesity_statistics$CellSubset %in% "Monocyte (CD14+ CD16+)"
# ii <- channel_intesity_statistics$CellSubset %in% "CD8+ T Cell (Effector Memory)"

###
generic_var <- dcast(data = channel_intesity_statistics[ii,],formula = SampleId~ChannelName, value.var = "Q0_95")
###
# dim(generic_var)
rownames(generic_var) <- generic_var$SampleId
generic_var <- generic_var [ , colnames(generic_var) %in% unique(sort(lm_cytof_markers_all$Marker[lm_cytof_markers_all$CellSubset %in% "CD8+ T Cell (Effector Memory)"])) ]
###
generic_var$astrolabe_ID <- rownames(generic_var)
generic_var <- (melt(generic_var))
colnames(generic_var)[colnames(generic_var) %in% "variable"] <- "Marker"
generic_var <- merge(generic_var,cimac_e4412_metadata,by="astrolabe_ID")

ii <- c(grep("Baseline",generic_var$Condition),grep("Restaging",generic_var$Condition),
        grep("Pre_Day",generic_var$Condition),grep("End_of_",generic_var$Condition))

generic_var <- generic_var[ii,]

table( names(table(generic_var$condition)) %in% c( "Baseline_1_BV+ipi","Pre_Day_1_Cycle_2_1_BV+ipi","Restaging_1_BV+ipi","End_of_treatment_1_BV+ipi",
                                                   "Baseline_2_BV+nivo","Pre_Day_1_Cycle_2_2_BV+nivo ","Restaging_2_BV+nivo","End_of_treatment_2_BV+nivo ",
                                                   "Baseline_3_BV+nivo+ipi","Pre_Day_1_Cycle_2_3_BV+nivo+ipi","Restaging_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi" ) )

names(table(generic_var$condition))[ ! names(table(generic_var$condition)) %in% c( "Baseline_1_BV+ipi","Pre_Day_1_Cycle_2_1_BV+ipi","Restaging_1_BV+ipi","End_of_treatment_1_BV+ipi",
                                                                                   "Baseline_2_BV+nivo","Pre_Day_1_Cycle_2_2_BV+nivo","Restaging_2_BV+nivo","End_of_treatment_2_BV+nivo",
                                                                                   "Baseline_3_BV+nivo+ipi","Pre_Day_1_Cycle_2_3_BV+nivo+ipi","Restaging_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi" )]

generic_var$condition <- factor(generic_var$condition, levels = c( "Baseline_1_BV+ipi","Pre_Day_1_Cycle_2_1_BV+ipi","Restaging_1_BV+ipi","End_of_treatment_1_BV+ipi",
                                                                   "Baseline_2_BV+nivo","Pre_Day_1_Cycle_2_2_BV+nivo","Restaging_2_BV+nivo","End_of_treatment_2_BV+nivo",
                                                                   "Baseline_3_BV+nivo+ipi","Pre_Day_1_Cycle_2_3_BV+nivo+ipi","Restaging_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi" ))

generic_var$condition <- fct_rev(generic_var$condition)

my_comparisons <- list( #c("Baseline_1_BV+ipi","Pre_Day_1_Cycle_2_1_BV+ipi"),
  #c("Baseline_1_BV+ipi","Restaging_1_BV+ipi"),
  #c("Baseline_1_BV+ipi","End_of_treatment_1_BV+ipi"),
  #c("Restaging_1_BV+ipi","End_of_treatment_1_BV+ipi"),
  c("End_of_treatment_2_BV+nivo","End_of_treatment_1_BV+ipi"),
  #c("End_of_treatment_2_BV+nivo","End_of_treatment_3_BV+nivo+ipi"),
  c("End_of_treatment_1_BV+ipi","End_of_treatment_3_BV+nivo+ipi"),
  #c("Restaging_1_BV+ipi","Restaging_2_BV+nivo"),
  #c("Restaging_3_BV+nivo+ipi","Restaging_2_BV+nivo"),
  c("Restaging_1_BV+ipi","Restaging_3_BV+nivo+ipi")
  #c("Baseline_2_BV+nivo","Pre_Day_1_Cycle_2_2_BV+nivo"),
  #c("Baseline_2_BV+nivo","Restaging_2_BV+nivo"),
  #c("Baseline_2_BV+nivo","End_of_treatment_2_BV+nivo"),
  #c("Restaging_2_BV+nivo","End_of_treatment_2_BV+nivo"),
  
  #c("Baseline_3_BV+nivo+ipi","Pre_Day_1_Cycle_2_3_BV+nivo+ipi"),
  #c("Baseline_3_BV+nivo+ipi","Restaging_3_BV+nivo+ipi"),
  #c("Baseline_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi"),
  #c("Restaging_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi")
)
#############################################################################################################################################



#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
###
### cell frequency abundances figure
###
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

generic_var <- lm_frequencies_all[lm_frequencies_all$adj.P.Val < 0.15,]
# levels(as.factor(generic_var$Comparison))
generic_var$Comparison <- as.character(generic_var$Comparison)
###
generic_var <- dcast(data = generic_var,formula = Comparison~Celltype, value.var = "logFC")
generic_var[is.na(generic_var)]<-0
rownames(generic_var) <- generic_var$Comparison
###
# ord <- hclust( dist(scale((generic_var[,-1])), method = "euclidean"), method = "ward.D" )$order
# ord
out <- pheatmap(generic_var[,-1])
dev.off()
###
generic_var <- lm_frequencies_all[lm_frequencies_all$adj.P.Val<0.15,]
generic_var$Comparison <- as.character(generic_var$Comparison)
###
generic_var$Celltype <- factor( generic_var$Celltype, levels = out$tree_col$labels[out$tree_col$order] )
generic_var$Comparison <- factor( generic_var$Comparison, levels = out$tree_row$labels[out$tree_row$order] )
generic_var$nLogFDR <- -log10(generic_var$adj.P.Val)
###
ix <- generic_var$Comparison %in% c("1A_vs_1B","3A_vs_3B","1A_vs_2A","1D_vs_3D","1A_vs_3A","1B_vs_1D","1A_vs_1C","1A_vs_1D","3B_vs_3D","3A_vs_3D","2A_vs_2D","2B_vs_2D")
###
iy <- generic_var$Comparison %in% c("Baseline_vs_1B","Baseline_vs_1C","Baseline_vs_3B","Baseline_vs_1D","Baseline_vs_3D","Baseline_vs_2D","2B_vs_3B", "1C_vs_3C","2D_vs_3D","1C_vs_2C","1D_vs_2D")
###
pdf(file="figures/cytof_cell_frequencies_summary_of_results_within.pdf",width = 6.5,height = 4.5)
ggplot(generic_var[ix,], aes(Comparison, Celltype, color= logFC, size=nLogFDR )) + 
  #geom_tile(color='white') + 
  geom_point(shape=15)+
  theme_bw() +
  rotate_x_text(angle=45) +
  # scale_fill_viridis(name='LogFC', limits = c(min(generic_var$value), max(generic_var$value)) ) +
  # scale_fill_gradientn(limits = c(min(generic_var$value), max(generic_var$value))) +
  scale_color_gradient2(low='firebrick', mid="white",high='steelblue', name='LogFC') +
  labs(x ='', y='', title='DA Celltypes') #+
  # theme(legend.position = "bottom") 
dev.off()
###
###
pdf(file="figures/cytof_cell_frequencies_summary_of_results_treatment.pdf",width = 6.5,height = 4.5)
ggplot(generic_var[iy,], aes(Comparison, Celltype, color= logFC, size=nLogFDR )) + 
  #geom_tile(color='white') + 
  geom_point(shape=15)+
  theme_bw() +
  rotate_x_text(angle=45) +
  # scale_fill_viridis(name='LogFC', limits = c(min(generic_var$value), max(generic_var$value)) ) +
  # scale_fill_gradientn(limits = c(min(generic_var$value), max(generic_var$value))) +
  scale_color_gradient2(low='firebrick', mid="white",high='steelblue', name='LogFC') +
  labs(x ='', y='', title='DA Celltypes') #+
# theme(legend.position = "bottom") 
dev.off()
###
###
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
###
### DAA example figures!
###
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#table(generic_var$Celltype,generic_var$Comparison)
### 
generic_var <- melt(v_abundance$E)
### 
colnames(generic_var) <- c("Celltype","astrolabe_ID","Log2Cells")
generic_var <- merge(generic_var,cimac_e4412_metadata,by="astrolabe_ID")
### 
generic_var <- generic_var[ ! generic_var$condition %in% c("Other_1_BV+ipi","Off_study_1_BV+ipi","Off_study_2_BV+nivo") , ]
###
levels(as.factor(generic_var$condition)) 
generic_var$condition <- as.factor(generic_var$condition)
###
generic_var$condition <- factor(generic_var$condition, levels = c("Baseline_1_BV+ipi","Pre_Day_1_Cycle_2_1_BV+ipi","Restaging_1_BV+ipi","End_of_treatment_1_BV+ipi",
                                                                  "Baseline_2_BV+nivo","Pre_Day_1_Cycle_2_2_BV+nivo","Restaging_2_BV+nivo","End_of_treatment_2_BV+nivo",
                                                                  "Baseline_3_BV+nivo+ipi","Restaging_3_BV+nivo+ipi","Pre_Day_1_Cycle_2_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi"))
###
#############################################################################################################################################
my_comparisons <- list( c("Baseline_1_BV+ipi","Pre_Day_1_Cycle_2_1_BV+ipi"),
                        c("Baseline_1_BV+ipi","Restaging_1_BV+ipi"),
                        c("Baseline_1_BV+ipi","End_of_treatment_1_BV+ipi"),
                        c("Restaging_1_BV+ipi","End_of_treatment_1_BV+ipi"),
                        c("Pre_Day_1_Cycle_2_1_BV+ipi","End_of_treatment_1_BV+ipi"),
                        
                        c("Baseline_2_BV+nivo","Pre_Day_1_Cycle_2_2_BV+nivo"),
                        c("Baseline_2_BV+nivo","Restaging_2_BV+nivo"),
                        c("Baseline_2_BV+nivo","End_of_treatment_2_BV+nivo"),
                        c("Restaging_2_BV+nivo","End_of_treatment_2_BV+nivo"),
                        c("Pre_Day_1_Cycle_2_2_BV+nivo","End_of_treatment_2_BV+nivo"),
                        
                        c("Baseline_3_BV+nivo+ipi","Pre_Day_1_Cycle_2_3_BV+nivo+ipi"),
                        c("Baseline_3_BV+nivo+ipi","Restaging_3_BV+nivo+ipi"),
                        c("Baseline_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi"),
                        c("Restaging_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi"),
                        c("Pre_Day_1_Cycle_2_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi")
)
#############################################################################################################################################
table( levels(as.factor(generic_var$condition)) %in% unlist(my_comparisons) )
levels(generic_var$Condition)[! levels(generic_var$Condition) %in% unlist(my_comparisons)]
#############################################################################################################################################

generic_var$condition <- fct_rev(generic_var$condition)

pdf(file="figures/examples_DA_celltype_diff_abu.pdf",width = 7, height = 4.5)

ii <- generic_var$Celltype %in% "CD4+ T Cell (Effector Memory)"

ggboxplot(data=generic_var[ii,], 
          x="condition", y="Log2Cells",fill="Cohort") + theme_minimal() +
  facet_wrap(~Cohort,nrow = 3,scales = 'free') +
  geom_quasirandom(size=1,color="grey50")+
  stat_compare_means(comparisons = my_comparisons, 
                     method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 0,hjust = 1)) +
  coord_flip()+ 
  ggtitle("CD4+ T Cell (Effector Memory)") + labs(x = "", y = "Log2Cells") + theme(legend.position="none")

ii <- generic_var$Celltype %in% "CD4+ T Cell (Central Memory)"

ggboxplot(data=generic_var[ii,], 
          x="condition", y="Log2Cells",fill="Cohort") + theme_minimal() +
  facet_wrap(~Cohort,nrow = 3,scales = 'free') +
  geom_quasirandom(size=1,color="grey50")+
  #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 0,hjust = 1)) +
  coord_flip()+ 
  ggtitle("CD4+ T Cell (Central Memory)") + labs(x = "", y = "Log2Cells") + theme(legend.position="none")

ii <- generic_var$Celltype %in% "CD8+ T Cell (Effector Memory)"

ggboxplot(data=generic_var[ii,], 
          x="condition", y="Log2Cells",fill="Cohort") + theme_minimal() +
  facet_wrap(~Cohort,nrow = 3,scales = 'free') +
  geom_quasirandom(size=1,color="grey50")+
  geom_quasirandom(size=1,color="grey50")+
  #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("CD8+ T Cell (Effector Memory)") + labs(x = "", y = "Log2Cells") + theme(legend.position="none")

ii <- generic_var$Celltype %in% "CD8+ T Cell (Central Memory)"

ggboxplot(data=generic_var[ii,], 
          x="condition", y="Log2Cells",fill="Cohort") + theme_minimal() +
  facet_wrap(~Cohort,nrow = 3,scales = 'free') +
  geom_quasirandom(size=1,color="grey50")+
  #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("CD8+ T Cell (Central Memory)") + labs(x = "", y = "Log2Cells") + theme(legend.position="none")

ii <- generic_var$Celltype %in% "CD4+ T Cell (EMRA)"

ggboxplot(data=generic_var[ii,], 
          x="condition", y="Log2Cells",fill="Cohort") + theme_minimal() +
  facet_wrap(~Cohort,nrow = 3,scales = 'free') +
  geom_quasirandom(size=1,color="grey50")+
  #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("CD4+ T Cell (EMRA)") + labs(x = "", y = "Log2Cells") + theme(legend.position="none")

ii <- generic_var$Celltype %in% "B Cell (Plasmablast)" 

ggboxplot(data=generic_var[ii,], 
          x="condition", y="Log2Cells",fill="Cohort") + theme_minimal() +
  facet_wrap(~Cohort,nrow = 3,scales = 'free') +
  geom_quasirandom(size=1,color="grey50")+
  #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("B Cell (Plasmablast)" ) + labs(x = "", y = "Log2Cells") + theme(legend.position="none")

ii <- generic_var$Celltype %in% "CM-_unassigned" 

ggboxplot(data=generic_var[ii,], 
          x="condition", y="Log2Cells",fill="Cohort") + theme_minimal() +
  facet_wrap(~Cohort,nrow = 3,scales = 'free') +
  geom_quasirandom(size=1,color="grey50")+
  #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("CM-_unassigned" ) + labs(x = "", y = "Log2Cells") + theme(legend.position="none")

ii <- generic_var$Celltype %in% "Monocyte (CD14+ CD16-)" 

ggboxplot(data=generic_var[ii,], 
          x="condition", y="Log2Cells",fill="Cohort") + theme_minimal() +
  facet_wrap(~Cohort,nrow = 3,scales = 'free') +
  geom_quasirandom(size=1,color="grey50")+
  #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("Monocyte (CD14+ CD16-)" ) + labs(x = "", y = "Log2Cells") + theme(legend.position="none")

ii <- generic_var$Celltype %in% "CD4+ CD8+ T Cell" 

ggboxplot(data=generic_var[ii,], 
          x="condition", y="Log2Cells",fill="Cohort") + theme_minimal() +
  facet_wrap(~Cohort,nrow = 3,scales = 'free') +
  geom_quasirandom(size=1,color="grey50")+
  #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("CD4+ CD8+ T Cell" ) + labs(x = "", y = "Log2Cells") + theme(legend.position="none")

dev.off()

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
###
### CYTOF MARKERS  FIGUREs
###
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

#############################################################################################################################################
### 
#############################################################################################################################################
ii <- channel_intesity_statistics$CellSubset %in% c("CD4+ T Cell (Effector Memory)")
###
generic_var <- dcast(data = channel_intesity_statistics[ii,],formula = SampleId~ChannelName, value.var = "Q0_95")
rownames(generic_var) <- generic_var$SampleId
generic_var <- generic_var [ , colnames(generic_var) %in% unique(sort(lm_cytof_markers_all$Marker[lm_cytof_markers_all$adj.P.Val<0.05])) ]
###
generic_var$astrolabe_ID <- rownames(generic_var)
generic_var <- melt(generic_var)
colnames(generic_var)[colnames(generic_var) %in% "variable"] <- "Marker"
generic_var <- merge(generic_var,cimac_e4412_metadata,by="astrolabe_ID")

ii <- c(grep("Baseline",generic_var$Condition),grep("Restaging",generic_var$Condition),
        grep("Pre_Day",generic_var$Condition),grep("End_of_",generic_var$Condition))

generic_var <- generic_var[ii,]

table( names(table(generic_var$condition)) %in% c( "Baseline_1_BV+ipi","Pre_Day_1_Cycle_2_1_BV+ipi","Restaging_1_BV+ipi","End_of_treatment_1_BV+ipi",
                                                   "Baseline_2_BV+nivo","Pre_Day_1_Cycle_2_2_BV+nivo ","Restaging_2_BV+nivo","End_of_treatment_2_BV+nivo ",
                                                   "Baseline_3_BV+nivo+ipi","Pre_Day_1_Cycle_2_3_BV+nivo+ipi","Restaging_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi" ) )
names(table(generic_var$condition))[ ! names(table(generic_var$condition)) %in% c( "Baseline_1_BV+ipi","Pre_Day_1_Cycle_2_1_BV+ipi","Restaging_1_BV+ipi","End_of_treatment_1_BV+ipi",
                                                                                   "Baseline_2_BV+nivo","Pre_Day_1_Cycle_2_2_BV+nivo ","Restaging_2_BV+nivo","End_of_treatment_2_BV+nivo ",
                                                                                   "Baseline_3_BV+nivo+ipi","Pre_Day_1_Cycle_2_3_BV+nivo+ipi","Restaging_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi" ) ]

generic_var$condition <- factor(generic_var$condition, levels = c( "Baseline_1_BV+ipi","Pre_Day_1_Cycle_2_1_BV+ipi","Restaging_1_BV+ipi","End_of_treatment_1_BV+ipi",
                                                                   "Baseline_2_BV+nivo","Pre_Day_1_Cycle_2_2_BV+nivo","Restaging_2_BV+nivo","End_of_treatment_2_BV+nivo",
                                                                   "Baseline_3_BV+nivo+ipi","Pre_Day_1_Cycle_2_3_BV+nivo+ipi","Restaging_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi" ))

generic_var$condition <- fct_rev(generic_var$condition)

my_comparisons <- list( 
  c("Baseline_1_BV+ipi","Pre_Day_1_Cycle_2_1_BV+ipi"),
  c("Baseline_1_BV+ipi","Restaging_1_BV+ipi"),
  c("Baseline_1_BV+ipi","End_of_treatment_1_BV+ipi"),
  c("Restaging_1_BV+ipi","End_of_treatment_1_BV+ipi"),
  c("End_of_treatment_2_BV+nivo","End_of_treatment_1_BV+ipi"),
  c("End_of_treatment_2_BV+nivo","End_of_treatment_3_BV+nivo+ipi"),
  c("End_of_treatment_1_BV+ipi","End_of_treatment_3_BV+nivo+ipi"),
  c("Restaging_1_BV+ipi","Restaging_2_BV+nivo"),
  c("Restaging_3_BV+nivo+ipi","Restaging_2_BV+nivo"),
  c("Restaging_1_BV+ipi","Restaging_3_BV+nivo+ipi"),
  c("Baseline_2_BV+nivo","Pre_Day_1_Cycle_2_2_BV+nivo"),
  c("Baseline_2_BV+nivo","Restaging_2_BV+nivo"),
  c("Baseline_2_BV+nivo","End_of_treatment_2_BV+nivo"),
  c("Restaging_2_BV+nivo","End_of_treatment_2_BV+nivo"),
  c("Baseline_3_BV+nivo+ipi","Pre_Day_1_Cycle_2_3_BV+nivo+ipi"),
  c("Baseline_3_BV+nivo+ipi","Restaging_3_BV+nivo+ipi"),
  c("Baseline_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi"),
  c("Restaging_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi")
)
#############################################################################################################################################

#############################################################################################################################################
examples_cytof_markers_diff_expr <- list()
#############################################################################################################################################

#############################################################################################################################################
ii <- generic_var$Marker %in% "ICOS"
examples_cytof_markers_diff_expr[[1]] <- ggboxplot(data=generic_var[ii,], 
          x="condition", y="value",fill="Cohort") + theme_minimal() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("ICOS - CD4+ T Cell (Effector Memory)") + labs(x = "", y = "Marker Exp. Values") + theme(legend.position="none")

ii <- generic_var$Marker %in% "PD_1"
examples_cytof_markers_diff_expr[[2]] <- ggboxplot(data=generic_var[ii,], 
          x="condition", y="value",fill="Cohort") + theme_minimal() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("PD1 - CD4+ T Cell (Effector Memory)") + labs(x = "", y = "Marker Exp. Values") + theme(legend.position="none")

ii <- generic_var$Marker %in% "CD39"
examples_cytof_markers_diff_expr[[3]] <- ggboxplot(data=generic_var[ii,], 
                                                   x="condition", y="value",fill="Cohort") + theme_minimal() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("CD39 - CD4+ T Cell (Effector Memory)") + labs(x = "", y = "Marker Exp. Values") + theme(legend.position="none")
#############################################################################################################################################

#############################################################################################################################################
ii <- channel_intesity_statistics$CellSubset %in% c("CD8+ T Cell (Effector Memory)")
generic_var <- dcast(data = channel_intesity_statistics[ii,],formula = SampleId~ChannelName, value.var = "Q0_95")
rownames(generic_var) <- generic_var$SampleId
generic_var <- generic_var [ , colnames(generic_var) %in% unique(sort(lm_cytof_markers_all$Marker[lm_cytof_markers_all$adj.P.Val<0.05])) ]
generic_var$astrolabe_ID <- rownames(generic_var)
generic_var <- melt(generic_var)
colnames(generic_var)[colnames(generic_var) %in% "variable"] <- "Marker"
generic_var <- merge(generic_var,cimac_e4412_metadata,by="astrolabe_ID")
ii <- c(grep("Baseline",generic_var$Condition),grep("Restaging",generic_var$Condition),
        grep("Pre_Day",generic_var$Condition),grep("End_of_",generic_var$Condition))
generic_var <- generic_var[ii,]
###
generic_var$condition <- factor(generic_var$condition, levels = c( "Baseline_1_BV+ipi","Pre_Day_1_Cycle_2_1_BV+ipi","Restaging_1_BV+ipi","End_of_treatment_1_BV+ipi",
                                                                   "Baseline_2_BV+nivo","Pre_Day_1_Cycle_2_2_BV+nivo","Restaging_2_BV+nivo","End_of_treatment_2_BV+nivo",
                                                                   "Baseline_3_BV+nivo+ipi","Pre_Day_1_Cycle_2_3_BV+nivo+ipi","Restaging_3_BV+nivo+ipi","End_of_treatment_3_BV+nivo+ipi" ))

generic_var$condition <- fct_rev(generic_var$condition)
#############################################################################################################################################

ii <- generic_var$Marker %in% "ICOS"
examples_cytof_markers_diff_expr[[4]] <- ggboxplot(data=generic_var[ii,], 
                                                   x="condition", y="value",fill="Cohort") + theme_minimal() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("ICOS - CD8+ T Cell (Effector Memory)") + labs(x = "", y = "Marker Exp. Values") + theme(legend.position="none")

ii <- generic_var$Marker %in% "PD_1"
examples_cytof_markers_diff_expr[[5]] <- ggboxplot(data=generic_var[ii,], 
                                                   x="condition", y="value",fill="Cohort") + theme_minimal() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("PD1 - CD8+ T Cell (Effector Memory)") + labs(x = "", y = "Marker Exp. Values") + theme(legend.position="none")

ii <- generic_var$Marker %in% "CD39"
examples_cytof_markers_diff_expr[[6]] <- ggboxplot(data=generic_var[ii,], 
                                                   x="condition", y="value",fill="Cohort") + theme_minimal() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1)) +
  coord_flip()+ 
  ggtitle("CD39 - CD8+ T Cell (Effector Memory)") + labs(x = "", y = "Marker Exp. Values") + theme(legend.position="none")

#############################################################################################################################################
#############################################################################################################################################

pdf(file="figures/examples_cytof_markers_diff_expr.pdf",width = 7, height = 4)
print(examples_cytof_markers_diff_expr)
dev.off()

#############################################################################################################################################
#############################################################################################################################################


#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################


#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

### saveRDS(file="/Users/gonzae34/Downloads/cytof.RDS",channel_intesity_statistics)
### saveRDS(file="/Users/gonzae34/Downloads/cytof_meta.RDS",cimac_e4412_metadata)

#############################################################################################################################################
### write.csv(file="cimac_e4412_metadata.csv",cimac_e4412_metadata) Latest Organized Metadata
#############################################################################################################################################


#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################


#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
###
### SEROLOGY ANALYSIS
###
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
### "Summary Grand Serology 100.xlsx"
###
grand_serology_od <- as.data.frame(readxl::read_excel("ELISA CIMAC E4412/Summary Grand Serology 100.xlsx", sheet = 1))
###
grand_serology_titers <- as.data.frame(readxl::read_excel("ELISA CIMAC E4412/Summary Grand Serology 100.xlsx", sheet = 3))
###
dim(grand_serology_titers)
###
which(colnames(grand_serology_titers) %in% "Original+Forecast")
###
ncol(grand_serology_titers)
###
serology_OF <- grand_serology_titers[,c(96:128)]
###
serology_OF_meta <- serology_OF[,1:9]
###
serology_OF <- serology_OF[,-c(1:9)]
###
colnames(serology_OF) <- c("NYESO1 (N)","NYESO1 (ID)","P53","MAGE-A1","MAGE-A3 (B)","MAGE-A4","MAGE-A10","SOX2","SSX2","SSX4","CT10","CT47",
                           "MELAN-A","HORMAD1","SURVIVIN-Δ","SURVIVIN 2b","UBTD2","XAGE","XAGE-1b","WT1","PRAME","ERG","GAGE7","DHFR")
###
serology_OF_meta_raw <- grand_serology_titers[,1:9]
###
dim(serology_OF)
dim(serology_OF_meta_raw)
###
serology_OF_meta_raw <- serology_OF_meta_raw[1:167,]
serology_OF <- serology_OF[1:167,]
###
rownames(serology_OF) <- serology_OF_meta_raw$`Processed sample id`
my_labels <- serology_OF_meta_raw$`Processed sample id`
###
###
table(cimac_e4412_metadata$Arm,cimac_e4412_metadata$Cohort)
###
###
serology_OF_meta_raw <- serology_OF_meta_raw[,c(1,2,3)]
colnames(serology_OF_meta_raw) <- c("Arm","ParticipantID","Collection")
test <- cimac_e4412_metadata [ , colnames(cimac_e4412_metadata) %in% c("Arm","Cohort") ]
serology_OF_meta_raw <- merge(serology_OF_meta_raw,unique(test),by="Arm")
###
rownames(serology_OF_meta_raw) <- my_labels
###

###
test <- as.matrix(serology_OF)
test[is.na(test)] <- 0
test <- log10(test+1)
###
dim(serology_OF)
dim(serology_OF_meta_raw)
dim(test)
###

###
pdf(file="figures/serology_heatmap_all_cohorts.pdf",width = 8, height = 6)
pheatmap(test,scale = 'none',show_rownames = FALSE,annotation_row =  serology_OF_meta_raw[,-2], angle_col = 90, color = viridis(255) )
dev.off()

###
pdf(file="figures/serology_heatmap_cohort_1.pdf",width = 8, height = 6)
ii <- serology_OF_meta_raw$Arm %in% c("Arm_A","Arm_B","Arm_C")
pheatmap(test[ii,],scale = 'none',show_rownames = FALSE,annotation_row =  serology_OF_meta_raw[,-2], angle_col = 90, color = viridis(255) )
dev.off()
###

###
pdf(file="figures/serology_heatmap_cohorts_2.pdf",width = 8, height = 6)
ii <- serology_OF_meta_raw$Arm %in% c("Arm_D","Arm_E","Arm_F")
pheatmap(test[ii,],scale = 'none',show_rownames = FALSE,annotation_row =  serology_OF_meta_raw[,-2], angle_col = 90, color = viridis(255) )
dev.off()
###

###
pdf(file="figures/serology_heatmap_cohorts_3.pdf",width = 8, height = 6)
ii <- serology_OF_meta_raw$Arm %in% c("Arm_G","Arm_H","Arm_I")
pheatmap(test[ii,],scale = 'none',show_rownames = FALSE,annotation_row =  serology_OF_meta_raw[,-2], angle_col = 90, color = viridis(255) )
dev.off()
###

### for OD
grand_serology_od <- grand_serology_od[1:167,]
rownames(grand_serology_od) <- grand_serology_od$`Processed sample id`
grand_serology_od <- grand_serology_od[,-c(1:9)]
grand_serology_od <- grand_serology_od[,1:24]
###

###
grand_serology_od_log10 <- log10(grand_serology_od)
###

###
pdf(file="figures/serology_heatmap_all_cohorts_OD.pdf",width = 8, height = 6)
pheatmap(grand_serology_od_log10,scale = 'none',show_rownames = FALSE,annotation_row =  serology_OF_meta_raw, angle_col = 90, color = viridis(255) )
dev.off()

###
pdf(file="figures/serology_heatmap_cohort_1_OD.pdf",width = 8, height = 6)
ii <- serology_OF_meta_raw$Arm %in% c("Arm_A","Arm_B","Arm_C")
pheatmap(grand_serology_od_log10[ii,],scale = 'none',show_rownames = FALSE,annotation_row =  serology_OF_meta_raw, angle_col = 90, color = viridis(255) )
dev.off()
###

###
pdf(file="figures/serology_heatmap_cohorts_2_OD.pdf",width = 8, height = 6)
ii <- serology_OF_meta_raw$Arm %in% c("Arm_D","Arm_E","Arm_F")
pheatmap(grand_serology_od_log10[ii,],scale = 'none',show_rownames = FALSE,annotation_row =  serology_OF_meta_raw, angle_col = 90, color = viridis(255) )
dev.off()
###

###
pdf(file="figures/serology_heatmap_cohorts_3_OD.pdf",width = 8, height = 6)
ii <- serology_OF_meta_raw$Arm %in% c("Arm_G","Arm_H","Arm_I")
pheatmap(grand_serology_od_log10[ii,],scale = 'none',show_rownames = FALSE,annotation_row =  serology_OF_meta_raw, angle_col = 90, color = viridis(255) )
dev.off()
###


#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
###
### Serology Variance
###
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################


#############################################################################################################################################
### For Mac with 12 Cores
libs<-c("variancePartition","doParallel")
lapply(libs, require, character.only = TRUE)
registerDoParallel(makeCluster(10))
rm(libs)
#############################################################################################################################################
colnames(serology_OF_meta_raw)[3] <- "Collection_Event"
serology_OF_meta_raw$ParticipantID <- as.factor(as.character(serology_OF_meta_raw$ParticipantID))
### Form for model
form <- ~  (1|ParticipantID) + (1|Cohort) + (1|Collection_Event)
### run the model
variance_sero_od <- fitExtractVarPartModel(t(grand_serology_od_log10), form, serology_OF_meta_raw)
###
#############################################################################################################################################
pdf(file="figures/VP_OD_ee4412_serology.pdf",width = 4, height = 4)
plotVarPart( sortCols( variance_sero_od ) ) + ggtitle( 'OD Data') + theme(axis.text.x = element_text(angle = 60,hjust = 1))
dev.off()
#############################################################################################################################################

grand_serology_OF_log10 <- test
grand_serology_OF_log10_v2 <- grand_serology_OF_log10
grand_serology_OF_log10_v2[grand_serology_OF_log10_v2 < 2] <- 0

#############################################################################################################################################
### Form for model
form <- ~  (1|ParticipantID) + (1|Cohort) + (1|Collection_Event)
### run the model
variance_sero_OF <- fitExtractVarPartModel(t(grand_serology_OF_log10), form, serology_OF_meta_raw)
###
#############################################################################################################################################
pdf(file="figures/VP_OF_ee4412_serology.pdf",width = 4, height = 4)
plotVarPart( sortCols( variance_sero_OF ) ) + ggtitle( 'Original+Forecast Data') + theme(axis.text.x = element_text(angle = 60,hjust = 1))
dev.off()
#############################################################################################################################################

#############################################################################################################################################
identical( rownames(grand_serology_OF_log10), rownames(serology_OF_meta_raw) )
###
design <- model.matrix(~ 0 + Cohort:Collection , data = serology_OF_meta_raw )
my_design_names <- colnames(design)
colnames(design) <- c("C1_B","C2_B","C3_B","C1_Pre","C2_Pre","C3_Pre","C1_Res","C2_Res","C3_Res","C1_OS","C2_OS","C3_OS")
### model accounting for replicates
vfit <- lmFit(t(grand_serology_OF_log10_v2), design)
efit <- eBayes(vfit)
summary(decideTests(efit))
#############################################################################################################################################
contr.matrix <- makeContrasts( "C1vsC2_B" = C1_B - C2_B,
                               "C1vsC3_B" = C1_B - C3_B,
                               "C2vsC3_B" = C2_B - C3_B,
                               
                               "C1vsC2_OS" = C1_OS - C2_OS,
                               "C1vsC3_OS" = C1_OS - C3_OS,
                               "C2vsC3_OS" = C2_OS - C3_OS,
                               
                               "C1vsC2_Pre" = C1_Pre - C2_Pre,
                               "C1vsC3_Pre" = C1_Pre - C3_Pre,
                               "C2vsC3_Pre" = C2_Pre - C3_Pre,
                               
                               "C1vsC2_Res" = C1_Res - C2_Res,
                               "C1vsC3_Res" = C1_Res - C3_Res,
                               "C2vsC3_Res" = C2_Res - C3_Res,
                               
                               "C1vsC1_Base_pre" = C1_B - C1_Pre,
                               "C2vsC2_Base_pre" = C2_B - C2_Pre,
                               "C3vsC3_Base_pre" = C3_B - C3_Pre,
                               
                               "C1vsC1_Base_OS" = C1_B - C1_OS,
                               "C2vsC2_Base_OS" = C2_B - C2_OS,
                               "C3vsC3_Base_OS" = C3_B - C3_OS,
                               
                               "C1vsC1_Base_Res" = C1_B - C1_Res,
                               "C2vsC2_Base_Res" = C2_B - C2_Res,
                               "C3vsC3_Base_Res" = C3_B - C3_Res,
                               
                               "C1vsC1_Pre_OS" = C1_Pre - C1_OS,
                               "C2vsC2_Pre_OS" = C2_Pre - C2_OS,
                               "C3vsC3_Pre_OS" = C3_Pre - C3_OS,
                               
                               "C1vsC1_Res_OS" = C1_Res - C1_OS,
                               "C2vsC2_Res_OS" = C2_Res - C2_OS,
                               "C3vsC3_Res_OS" = C3_Res - C3_OS ,
                               
                               "C1vsC1_Pre_Res" = C1_Pre - C1_Res,
                               "C2vsC2_Pre_Res" = C2_Pre - C2_Res,
                               "C3vsC3_Pre_Res" = C3_Pre - C3_Res,
                               
                               "BaselinevsC1_Pre" = (C1_B+C2_B+C3_B)/3 - C1_Pre,
                               "BaselinevsC1_Res" = (C1_B+C2_B+C3_B)/3 - C1_Res,
                               "BaselinevsC1_OS" = (C1_B+C2_B+C3_B)/3 - C1_OS,
                               "BaselinevsC2_Pre" = (C1_B+C2_B+C3_B)/3 - C2_Pre,
                               "BaselinevsC2_Res" = (C1_B+C2_B+C3_B)/3 - C2_Res,
                               "BaselinevsC2_OS" = (C1_B+C2_B+C3_B)/3 - C2_OS,
                               "BaselinevsC3_Pre" = (C1_B+C2_B+C3_B)/3 - C3_Pre,
                               "BaselinevsC3_Res" = (C1_B+C2_B+C3_B)/3 - C3_Res,
                               "BaselinevsC3_OS" = (C1_B+C2_B+C3_B)/3 - C3_OS,
                               
                               levels = colnames(design) )
#############################################################################################################################################
vfit <- lmFit(t(grand_serology_OF_log10), design, method = "robust")
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE)
summary(decideTests(efit))
#############################################################################################################################################
lm_serology <- list()
for ( i in 1:length( colnames(summary(decideTests(efit))) ) ) {
  lm_serology[[i]] <- topTable(efit,coef= colnames(summary(decideTests(efit)))[i] ,n=Inf)  
  lm_serology[[i]]$Comparison <- colnames(summary(decideTests(efit)))[i]
  lm_serology[[i]]$Marker <- rownames(lm_serology[[i]])
}
lm_serology_all <- do.call(rbind,lm_serology)
#############################################################################################################################################

#############################################################################################################################################
### TEST FOR OD
#############################################################################################################################################
vfit <- lmFit(t(grand_serology_od_log10), design, method = "robust")
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE)
summary(decideTests(efit))
#############################################################################################################################################
lm_serology_od <- list()
for ( i in 1:length( colnames(summary(decideTests(efit))) ) ) {
  lm_serology_od[[i]] <- topTable(efit,coef= colnames(summary(decideTests(efit)))[i] ,n=Inf)  
  lm_serology_od[[i]]$Comparison <- colnames(summary(decideTests(efit)))[i]
  lm_serology_od[[i]]$Marker <- rownames(lm_serology_od[[i]])
}
lm_serology_od <- do.call(rbind,lm_serology_od)
#############################################################################################################################################
lm_serology_od[lm_serology_od$adj.P.Val<0.05,] ###  Much lower values ...
#############################################################################################################################################


#############################################################################################################################################
dim(lm_serology_all)
table(lm_serology_all$adj.P.Val<0.05) ## 20
# lm_serology_all[lm_serology_all$adj.P.Val<0.05,]
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
###
### HEATMAP for SEROLOGY
###
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

annotation_row <- lm_serology_all[lm_serology_all$adj.P.Val<0.05,]

annotation_row <- as.matrix(table(lm_serology_all$Marker[lm_serology_all$adj.P.Val<0.05],lm_serology_all$Comparison[lm_serology_all$adj.P.Val<0.05]))
class(annotation_row) <- "matrix"
annotation_row <- as.data.frame(annotation_row)

#iy <- colnames(grand_serology_OF_log10) %in% unique(sort(lm_serology_all[lm_serology_all$adj.P.Val<0.05,]$Marker))

annotation_colors_serology <- list( Cohort = c( `1_BV+ipi`= brewer.pal(7, "Dark2")[1], 
                                       `2_BV+nivo` = brewer.pal(7, "Dark2")[2], 
                                       `3_BV+nivo+ipi` = brewer.pal(7, "Dark2")[3] ) , 
                           Collection = c( Baseline = "grey30",
                                           Pre_D1_C2 = "grey50",
                                           Restaging = "grey70",
                                           Off_study = "grey90"),
                           BaselinevsC1_OS = c( `1`="white", `0` = "black"),
                           BaselinevsC1_Pre = c( `1`="white", `0` = "black"),
                           BaselinevsC1_Res = c( `1`="white", `0` = "black"),
                           BaselinevsC2_OS = c( `1`="white", `0` = "black"),
                           BaselinevsC2_Pre = c( `1`="white", `0` = "black"),
                           BaselinevsC2_Res = c( `1`="white", `0` = "black"),
                           BaselinevsC3_OS = c( `1`="white", `0` = "black"),
                           BaselinevsC3_Pre = c( `1`="white", `0` = "black"),
                           BaselinevsC3_Res = c( `1`="white", `0` = "black"),
                           C1vsC3_B = c( `1`="white", `0` = "black"),
                           C1vsC3_Pre = c( `1`="white", `0` = "black"),
                           C2vsC3_B = c( `1`="white", `0` = "black"),
                           C2vsC3_Pre = c( `1`="white", `0` = "black"),
                           C3vsC3_Base_OS = c( `1`="white", `0` = "black"),
                           C3vsC3_Base_pre = c( `1`="white", `0` = "black"),
                           C3vsC3_Base_Res  = c( `1`="white", `0` = "black"),
                           C3vsC3_Pre_OS = c( `1`="white", `0` = "black"),
                           C3vsC3_Pre_Res = c( `1`="white", `0` = "black")
                           )

annotation_colors_serology$Arm <- c( Arm_A = "grey90", Arm_B = "grey80", Arm_C = "grey70",
                                     Arm_D = "grey60", Arm_E = "grey50", Arm_F = "grey40",
                                     Arm_G = "grey30", Arm_H = "grey20", Arm_I = "grey10")

names(annotation_colors_serology$Collection) <- names(table(serology_OF_meta_raw$Collection))

identical(rownames(grand_serology_OF_log10),rownames(serology_OF_meta_raw))

iy <- order(paste(serology_OF_meta_raw$Arm,serology_OF_meta_raw$Collection))

pdf(file="figures/serology_heatmap_all_cohorts_DE.pdf",width = 14, height = 6)
pheatmap(t(grand_serology_OF_log10_v2[iy,]),
         scale = 'none',
         cluster_rows = TRUE, 
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_row = annotation_row, 
         annotation_colors = annotation_colors_serology,
         annotation_col =  serology_OF_meta_raw[,-2], 
         angle_col = 90, color = viridis(255))
dev.off()

table(is.na(t(grand_serology_OF_log10_v2[iy,])))
  
pdf(file="figures/serology_heatmap_all_cohorts_DE_col.pdf",width = 14, height = 6)
pheatmap(t(grand_serology_OF_log10_v2[iy,]),
         scale = 'column',
         cluster_rows = TRUE, 
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_row = annotation_row, 
         annotation_colors = annotation_colors_serology,
         annotation_col =  serology_OF_meta_raw[,-2], angle_col = 90, 
         color = colorRampPalette(c("steelblue","white","firebrick"))(10))
#color = viridis(255) )
dev.off()


#############################################################################################################################################
#############################################################################################################################################

###
generic_meta <- serology_OF_meta_raw
###
generic_meta$Sample <- rownames(generic_meta)
###
generic_var <- melt(grand_serology_OF_log10_v2)
colnames(generic_var)<- c("Sample","Antigen","Log10titers")
###
generic_var <- merge(generic_var,generic_meta,by="Sample")
###
#head(generic_var)
#table(lm_serology_all$adj.P.Val<0.05)
table( generic_var$Antigen %in% unique(sort(lm_serology_all$Marker[lm_serology_all$adj.P.Val<0.05])) )
ii <- generic_var$Antigen %in% unique(sort(lm_serology_all$Marker[lm_serology_all$adj.P.Val<0.05])) & generic_var$Log10titers > 1
###
generic_var$Collection <- as.factor(generic_var$Collection)
generic_var$Time_point <- as.numeric(generic_var$Collection)
### table(generic_var$Collection,generic_var$Time_point)
table(generic_var$Collection,generic_var$Time_point)
head(generic_var)
###
library(beeswarm)

pdf(file="figures/serology_line_plot_DE_markers.pdf",width = 10, height = 10)
ggarrange(
ggplot(data=generic_var[ii,]) + aes(x=Collection,y=Log10titers,color=Cohort) + theme_bw() +
  facet_wrap(~Antigen,ncol=1,nrow=2)+
  geom_boxplot() +
  geom_quasirandom()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  theme(legend.position="bottom")
,
ggplot(data=generic_var[ii,]) + aes(x=Time_point,y=Log10titers,color=Cohort) + theme_bw() +
  facet_wrap(~Antigen,ncol=1,nrow=2)+
  stat_smooth(method="glm",fullrange = TRUE,alpha=0,show.legend = FALSE) + 
  #geom_boxplot() +
  geom_quasirandom()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  theme(legend.position="bottom")
,
ncol=2,nrow=1,align = 'hv')
dev.off()
#lm, glm, gam, loess, rlm

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
###
### Variance Analysis cytoff
###
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
###
cimac_e4412_metadata$Collection_Event <- cimac_e4412_metadata$Timepoint
###
cytof_big_matrix <- dcast(data = channel_intesity_statistics,formula = ChannelName~sorter, value.var = "Q0_95")
rownames(cytof_big_matrix) <- cytof_big_matrix$ChannelName
cytof_big_matrix <- cytof_big_matrix[,-1]
###
cytof_big_matrix_metadata <- data.frame(  ParticipantID = tidyr::separate(data.frame( colnames(cytof_big_matrix) ), 1, sep="_", c("a","b"))$a , 
                  CellSubset = tidyr::separate(data.frame( colnames(cytof_big_matrix) ), 1, sep="_", c("a","b"))$b )
###
test <- as.data.frame(channel_intesity_statistics[,colnames(channel_intesity_statistics) %in% c("SampleId","Name")])
colnames(test)[1]<-"ParticipantID"
test <- unique(test)
cytof_big_matrix_metadata <- merge(cytof_big_matrix_metadata,test,by="ParticipantID")
###
test <- cimac_e4412_metadata [ , colnames(cimac_e4412_metadata) %in% c("astrolabe_ID","Cohort","Collection_Event") ]
test <- unique(test)
colnames(test) <- c("Participant_ID","Cohort","Collection_Event")
dim(cytof_big_matrix_metadata)
###
ii <- cytof_big_matrix_metadata$ParticipantID %in% test$Participant_ID
cytof_big_matrix <- cytof_big_matrix[,ii]
cytof_big_matrix_metadata <- cytof_big_matrix_metadata[ii,]
cytof_big_matrix_metadata$Cohort <- "empty"
cytof_big_matrix_metadata$Collection_Event <- "empty"
###
for (i in 1:nrow(cytof_big_matrix_metadata)){
  cytof_big_matrix_metadata$Cohort[i] <- test$Cohort[ test$Participant_ID %in% cytof_big_matrix_metadata$ParticipantID[i] ]
  cytof_big_matrix_metadata$Collection_Event[i] <- test$Collection_Event[ test$Participant_ID %in% cytof_big_matrix_metadata$ParticipantID[i] ]
}
#############################################################################################################################################
### Form for model
form <- ~ (1|ParticipantID) + (1|CellSubset) + (1|Cohort) + (1|Collection_Event) 
#############################################################################################################################################
cytof_big_matrix[is.na(cytof_big_matrix)] <- 0
### run the model
vp_cytof_exprs_matrix <- fitExtractVarPartModel(cytof_big_matrix, form, cytof_big_matrix_metadata)
#############################################################################################################################################
pdf(file="figures/VP_ee4412_CYTOF.pdf",width = 4, height = 4)
plotVarPart( sortCols( vp_cytof_exprs_matrix ) ) + ggtitle( 'Cytof Data') + theme(axis.text.x = element_text(angle = 60,hjust = 1))
dev.off()
#############################################################################################################################################


#############################################################################################################################################
#############################################################################################################################################

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
###
### Correlations
###
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

### Olink
### cimac_e4412_olink
### cimac_e4412_olink_clustergram_meta

## cytof
## cytof_big_matrix
## cytof_big_matrix_metadata

# Serology
# grand_serology_OF_log10_v2
# serology_OF_meta_raw

#############################################################################################################################################

colnames(cimac_e4412_olink_clustergram_meta) %in% c("Cohort","Collection_Event")

table(cytof_big_matrix_metadata$ParticipantID %in% cimac_e4412_olink_clustergram_meta$ParticipantID) 
table(cytof_big_matrix_metadata$ParticipantID %in% as.character(serology_OF_meta_raw$ParticipantID) )
cytof_big_matrix_metadata$something <- tidyr::separate(data.frame( cytof_big_matrix_metadata$Name ), 1, sep="_", c("a","b"))$a

cytof_big_matrix_metadata$ParticipantID <- NULL
colnames(cytof_big_matrix_metadata)[5] <- "ParticipantID"

table(cytof_big_matrix_metadata$ParticipantID %in% as.character(serology_OF_meta_raw$ParticipantID) )
table( as.character(serology_OF_meta_raw$ParticipantID) %in% cimac_e4412_olink_clustergram_meta$ParticipantID )
table( as.character(cytof_big_matrix_metadata$ParticipantID) %in% cimac_e4412_olink_clustergram_meta$ParticipantID )

grand_serology_OF_log10_v3 <- t(grand_serology_OF_log10_v2)

identical(colnames(grand_serology_OF_log10_v3),rownames(serology_OF_meta_raw)) ### check
colnames(grand_serology_OF_log10_v3) <- serology_OF_meta_raw$ParticipantID

#############################################################################################################################################
#############################################################################################################################################

dim(cytof_big_matrix)
dim(cimac_e4412_olink)
dim(grand_serology_OF_log10_v3)

#############################################################################################################################################
### remove duplicates
ii <- ! duplicated(serology_OF_meta_raw)
serology_OF_meta_raw[!ii,]
serology_OF_meta_raw[serology_OF_meta_raw$ParticipantID %in% "8269916",]


grand_serology_OF_log10_v3 <- grand_serology_OF_log10_v3[,ii]
serology_OF_meta_raw <- serology_OF_meta_raw[ii,]

identical(colnames(grand_serology_OF_log10_v3),as.character(serology_OF_meta_raw$ParticipantID)) # check
###
serology_OF_meta_raw$ParticipantID <- as.character(serology_OF_meta_raw$ParticipantID) 
#############################################################################################################################################

rownames(cytof_big_matrix_metadata) <- colnames(cytof_big_matrix)

ii <- ! cytof_big_matrix_metadata$Collection_Event %in% c("Off_study","Other")
cytof_big_matrix_metadata <- cytof_big_matrix_metadata[ii,]
cytof_big_matrix <- cytof_big_matrix[,ii]

#############################################################################################################################################

ii <- cytof_big_matrix_metadata$CellSubset %in% cytof_big_matrix_metadata$CellSubset[1]
dim(cytof_big_matrix[,ii])
dim(serology_OF_meta_raw)

table(serology_OF_meta_raw$ParticipantID %in% cytof_big_matrix_metadata$ParticipantID[ii])
table(cytof_big_matrix_metadata$ParticipantID[ii] %in% serology_OF_meta_raw$ParticipantID)

ix <- serology_OF_meta_raw$ParticipantID %in% cytof_big_matrix_metadata$ParticipantID[ii]

serology_OF_meta_raw$ParticipantID[!ix]

serology_OF_meta_raw <- serology_OF_meta_raw[ix,]
grand_serology_OF_log10_v3 <- grand_serology_OF_log10_v3[,ix]

table(serology_OF_meta_raw$Collection_Event)
table(cytof_big_matrix_metadata$Collection_Event[ii])
table(cimac_e4412_olink_clustergram_meta$Collection_Event)

serology_OF_meta_raw$Collection_Event[serology_OF_meta_raw$Collection_Event %in% "1.Baseline"] <- "Baseline"
serology_OF_meta_raw$Collection_Event[serology_OF_meta_raw$Collection_Event %in% "2.Pre_Day_1_Cycle_2"] <- "Pre_Day_1_Cycle_2"
serology_OF_meta_raw$Collection_Event[serology_OF_meta_raw$Collection_Event %in% "4.Off_study"] <- "End_of_treatment"
serology_OF_meta_raw$Collection_Event[serology_OF_meta_raw$Collection_Event %in% "3.Restaging"] <- "Restaging"

cimac_e4412_olink_clustergram_meta$Collection_Event[cimac_e4412_olink_clustergram_meta$Collection_Event %in% "Pre_D1_C2"] <- "Pre_Day_1_Cycle_2"
cimac_e4412_olink_clustergram_meta$Collection_Event[cimac_e4412_olink_clustergram_meta$Collection_Event %in% "Off_study"] <- "End_of_treatment"

table(serology_OF_meta_raw$Collection_Event)
table(cytof_big_matrix_metadata$Collection_Event[ii])
table(cimac_e4412_olink_clustergram_meta$Collection_Event)

table(serology_OF_meta_raw$Cohort)
table(cytof_big_matrix_metadata$Cohort[ii])
table(cimac_e4412_olink_clustergram_meta$Cohort)

serology_OF_meta_raw$max_sorter <- paste(serology_OF_meta_raw$ParticipantID,serology_OF_meta_raw$Collection_Event,serology_OF_meta_raw$Cohort,sep=".")
cytof_big_matrix_metadata$max_sorter <- paste(cytof_big_matrix_metadata$ParticipantID,cytof_big_matrix_metadata$Collection_Event,cytof_big_matrix_metadata$Cohort,sep=".")
cimac_e4412_olink_clustergram_meta$max_sorter <- paste(cimac_e4412_olink_clustergram_meta$ParticipantID,cimac_e4412_olink_clustergram_meta$Collection_Event,cimac_e4412_olink_clustergram_meta$Cohort,sep=".")

big_list <- list( serology= unique(sort(serology_OF_meta_raw$max_sorter)),
                  cytof= unique(sort(cytof_big_matrix_metadata$max_sorter)),
                  olink= unique(sort(cimac_e4412_olink_clustergram_meta$max_sorter)) )

my_big_list <- Reduce(intersect,big_list)

ii <- serology_OF_meta_raw$max_sorter %in% my_big_list
serology_OF_meta_raw <- serology_OF_meta_raw[ii,]
grand_serology_OF_log10_v3 <- grand_serology_OF_log10_v3[,ii]

dim(serology_OF_meta_raw)
dim(grand_serology_OF_log10_v3)

ii <- cimac_e4412_olink_clustergram_meta$max_sorter %in% my_big_list
onlink_metadata <- cimac_e4412_olink_clustergram_meta[ii,]
olink_npx_matrix <- cimac_e4412_olink[,ii]

ii <- ! duplicated(onlink_metadata$max_sorter)
onlink_metadata <- onlink_metadata[ii,]
olink_npx_matrix <- olink_npx_matrix[,ii]

ii <- cytof_big_matrix_metadata$CellSubset %in% cytof_big_matrix_metadata$CellSubset[1]
dim(cytof_big_matrix[,ii])
dim(cytof_big_matrix_metadata[ii,])

ix <- cytof_big_matrix_metadata$max_sorter %in% my_big_list
table(ix)
cytof_big_matrix <- cytof_big_matrix[,ix]
cytof_big_matrix_metadata <- cytof_big_matrix_metadata[ix,]

### 145 !!!
ii <- cytof_big_matrix_metadata$CellSubset %in% cytof_big_matrix_metadata$CellSubset[1]
dim(cytof_big_matrix[,ii])
dim(cytof_big_matrix_metadata[ii,])

colnames(cytof_big_matrix)[ii]
cytof_big_matrix_metadata$ParticipantID[ii]

my_levels <- levels(as.factor(cytof_big_matrix_metadata$CellSubset))

my_new_matrix <- list()
for (i in 1:length(my_levels)){
  ix <- cytof_big_matrix_metadata$CellSubset %in%  my_levels[i]
  my_new_matrix[[i]] <- cytof_big_matrix[,ix]
  names(my_new_matrix)[i] <- gsub(" ",".",my_levels[i])
  rownames(my_new_matrix[[i]]) <- paste(gsub(" ",".",my_levels[i]),rownames(my_new_matrix[[i]]),sep="-")
  colnames(my_new_matrix[[i]]) <- tidyr::separate(data.frame( colnames(my_new_matrix[[i]]) ), 1, sep="_", c("a","b"))$a
}

#############################################################################################################################################
names(my_new_matrix) <- NULL
my_new_matrix <- do.call(rbind,my_new_matrix)

dim(cytof_big_matrix_metadata) #check

unique(cytof_big_matrix_metadata$max_sorter)[1]
unique(rownames(cytof_big_matrix_metadata))[1]

colnames(my_new_matrix) <- unique(cytof_big_matrix_metadata$max_sorter)

identical(colnames(olink_npx_matrix), onlink_metadata$Sample_ID)
colnames(olink_npx_matrix)  <- onlink_metadata$max_sorter

identical(colnames(grand_serology_OF_log10_v3),serology_OF_meta_raw$ParticipantID)

colnames(grand_serology_OF_log10_v3) <- serology_OF_meta_raw$max_sorter

dim(onlink_metadata)
dim(serology_OF_meta_raw) #check

#############################################################################################################################################
#############################################################################################################################################
###
### save(file="big_matrix_cross_experiment.RData",my_new_matrix,olink_npx_matrix,grand_serology_OF_log10_v3,onlink_metadata,serology_OF_meta_raw)
load("big_matrix_cross_experiment.RData")
###
#############################################################################################################################################
#############################################################################################################################################
dim(my_new_matrix)
rownames(my_new_matrix) <- paste("cytof",rownames(my_new_matrix),sep="-")
rownames(olink_npx_matrix) <- paste("olink",rownames(olink_npx_matrix),sep="-")
rownames(grand_serology_OF_log10_v3) <- paste("serology",rownames(grand_serology_OF_log10_v3),sep="-")
#############################################################################################################################################
#############################################################################################################################################

my_new_matrix <- my_new_matrix[ , order(colnames(my_new_matrix)) ]
olink_npx_matrix <- olink_npx_matrix[ , order(colnames(olink_npx_matrix)) ]
grand_serology_OF_log10_v3 <- grand_serology_OF_log10_v3[ , order(colnames(grand_serology_OF_log10_v3)) ]

identical( colnames(my_new_matrix) , colnames(olink_npx_matrix) )
identical( colnames(my_new_matrix) , colnames(grand_serology_OF_log10_v3) )

big_cross_matrix <- rbind(my_new_matrix,olink_npx_matrix,grand_serology_OF_log10_v3)

dim(big_cross_matrix)
#############################################################################################################################################
write.csv(file="big_cross_matrix.csv",big_cross_matrix)
#############################################################################################################################################
###
big_cross_matrix[is.na(big_cross_matrix)] <- 0

### Correlation across samples
p_mat_df <- cor.mtest(big_cross_matrix,method = "spearman")
p_mat_df$r <- cor(big_cross_matrix,method = "spearman")
###

dim(p_mat_df$r)

####
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame( row = rownames(cormat)[row(cormat)[ut]],
              column = rownames(cormat)[col(cormat)[ut]],
              cor  =(cormat)[ut],
              p = pmat[ut] ) }
####
flat_p_mat_df <- flattenCorrMatrix(p_mat_df$r,p_mat_df$p)


####
annotation_p_mat_df <- data.frame(
  #Patient = tidyr::separate(data.frame( colnames(p_mat_df$r) ), 1, sep="\\.", c("a","b","c"))$a
  Collection_Event = tidyr::separate(data.frame( colnames(p_mat_df$r) ), 1, sep="\\.", c("a","b","c"))$b,
  Cohort = tidyr::separate(data.frame( colnames(p_mat_df$r) ), 1, sep="\\.", c("a","b","c"))$c
)
rownames(annotation_p_mat_df) <- colnames(p_mat_df$r)
####

annotation_color_p_mat_df <- list( 
  Cohort = c( `1_BV+ipi`= brewer.pal(7, "Dark2")[1],`2_BV+nivo` = brewer.pal(7, "Dark2")[2],  `3_BV+nivo+ipi` = brewer.pal(7, "Dark2")[3] ) , 
  Collection_Event = c( Baseline = "grey30", Pre_Day_1_Cycle_2 = "grey50", Restaging = "grey70", End_of_treatment = "grey90")
)

pdf(file="figures/heatmap_spearman_corr_samples.pdf",width = 10, height = 8)
pheatmap( p_mat_df$r, 
          main = "",
          color = colorRampPalette(c("steelblue","white","firebrick"))(255),         
          cluster_cols = T, 
          cluster_rows = T, 
          show_colnames = FALSE,
          border_color = FALSE, 
          show_rownames = FALSE, 
          annotation_col = annotation_p_mat_df,
          annotation_row = annotation_p_mat_df,
          annotation_colors = annotation_color_p_mat_df,
          clustering_distance_rows ="euclidean", 
          scale="none")
dev.off()

#############################################################################################################################################
#############################################################################################################################################

### Correlation across features
p_mat_feat_df <- cor.mtest(t(big_cross_matrix),method = "spearman")
p_mat_feat_df$r <- cor(t(big_cross_matrix),method = "spearman")
###

dim(p_mat_feat_df$r)
dim(p_mat_feat_df$p)
table(p_mat_feat_df$p<0.05)

###
flat_p_mat_feat_df <- flattenCorrMatrix(p_mat_feat_df$r,p_mat_feat_df$p)
flat_p_mat_feat_df$fdr <- p.adjust(flat_p_mat_feat_df$p,method="BH")
###
table(is.na(p_mat_feat_df$r))
p_mat_feat_df$r[is.na(p_mat_feat_df$r)]<- 0
###
annotation_p_mat_feat_df <- data.frame( Experiment = tidyr::separate(data.frame( colnames(p_mat_feat_df$r) ), 1, sep="-", c("a","b","c"))$a )
annotation_color_p_mat_feat_df <- list( Experiment = c( olink= brewer.pal(7, "Dark2")[1], cytof= brewer.pal(7, "Dark2")[2],  serology = brewer.pal(7, "Dark2")[3] ))
###
rownames(annotation_p_mat_feat_df) <- colnames(p_mat_feat_df$r)

generic_var <- p_mat_feat_df$r
generic_var[ p_mat_feat_df$p < 0.05 ] <- NA

###
pdf(file="figures/heatmap_spearman_corr_features_hclust.pdf",width = 10, height = 8)
out <- pheatmap( p_mat_feat_df$r, 
          main = "",
          color = colorRampPalette(c("steelblue","white","firebrick"))(255),         
          cluster_cols = T, 
          cluster_rows = T, 
          show_colnames = FALSE,
          border_color = FALSE, 
          show_rownames = FALSE, 
          annotation_col = annotation_p_mat_feat_df,
          annotation_row = annotation_p_mat_feat_df,
          annotation_colors = annotation_color_p_mat_feat_df,
          clustering_distance_rows ="euclidean", 
          scale="none")
dev.off()
###

###
generic_var <- flat_p_mat_feat_df[flat_p_mat_feat_df$fdr<0.05 & abs(flat_p_mat_feat_df$cor) > 0.5 ,]
###
#generic_var$row <- factor( generic_var$row, levels = out$tree_col$labels[out$tree_col$order] )
#generic_var$column <- factor( generic_var$column, levels = out$tree_row$labels[out$tree_row$order] )
###
generic_var$nLogFDR <- -log10(generic_var$fdr)
###
pdf(file="figures/heatmap_spearman_corr_features.pdf",width = 36, height = 32)
ggplot(generic_var, aes(row, column, color= cor, size=nLogFDR )) + 
  geom_point(shape=15)+
  scale_colour_gradient2(low='firebrick', mid="white",high='steelblue', name='Rho') +
  theme_bw() +
  rotate_x_text(angle=90) + 
  labs(x ='', y='', title='') 
dev.off()

#############################################################################################################################################
#############################################################################################################################################

### set feature x sample

#############################################################################################################################################
#############################################################################################################################################

####
annotation_col_bcm <- data.frame( Collection_Event = tidyr::separate(data.frame( colnames(big_cross_matrix) ), 1, sep="\\.", c("a","b","c"))$b,
                                  Cohort = tidyr::separate(data.frame( colnames(big_cross_matrix) ), 1, sep="\\.", c("a","b","c"))$c )
####
rownames(annotation_col_bcm) <- colnames(big_cross_matrix)
####
annotation_row_bcm <- data.frame( Experiment = tidyr::separate(data.frame( rownames(big_cross_matrix) ), 1, sep="-", c("a","b","c"))$a )
rownames(annotation_row_bcm) <- rownames(big_cross_matrix)
####
annotation_color_bcm <- list( 
  Experiment = c( olink= "blue", cytof= "red",  serology = "green" ) ,
  Collection_Event = c(Baseline = "grey30", Pre_Day_1_Cycle_2 = "grey50", Restaging = "grey70", End_of_treatment = "grey90") ,
  Cohort = c( `1_BV+ipi`= brewer.pal(7, "Dark2")[1],`2_BV+nivo` = brewer.pal(7, "Dark2")[2],  `3_BV+nivo+ipi` = brewer.pal(7, "Dark2")[3] )  
  )
####

# big_cross_matrix #
###
generic_var <- as.matrix(big_cross_matrix)
###
for ( i in 1:nrow(big_cross_matrix)){
  generic_var[i,] <- as.numeric(scale(as.numeric(big_cross_matrix[i,])))
}
dim(generic_var)
summary(c(generic_var))
generic_var[is.na(generic_var)] <- 0
###
my_ids <- list()
for (i in 1:length(c(astrolabe_markers,"olink","serology") ) ) { my_ids[[i]] <- grep(c(astrolabe_markers,"olink","serology")[i], rownames(big_cross_matrix)) }

ii <- unique(sort(unlist(my_ids)))
generic_var <- as.matrix(big_cross_matrix[ii,])

my_sds <- list()
for( i in 1:nrow(generic_var)){
  my_sds[[i]] <- sd(generic_var[i,])  
}
my_sds <- unlist(my_sds)

ii <- which(my_sds > 0)
generic_var <- generic_var[ii,]


quantile_breaks <- function(xs, n = 20) { breaks <- quantile(xs, probs = seq(0, 1, length.out = n)) ; breaks[!duplicated(breaks)] }
mat_breaks1 <- quantile_breaks( as.matrix(generic_var) , n = 10)
mat_breaks <- seq(min(generic_var), max(generic_var), length.out = 20)

generic_var2 <- t(scale(t(generic_var)))

#############################################################################################################################################
pdf(file="figures/heatmap_spearman_corr_big_cross_matrix_hclust_mini.pdf",width = 8, height = 6)
my_hclust_big_cross_matrix <- pheatmap( generic_var2 , main = "",
                 color = colorRampPalette(c("steelblue","white","firebrick"))(10),         
                 cluster_cols = T, 
                 cluster_rows = T, 
                 show_colnames = FALSE,
                 border_color = FALSE, 
                 show_rownames = FALSE,
                 breaks = mat_breaks1,
                 annotation_col = annotation_col_bcm,
                 annotation_row = annotation_row_bcm,
                 annotation_colors = annotation_color_bcm,
                 clustering_distance_rows ="euclidean", 
                 scale="none")
dev.off()
#############################################################################################################################################

### write.csv(file="hclust_matrix_mini.csv",generic_var)
###
mydata <- prcomp( t(generic_var) , scale=TRUE, center = TRUE ) 
###
fviz_eig(mydata, addlabels = T, barfill = "lightsteelblue3", barcolor = "lightsteelblue3") + theme_bw() + ggtitle("")
###
pca_mts <- as.data.frame(mydata$x)
pca_mts$Sample_ID <- rownames(pca_mts)
###
umap_clustering <- umap(t(generic_var))
pca_mts$UMAP1 <- umap_clustering$layout[,1]
pca_mts$UMAP2 <- umap_clustering$layout[,2]
###
tsne_clustering <- Rtsne(t(generic_var), dims = 2, initial_dims = 30,
                         perplexity = 15, theta = 1, check_duplicates = TRUE,
                         pca = TRUE, partial_pca = FALSE, max_iter = 5000,
                         normalize = TRUE, momentum = 0.5, final_momentum = 0.8, eta = 200, exaggeration_factor = 12, num_threads = 1)
###      
pca_mts$tsne1 <- tsne_clustering$Y[,1]
pca_mts$tsne2 <- tsne_clustering$Y[,2]
###
### over PCA
umap_clustering <- umap(pca_mts[,1:30])
pca_mts$UMAP1_s <- umap_clustering$layout[,1]
pca_mts$UMAP2_s <- umap_clustering$layout[,2]
### over PCA
tsne_clustering <- Rtsne(pca_mts[,1:30], dims = 2, initial_dims = 30,
                         perplexity = 15, theta = 1, check_duplicates = FALSE,
                         pca = FALSE, partial_pca = FALSE, max_iter = 5000,
                         normalize = TRUE, momentum = 0.5, final_momentum = 0.8, eta = 200, 
                         exaggeration_factor = 12, num_threads = 1)
### 
pca_mts$tsne1_s <- tsne_clustering$Y[,1]
pca_mts$tsne2_s <- tsne_clustering$Y[,2]

###
test <- annotation_col_bcm
test$Sample_ID <- rownames(annotation_col_bcm)

### merge to add metadata columns
pca_mts <- merge(pca_mts, test, by="Sample_ID")
pca_mts$Condition <- paste(pca_mts$Cohort,pca_mts$Collection_Event,sep="-")

#############################################################################################################################################
#pdf(file="figures/reductions_all_data_ee4412_olink.pdf",width = 14, height = 10)
ggarrange( ggplot(pca_mts, aes(PC1, PC2, color = Condition)) + theme_bw() + geom_point(alpha=0.69) ,
           ggplot(pca_mts, aes(UMAP1, UMAP2, color = Condition)) + theme_bw() + geom_point(alpha=0.69) , 
           fviz_eig(mydata, addlabels = T, barfill = "lightsteelblue3", barcolor = "lightsteelblue3") + theme_bw() + ggtitle(""), 
           ggplot(pca_mts, aes(tsne1, tsne2, color = Condition)) + theme_bw() + geom_point(alpha=0.69) , 
           ncol=2,nrow=2 )
dev.off()
#############################################################################################################################################
#############################################################################################################################################
big_mini_cross_matrix <- generic_var
#############################################################################################################################################
#############################################################################################################################################

#save(file="big_matrix_cross_experiment.RData",my_new_matrix,olink_npx_matrix,grand_serology_OF_log10_v3,onlink_metadata,serology_OF_meta_raw,
#my_hclust_big_cross_matrix,big_mini_cross_matrix)
#load("big_matrix_cross_experiment.RData")

pdf(file="figures/tree_mini_big_cross_matrix_tree.pdf",width = 24, height = 12)
plot(my_hclust_big_cross_matrix$tree_col)
dev.off()

#############################################################################################################################################
#############################################################################################################################################
rownames(lm_serology_all)<-NULL
###
setwd("/Users/gonzae34/Documents/projects_gnjatic/U24/U24/E4412_Diefenbach")
write.csv(file="olink_e4412_results_table_eegk.csv",lm_olink_all[lm_olink_all$adj.P.Val<0.05,])
# write.csv(file="cytof_cell_frequencies_e4412_results_table_eegk.csv",lm_frequencies_all)
write.csv(file="cytof_cell_frequencies_e4412_results_table_eegk.csv",lm_frequencies_all[lm_frequencies_all$adj.P.Val<0.15,])
#write.csv(file="cytof_cell_markers_e4412_results_table_eegk.csv",lm_cytof_markers_all)
write.csv(file="cytof_cell_markers_e4412_results_table_eegk.csv",lm_cytof_markers_all[lm_cytof_markers_all$adj.P.Val<0.05,])
#write.csv(file="cytof_cell_markers_e4412_cohort_results_table_eegk.csv",lm_cytof_markers_cohort_all)
write.csv(file="cytof_cell_markers_e4412_cohort_results_table_eegk.csv",lm_cytof_markers_cohort_all[lm_cytof_markers_cohort_all$adj.P.Val<0.05,])
###
write.csv(file="serology_e4412_results_table_eegk.csv",lm_serology_all[ lm_serology_all$adj.P.Val<0.05,])
###
#############################################################################################################################################
#############################################################################################################################################
###
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
rm(my_plot,my_cytof_plots,scree_pca,pca_mts,mydata)
### save.image("e4412_analysis.RData")
### load("e4412_analysis.RData")
### load("big_matrix_cross_experiment.RData")
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################