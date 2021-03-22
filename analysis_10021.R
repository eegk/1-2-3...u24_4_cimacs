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

### Load data
load("/Users/gonzae34/Documents/projects_gnjatic/U24/rdata_files/cimac_10021_raw_data.RData")
#############################################################################################################################################


### Preliminary analysis DE
#############################################################################################################################################

### files
dim(cimac_10021_samples) # samples data
dim(cimac_10021_npx_raw_controls) # control data
dim(cimac_10021_npx) # original data
dim(cimac_10021_protein_info) # protein info
dim(cimac_10021_npx_raw) # data matrix
################################################


### relevant Variables 
###################################################################################
### conditions
table(cimac_10021_samples$condition)
table(cimac_10021_samples$`Plate ID`)
table(cimac_10021_samples$`Plate ID`)
###head(cimac_10021_samples)
cimac_10021_samples$condition <- cimac_10021_samples$`Collection event name`
cimac_10021_samples$plate <- cimac_10021_samples$`Plate ID`
cimac_10021_samples$cohort <- cimac_10021_samples$`Cohort name`
cimac_10021_samples$patientID <- cimac_10021_samples$`Participant id`

table((cimac_10021_samples$`Cohort name`))
cimac_10021_samples$`Material storage condition`

cimac_10021_samples$condition[is.na(cimac_10021_samples$condition)] <- "Other"

### VariancePartition
###################################################################################

identical(colnames(cimac_10021_npx_raw) , cimac_10021_samples$`Cimac id`)
cimac_10021_npx_raw[is.na(cimac_10021_npx_raw)] <-0

ix <- grep("Ctrl",rownames(cimac_10021_npx_raw)) ### no patterns ### remove not select from data.

### Form for model
form <- ~ (1|patientID) + (1|plate) + (1|condition) + (1|cohort)

### run the model
variance_exprs_matrix <- fitExtractVarPartModel(cimac_10021_npx_raw[-ix,], form, cimac_10021_samples)

### Store the figure
###################################################################################
pdf(file="/Users/gonzae34/Documents/projects_gnjatic/U24/figures_u24/10021_figures/s10021_VP_plot.pdf",width = 4, height = 4)
plotVarPart( sortCols( variance_exprs_matrix ) ) + ggtitle( 'NPX Data') + theme(axis.text.x = element_text(angle = 60,hjust = 1))
dev.off()
###################################################################################



### DE analysis
###################################################################################

rownames(cimac_10021_samples) <- cimac_10021_samples$`Cimac id`
### Design matrix for the model
design <- model.matrix( ~ 0 + condition + cohort + plate, data = cimac_10021_samples )
colnames(design) <- gsub("condition","",colnames(design))
colnames(design) <- gsub("cohort","",colnames(design))
colnames(design) <- gsub("plate","",colnames(design))
colnames(design) <- make.names(colnames(design))

### contrast for comparisons
contr.matrix <- makeContrasts("Baseline vs On_Treatment" = Baseline - On_Treatment, 
                              "Baseline vs Other" = Baseline - Other, 
                              "On_Treatment vs Other" = On_Treatment - Other, 
                              
                              "NSCLC_A vs CRC_B" = NSCLC_A - CRC_B, 
                              "NSCLC_A vs NSCLC_B" = NSCLC_A - NSCLC_B, 
                              "NSCLC_A vs NSCLC_C" = NSCLC_A - NSCLC_C, 
                              
                              "CRC_B vs NSCLC_B" = CRC_B - NSCLC_B, 
                              "CRC_B vs NSCLC_C" = CRC_B - NSCLC_B, 
                              
                              "NSCLC_B vs NSCLC_C" = NSCLC_B - NSCLC_C, 
                              
                              levels = colnames(design) )


ix <- which(colnames(cimac_10021_npx_raw) %in% rownames(design))
identical(colnames(cimac_10021_npx_raw[,ix]),rownames(design))

### STATs
vfit <- lmFit(cimac_10021_npx_raw[,ix], design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE) ; summary(decideTests(efit))

### Store data
lmfreq_results <- list()
for ( i in 1:ncol(summary(decideTests(efit))) ) { lmfreq_results[[i]] <- topTable(efit, coef=i,n=Inf,adjust.method="fdr") 
lmfreq_results[[i]]$Comparison <- colnames(summary(decideTests(efit)))[i]
lmfreq_results[[i]]$Marker <- rownames( topTable(efit, coef=i,n=Inf,adjust.method="fdr") ) }
olink_de_results_10021 <- do.call(rbind,lmfreq_results)
rownames(olink_de_results_10021) <- NULL

olink_de_results_10021$nLogFDR <- -log10(olink_de_results_10021$adj.P.Val)
table(olink_de_results_10021$adj.P.Val<0.05)

### barplot(-log10(olink_de_results_10021$adj.P.Val))

### Make Heatmap of results
###################################################################################
pdf(file="/Users/gonzae34/Documents/projects_gnjatic/U24/figures_u24/10021_figures/10021_de_stats_heatmap.pdf",width = 6.5,height = 4)
ggplot(olink_de_results_10021[olink_de_results_10021$adj.P.Val<0.05,], aes(x=Marker, y=Comparison, color=logFC, size=nLogFDR )) + 
  geom_point(shape=16, colour = "black", aes(size = max(nLogFDR))) +
  geom_point(shape=16, colour = "white", aes(size = 0.8*max(nLogFDR))) +
  geom_point(shape=16, aes(size = 0.81*nLogFDR)) + 
  scale_colour_gradient2(low='firebrick', mid="white",high='steelblue', name='LogFC') + #,,limits=c(-1,1), breaks = c(-2.8, 0, 1, 5)) +
  theme_classic2() + rotate_x_text(angle=45) + 
  labs(x ='', y='', title='') 
dev.off()
dev.off()
###################################################################################

### copy varibles for easy handling of the names
cimac_10021_samples$cimac_id <- cimac_10021_samples$`Cimac id`
cimac_10021_samples$QC_Warning <- cimac_10021_samples$`QC Warning`
rownames(cimac_10021_samples) <- cimac_10021_samples$`Cimac id`

### Heatmap
################################################################################
ix <- which(rownames(cimac_10021_npx_raw) %in% olink_de_results_10021$Marker[olink_de_results_10021$adj.P.Val<0.05])

#iy <- c(grep("Baseline",cimac_10021_samples$condition),
#        grep("Cycle_2_Week_3",cimac_10021_samples$condition))

iy <- which(cimac_10021_samples$codition %in% "Other")

pdf(file="/Users/gonzae34/Documents/projects_gnjatic/U24/figures_u24/10021_figures/10021_de_markers_heatmap.pdf", width = 8, height =5)
pheatmap(cimac_10021_npx_raw[ix,],
         #color = colorRampPalette(c("steelblue","white","firebrick"))(255),
         color = viridis(255), 
         #clustering_method = "ward.D2",
         angle_col = 90,
         show_rownames = TRUE, 
         cluster_cols = TRUE,
         show_colnames = FALSE, 
         annotation_col = cimac_10021_samples[,c("condition","cohort","QC_Warning")],
         scale = "row")
dev.off()
dev.off()
################################################################################


### Violin plots
################################################################################
ix <- which(rownames(cimac_10021_npx_raw) %in% olink_de_results_10021$Marker[olink_de_results_10021$adj.P.Val<0.05])
melted_df <- as.data.frame(cimac_10021_npx_raw[ix,])
melted_df$Marker <- rownames(melted_df)
melted_df <- (melt(melted_df))
################################################################################

colnames(melted_df)[2] <- "cimac_id"
melted_df <- merge(melted_df,cimac_10021_samples[,c("condition","QC_Warning","cimac_id")],by="cimac_id")

### Comparisons
my_comparisons <- list( c("Baseline", "Other"),
                        c("Baseline", "On_Treatment"),
                        c("On_Treatment", "Other") )
                     
melted_df$condition <- factor(melted_df$condition, levels=c("Baseline","On_Treatment","Other"))

### boxplot
################################################################################
pdf(file="/Users/gonzae34/Documents/projects_gnjatic/U24/figures_u24/10021_figures/10021_de_markers_boxplot.pdf", width = 20, height = 14)
ggboxplot(data=melted_df, x="condition" , y="value",color='condition') +
  facet_wrap(~Marker,ncol=6, scales = 'free') +
  theme_bw() +  #geom_quasirandom() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  #scale_fill_manual(values = c("firebrick","steelblue")) +
  #scale_color_manual(values = c("firebrick","firebrick","steelblue","steelblue","steelblue")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE ) +
  theme(legend.position = "none") + 
  labs(x ='', y='', title='',color="",fill="") + 
  theme(text=element_text(size=rel(5.5)),strip.text=element_text(size=rel(2.5),face="bold"))
dev.off() 
dev.off() 


### PCA
################################################################################
### PCA
cimac_10021_npx_raw[is.na(cimac_10021_npx_raw)] <- 0
#ix <- grep("Ctrl",rownames(cimac_10021_npx_raw)) ### no patterns ### remove not select from data.
ix <- which(rownames(cimac_10021_npx_raw) %in% olink_de_results_10021$Marker[olink_de_results_10021$adj.P.Val<0.05])
mydata <- prcomp( t(cimac_10021_npx_raw[ix,]) , scale=TRUE, center = TRUE )

### SCREE
Scree.plot <- fviz_eig(mydata, addlabels = T, barfill = "lightsteelblue3", barcolor = "lightsteelblue3") + theme_bw()
eigen <- get_eigenvalue(mydata)
pca1 <- format(round(eigen$variance.percent[1], 2), nsmall = 2)
pca2 <- format(round(eigen$variance.percent[2], 2), nsmall = 2) ;Scree.plot ; rm(eigen)

prcomp_obj_mssm_pca <- as.data.frame(mydata$x)
prcomp_obj_mssm_pca$cimac_id <- rownames(mydata$x)
prcomp_obj_mssm_pca <- merge(prcomp_obj_mssm_pca,cimac_10021_samples,by="cimac_id")

###
pdf(file="/Users/gonzae34/Documents/projects_gnjatic/U24/figures_u24/10021_figures/10021_PCA.pdf",width = 8,height = 6)
ggplot(data=prcomp_obj_mssm_pca, aes(x=PC1, y=PC2, color = condition)) + 
  geom_point() + 
  #geom_label_repel(aes(label=PatientID),force = 10, max.iter = 5000,show.legend = FALSE) +
  theme_classic2() + 
  labs(x=paste("PC1:",pca1,"%",sep=""),title="", y=paste("PC1:", pca2,"%", sep=""), color="Timepoints") +
  theme(text = element_text(size=rel(4.8)),
        legend.text=element_text(size=rel(4.0)),
        plot.title=element_text(size=rel(4.8)) ) + 
  theme(legend.position="bottom") + 
  guides(color=guide_legend(ncol=2), shape = guide_legend(override.aes = list(size = rel(4.0)))) 
dev.off()
###

### STORE Analysis
################################################################################
#save.image(file="/Users/gonzae34/Documents/projects_gnjatic/U24/rdata_files/10021_analysis.RData")
################################################################################