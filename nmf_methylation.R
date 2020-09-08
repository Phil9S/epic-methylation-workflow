# https://bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.html#8_identifying_differentially_methylated_regions
# https://academic.oup.com/bioinformatics/article/33/4/558/2666344
# https://www.ncbi.nlm.nih.gov/pubmed/29469594
# https://www.bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html#identifying-dmrs-and-dmps
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6366976/
# https://www.ncbi.nlm.nih.gov/pubmed/23314698?dopt=Abstract
# http://www.epicenteredresearch.com/amlbinder/SBE2015/Day2/RMarkdown_for_RA_final.html
rm(list=ls())

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  library(BiocManager)
}
if (!require("reshape2", quietly = TRUE)){
  install.packages("reshape2")
  library(reshape2)
}

if (!require(IlluminaHumanMethylationEPICmanifest, quietly = TRUE)){
  BiocManager::install("IlluminaHumanMethylationEPICmanifest", version = "3.8")
  library(IlluminaHumanMethylationEPICmanifest)
}

if (!require(DMRcate, quietly = TRUE)){
  BiocManager::install("DMRcate")
  library(DMRcate)
}

if (!require(IlluminaHumanMethylationEPICanno.ilm10b4.hg19, quietly = TRUE)){
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", version = "3.8")
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
}

if (!require(GenomicFeatures, quietly = TRUE)){
  BiocManager::install("GenomicFeatures")
  library(GenomicFeatures)
}

if (!require(TxDb.Hsapiens.UCSC.hg19.knownGene, quietly = TRUE)){
  BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
}

if (!require(org.Hs.eg.db, quietly = TRUE)){
  BiocManager::install("org.Hs.eg.db")
  library(org.Hs.eg.db)
}

if(!require(minfi,quietly = TRUE)){
  BiocManager::install("minfi")
  library(minfi)
}


if(!require(NMF,quietly = TRUE)){
  BiocManager::install("NMF")
  library(NMF)
}

if (!require(topGO, quietly = TRUE)){
  BiocManager::install("topGO")
  library(topGO)
}

if (!require(limma, quietly = TRUE)){
  BiocManager::install("limma")
  library(limma)
}

if (!require(ggplot2, quietly = TRUE)){
  install.packages("ggplot2")
  library(ggplot2)
}

if (!require(cowplot, quietly = TRUE)){
  install.packages("cowplot")
  library(cowplot)
}

if (!require(circlize, quietly = TRUE)){
  install.packages("circlize")
  library(circlize)
}

if (!require(RColorBrewer, quietly = TRUE)){
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

if (!require(ggdendro, quietly = TRUE)){
  install.packages("ggdendro")
  library(ggdendro)
}

if (!require(RColorBrewer, quietly = TRUE)){
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

if (!require(gplots, quietly = TRUE)){
  install.packages("gplots")
  library(gplots)
}

load("filtered.methyl.data.RData")
pheno_data <- as.data.frame(pData(gRS.flt))

group <- c("PPGL","GLIOMA")

## correlation/Distribution - metabolomic - methylation

metabolomic <- sample_sheet[which(colnames(sample_sheet) %in% c("Sample_Name","X2HG_concentration_uM","Succinate_concentration_uM","Group","Mutation_status"))]
mean_methyl <- data.frame(mean_methyl=colMeans(getBeta(gRS.flt),na.rm = T))
meta.methyl.merge <- merge(mean_methyl,metabolomic,by.x=0,by.y="Sample_Name",all.y=T)

meth.dist <- ggplot(meta.methyl.merge) +
  geom_histogram(aes(mean_methyl),bins = 10,fill="springgreen",color="grey15") +
  xlab(label = "mean β value") +
  theme_linedraw() +
  theme(axis.title.y = element_blank(),plot.margin = unit(x = c(0.01,0.1,0.01,0.01),units = "npc"))

HG.dist <- ggplot(meta.methyl.merge) +
  geom_histogram(aes(X2HG_concentration_uM),bins = 10,fill="tomato",color="grey15") +
  xlab(label = "2HG concentration (µM)") +
  theme_linedraw() +
  theme(axis.title.y = element_blank(),plot.margin = unit(x = c(0.01,0.1,0.01,0.01),units = "npc"))
  
succ.dist <- ggplot(meta.methyl.merge) +
  geom_histogram(aes(Succinate_concentration_uM),bins = 5,fill="skyblue",color="grey15") +
  xlab(label = "Succinate concentration (µM)") +
  theme_linedraw() +
  theme(axis.title.y = element_blank(),plot.margin = unit(x = c(0.01,0.1,0.01,0.01),units = "npc"))  

X2HG.data <- meta.methyl.merge[meta.methyl.merge$Group %in% "GLIOMA" & !is.na(meta.methyl.merge$X2HG_concentration_uM) & !is.na(meta.methyl.merge$mean_methyl),]
Succinate.data <- meta.methyl.merge[meta.methyl.merge$Group == "PPGL" & !is.na(meta.methyl.merge$Succinate_concentration_uM) & !is.na(meta.methyl.merge$mean_methyl),]

HG <- ggplot(X2HG.data) +
  geom_point(aes(x = mean_methyl,y = X2HG_concentration_uM,color=as.factor(Mutation_status))) +
  geom_smooth(aes(mean_methyl,X2HG_concentration_uM),method = "lm",se = F) +
  xlab("mean β value") + ylab("2HG concentration (µM)") +
  theme_linedraw() +
  theme(legend.position = "none",plot.margin = unit(x = c(0.01,0.1,0.01,0.01),units = "npc"))
HG <- HG + annotate(geom = "text",x = 0.57,y = 5200,hjust = 0,label = paste("p.value =",
                                          signif(cor.test(x = X2HG.data$mean_methyl,y = X2HG.data$X2HG_concentration_uM,method = "spearman")$p.value,digits = 3))) +
     annotate(geom = "text",x = 0.57,y = 5000,hjust = 0,label = paste("spearman rho =",
                                          signif(cor.test(x = X2HG.data$mean_methyl,y = X2HG.data$X2HG_concentration_uM,method = "spearman")$estimate,digits = 3)))

succ <- ggplot(Succinate.data) +
  geom_point(aes(x = mean_methyl,y = Succinate_concentration_uM,color=as.factor(Mutation_status))) +
  geom_smooth(aes(mean_methyl,Succinate_concentration_uM),method = "lm",se = F) +
  xlab("mean β value") + ylab("Succinate concentration (µM)") +
  theme_linedraw() +
  theme(legend.position = "none",plot.margin = unit(x = c(0.01,0.1,0.01,0.01),units = "npc"))
succ <- succ + annotate(geom = "text",x = 0.45,y = 7000,hjust = 0,label = paste("p.value =",
                                          signif(cor.test(x = Succinate.data$mean_methyl,y = Succinate.data$Succinate_concentration_uM,method = "spearman")$p.value,digits = 3))) +
     annotate(geom = "text",x = 0.45,y = 6800,hjust = 0,label = paste("spearman rho =",
                                          signif(cor.test(x = Succinate.data$mean_methyl,y = Succinate.data$Succinate_concentration_uM,method = "spearman")$estimate,digits = 3)))

plot_grid(plotlist = list(plot_grid(meth.dist,succ.dist,HG.dist,ncol = 3,labels = c("A","B","C")),plot_grid(HG,succ,labels = c("D","E"))),nrow = 2)

## NMF

methyl.mat <- as.matrix(getBeta(gRS.flt[,which(colnames(gRS.flt) %in% pheno_data$Sample_Name[pheno_data$Group %in% group])]))
sds <- apply(methyl.mat,1,FUN = function(x) sd(x,na.rm = T))
most.var.rows <- names(sds[sds >= quantile(sds,probs = seq.int(0,1,0.01))[which(names(quantile(sds,probs = seq.int(0,1,0.01))) == "95%")]])
sgRS.flt <- methyl.mat[rownames(methyl.mat) %in% most.var.rows,]

results.nmf <- nmf(sgRS.flt,2:6,nrun=100,.options="t")
rand.RS <- randomize(sgRS.flt)
rand.nmf <- nmf(rand.RS,2:6,nrun=100,.options="t")

pca.data.sig <- prcomp(t(sgRS.flt))
pca.data <- prcomp(t(methyl.mat),center = TRUE,scale. = FALSE)

par(mfrow=c(1,2))
plot(pca.data$x[,1],pca.data$x[,2],
     col=as.factor(pheno_data$Group[pheno_data$Group %in% group]),
     pch=ifelse(as.factor(pheno_data$Mutation_status[pheno_data$Group %in% group]) == 0,15,17))
plot(pca.data.sig$x[,1],pca.data.sig$x[,2],
     col=as.factor(pheno_data$Group[pheno_data$Group %in% group]),
     pch=ifelse(as.factor(pheno_data$Mutation_status[pheno_data$Group %in% group]) == 0,15,17))
dev.off()


plot(results.nmf,rand.nmf)

png(filename = "nmf_concensus.png",width = 8,height = 6,units = "in",res = 600)
consensusmap(results.nmf$fit$`3`,
             annCol=list(Group=pheno_data$Group[pheno_data$Group %in% group],
                         Mutation=pheno_data$Mutation_status[pheno_data$Group %in% group],
                         Sex=pheno_data$Sex[pheno_data$Group %in% group],
                         Tissue=pheno_data$Tissue[pheno_data$Group %in% group]))
dev.off()
