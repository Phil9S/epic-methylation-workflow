# ## Install bioc manager
#
# if (!requireNamespace("BiocManager", quietly = TRUE)){
#   install.packages("BiocManager",dependencies = TRUE)
#   library(BiocManager)
# }
# ## CRAN packages
# CRANpackages <- c("reshape2","ggplot2","cowplot",
#                   "circlize","RColorBrewer","ggdendro")
# for(package in CRANpackages){
#   if (!sapply(package,require,character=TRUE)){
#     install.packages(package,dependencies = TRUE)
#     sapply(package,require,character=TRUE)
#
#   }
# }
#
# ## Bioc packages
# bicopackages <- c("IlluminaHumanMethylationEPICmanifest","DMRcate",
#                   "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
#                   "GenomicFeatures","TxDb.Hsapiens.UCSC.hg19.knownGene",
#                   "org.Hs.eg.db","minfi","sva","doParallel","topGO","limma")
# for(package in bicopackages){
#   if (!sapply(package,require,character=TRUE)){
#     BiocManager::install(package)
#     sapply(package,require,character=TRUE)
#   }
# }
#
# ## Functions
# ## Function for performing data QC
# sex.qc <- function(exp.data=Meth.data){
#   set.raw <- mapToGenome(preprocessRaw(exp.data))
#   set.raw <- addSex(set.raw)
#
#   if(length(which(!colData(set.raw)$`predictedSex` == sample_sheet$Sex[!sample_sheet$Sample_Name %in% excluded])) > 0){
#     conflicting_sex.raw <- as.data.frame(merge(x = colData(set.raw),y = sample_sheet[!sample_sheet$Sample_Name %in% excluded,]))[c("Sample_Name","Sex","predictedSex")]
#     conflicting_sex.raw <- conflicting_sex.raw[which(conflicting_sex.raw$Sex != conflicting_sex.raw$predictedSex),]
#     cat(paste("Conflicting Sex identified in ",nrow(conflicting_sex.raw)," cases (exluding NA values - preprocessRaw)\n",sep = ""))
#     write.table(x = conflicting_sex.raw,file = paste("SampleSex_",paste(test_group,collapse = "_"),"_conflicts.tsv",sep=""),sep = "\t",append = F,quote = F,row.names = F,col.names = T)
#   }
#   png(file = paste("methyl_sex_",paste(test_group,collapse = "_"),"_known.png",sep=""),height = 6,width = 8,units = "in",res = 600)
#   plot(x = set.raw$xMed, y = set.raw$yMed, type = "n", xlab = "X chr, median total intensity (log2)", ylab = "Y chr, median total intensity (log2)")
#   title(main = "No normalisation - Actual (preprocessRaw)",cex.main = 0.75)
#   text(x = set.raw$xMed, y = set.raw$yMed, labels = sample_sheet$Sample_Name,
#        col = ifelse(sample_sheet$Sex[!sample_sheet$Sample_Name %in% excluded] == "M", "deepskyblue",ifelse(sample_sheet$Sex[!sample_sheet$Sample_Name %in% excluded] == "F", "deeppink3","grey50")))
#   legend("bottomleft", c("M", "F","NA"), col = c("deepskyblue","deeppink3","grey50"), pch = 16)
#   dev.off()
#
#   png(file = paste("methyl_sex_",paste(test_group,collapse = "_"),"_pred.png",sep=""),height = 6,width = 8,units = "in",res = 600)
#   plot(x = set.raw$xMed, y = set.raw$yMed, type = "n", xlab = "X chr, median total intensity (log2)", ylab = "Y chr, median total intensity (log2)")
#   title(main = "No normalisation - Predicted (preprocessRaw)",cex.main = 0.75)
#   text(x = set.raw$xMed, y = set.raw$yMed, labels = set.raw$Sample_Name, col = ifelse(set.raw$predictedSex == "M", "deepskyblue", "deeppink3"))
#   legend("bottomleft", c("M", "F","NA"), col = c("deepskyblue","deeppink3","grey50"), pch = 16)
#   dev.off()
# }
#
# ## Normalisation function
# normalise <- function(exp.data=Meth.data){
#   if(preprocess.var == "Raw"){
#     set <- preprocessRaw(exp.data)
#     return(set)
#   } else if(preprocess.var == "Quantile"){
#     set <- preprocessQuantile(exp.data)
#     return(set)
#   } else if(preprocess.var == "Funnorm"){
#     set <- preprocessFunnorm(exp.data,verbose = T)
#     return(set)
#   } else if(preprocess.var == "Noob"){
#     set <- preprocessNoob(exp.data)
#     return(set)
#   } else if(preprocess.var == "Illumina"){
#     set <- preprocessIllumina(exp.data)
#     return(set)
#   } else if(preprocess.var == "SWAN"){
#     set <- preprocessSWAN(exp.data)
#     return(set)
#   } else {
#     print("[preprocessing] Unknown preprocessing method - exiting now")
#   }
# }
#
# ## DMR annotate function
# annotation <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene,annotation="org.Hs.eg.db")
# annotate_dmrs <- function(x){
#   dmrs.GRange <- makeGRangesFromDataFrame(as.data.frame(x$table),keep.extra.columns = T)
#   genes <- matchGenes(x = dmrs.GRange, subject = annotation,skipExons = F,promoterDist = 2500,type = "fiveprime")
#   values(dmrs.GRange) <- c(values(dmrs.GRange),genes)
#   dmrs.GRange <- dmrs.GRange[dmrs.GRange$distance < dmr_distance,]
#   return(dmrs.GRange)
# }
#
# ## TOPGo gene selection function
# geneSelction <- function(x){
#   x <- x < 0.05
#   return(x)
# }
#
# ## Probe exclusion function
# probe.QC <- function(x){
#   probe.p <- detectionP(Meth.data)
#   row.probe.p <- rowMeans(probe.p)
#   row.keep <- row.probe.p < 0.05
#   x <- x[row.keep,]
#   if(length(names(row.probe.p[row.probe.p > 0.05])) > 0){
#     write.table(x = paste(names(row.probe.p[row.probe.p > 0.05]),"dectionP",sep = "\t"),
#                 file = paste("excluded_probes_",paste(test_group,collapse = "_"),".txt",sep = ""),
#                 sep = "\t",
#                 append = T,
#                 col.names = F,
#                 row.names = F,
#                 quote = F)
#   }
#   return(x)
# }
#
# ## Sample exclusion fucntion
# sample.QC <- function(x){
#   sample.p <- detectionP(Meth.data)
#   col.probe.p <- colMeans(sample.p)
#   col.keep <- col.probe.p < 0.05
#   x <- x[,col.keep]
#   if(length(names(col.probe.p[col.probe.p > 0.05])) > 0){
#     write.table(x = paste(names(col.probe.p[col.probe.p > 0.05]),"dectionP",sep = "\t"),
#                 file = paste("excluded_samples_",paste(test_group,collapse = "_"),".txt",sep = ""),
#                 sep = "\t",
#                 append = T,
#                 col.names = F,
#                 row.names = F,
#                 quote = F)
#   }
#   return(x)
# }
#
# ## dendrogram colours function all
# colLab <- function(n){
#   if(is.leaf(n)){
#     a <- attributes(n)
#     group <- pData(gRS.flt)$Group[which(rownames(pData(gRS.flt)) == attributes(n)$label)]
#     point_type <- pData(gRS.flt)$Mutation_status[which(rownames(pData(gRS.flt)) == attributes(n)$label)]
#     cols <- unique(clrs.clst[names(clrs.clst) %in% group])
#     attr(n,"nodePar")<-c(a$nodePar,list(cex=1,lab.cex=1,pch=as.numeric(point_type)*15,col=cols,lab.font=1,lab.cex=1))
#   }
#   return(n)
# }
#
# ## dendrogram colours function sub
# colLab.sub <- function(n){
#   if(is.leaf(n)){
#     a <- attributes(n)
#     group <- pData(sgRS)$Group[which(rownames(pData(sgRS)) == attributes(n)$label)]
#     point_type <- pData(sgRS)$Mutation_status[which(rownames(pData(sgRS)) == attributes(n)$label)]
#     cols <- unique(clrs.clst.sub[names(clrs.clst.sub) %in% group])
#     attr(n,"nodePar")<-c(a$nodePar,list(cex=1,lab.cex=1,pch=as.numeric(point_type)*15,col=cols,lab.font=1,lab.cex=1))
#   }
#   return(n)
# }
#
# ## Vars
# preprocess.var <- "Funnorm"
# samplesheet <- "MethArray_Unknown_Unknown/850K_final_sample_sheet_recoded_NOV2019.txt"
# excluded <- c("D1234004","s065669","S164733") # ALWAYS EXCLUDE
# excluded <- c(excluded,"TB17.0821","gtb81831") # ADD ADDITIONAL SAMPLES HERE
# test_group <- c("PPGL","GLIOMA")
# cores <- 4
# bmphunter_cutoff <- 0.3
# bmphunter_permutation <- 1000
# test_pheno <- "Mutation_status"
# dmr_distance <- 2000
# SVA  <- TRUE
#
# ## Report excluded samples
# write.table(x = paste(excluded,"MANUAL",sep = "\t"),
#             file = paste("excluded_samples_",paste(test_group,collapse = "_"),".txt",sep = ""),
#             sep = "\t",
#             append = T,
#             col.names = F,
#             row.names = F,
#             quote = F)
#
# ## Load sample sheet & exclude known unwanted samples
# sample_sheet <- read.table(samplesheet,header = T,stringsAsFactors = F,sep = "\t",strip.white = T)
# write.table(x = paste(sample_sheet$Sample_Name[sample_sheet$Exclude == TRUE],"SAMPLE_SHEET_ERRORS",sep = "\t"),
#             file = paste("excluded_samples_",paste(test_group,collapse = "_"),".txt",sep = ""),
#             sep = "\t",
#             append = T,
#             col.names = F,
#             row.names = F,
#             quote = F)
# sample_sheet <- sample_sheet[sample_sheet$Exclude == FALSE,]
#
# ## Manifest ann
# annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
#
# ## Read methylation data & label
# Meth.data <- read.metharray.exp(targets = sample_sheet)
# sampleNames(Meth.data) <- as.character(sample_sheet$Sample_Name)
#
# ## Print manifest
# # getManifest(Meth.data)
#
# ## Manual exclusion of specific samples
# Meth.data <- Meth.data[,!colnames(Meth.data) %in% excluded]
#
# ## Sex QC
# sex.qc(Meth.data)
#
# ## Failed probe / sample QC
# Meth.data.QC <- probe.QC(Meth.data)
# Meth.data.QC <- sample.QC(Meth.data.QC)
#
# ## Perform normalisation method
# gRS <- normalise(exp.data = Meth.data.QC)
#
# ## Drop chrX and chrY probes - mixed sex subset
# write.table(x = paste(rownames(gRS[featureNames(gRS) %in% annEPIC$Name[annEPIC$chr %in% c("chrX","chrY")],]),"CHR_X/Y",sep = "\t"),
#             file = paste("excluded_probes_",paste(test_group,collapse = "_"),".txt",sep = ""),
#             sep = "\t",
#             append = T,
#             col.names = F,
#             row.names = F,
#             quote = F)
# gRS.flt <- gRS[!featureNames(gRS) %in% annEPIC$Name[annEPIC$chr %in% c("chrX","chrY")],]
#
# ## Filter probes overlaping with common SNPs
# write.table(x = paste(rownames(gRS.flt)[!rownames(gRS.flt) %in% rownames(dropLociWithSnps(gRS.flt, snps=c("SBE","CpG"), maf=0.05))],"SNP_MAF_5pct",sep = "\t"),
#             file = paste("excluded_probes_",paste(test_group,collapse = "_"),".txt",sep = ""),
#             sep = "\t",
#             append = T,
#             col.names = F,
#             row.names = F,
#             quote = F)
# gRS.flt <- dropLociWithSnps(gRS.flt, snps=c("SBE","CpG"), maf=0.05)
#
# ## Cross-reactivity data - as provided here https://doi.org/10.1186/s13059-016-1066-1
# x.reactive <- read.table("13059_2016_1066_MOESM1_ESM.csv",sep = ",",header = T,stringsAsFactors = F)
#
# ## Dropping cross reactive probes
# write.table(x = paste(rownames(gRS.flt[featureNames(gRS.flt) %in% x.reactive$Name,]),"cross_reactive",sep = "\t"),
#             file = paste("excluded_probes_",paste(test_group,collapse = "_"),".txt",sep = ""),
#             sep = "\t",
#             append = T,
#             col.names = F,
#             row.names = F,
#             quote = F)
# gRS.flt <- gRS.flt[!featureNames(gRS.flt) %in% x.reactive$Name,]
#
# ###############################################################################
#
# ## General plots
# col.probe.p <- colMeans(detectionP(Meth.data))
# png(paste("probe_sample_",paste(test_group,collapse = "_"),"_excluded.png",sep = ""),width = 16,height = 6,units = "in",res = 600)
# ggplot() + geom_col(aes(y = col.probe.p,x = names(col.probe.p)),color = "grey15",fill="cadetblue3") +
#   geom_hline(yintercept = 0.05,linetype=2) +
#   geom_hline(yintercept = 0.01,linetype=3) +
#   labs(title = "Methylation Detection significance - Cummulative per sample") +
#   xlab("Sample") +
#   ylab("ColMean(detectP-function)") +
#   scale_y_continuous(expand = c(0,0),limits = c(0,0.25)) +
#   theme(axis.text.x = element_text(angle = 90))
# dev.off()
#
# ## Density Plots -before and after normalisation and QC
# png(file = paste("beta_value_",paste(test_group,collapse = "_"),"_normalisation.png",sep = ""),height = 12,width = 6,units = "in",res = 600)
# par(mfrow=c(3,1))
# densityPlot(Meth.data.QC,main = "No normalisation")
# densityPlot(getBeta(gRS),main = "preprocessFunnorm normalisation")
# densityPlot(getBeta(gRS.flt),main = "Post-normalised filtering")
# dev.off()

# ##Clustering all data
# clrs.clst <- brewer.pal(nlevels(as.factor(pData(gRS.flt)$Group)),"Dark2")[as.factor(pData(gRS.flt)$Group)]
# names(clrs.clst) <- as.factor(pData(gRS.flt)$Group)
#
# pca.clst <- t(getBeta(gRS.flt))
# pca.clst.prc <- as.data.frame(prcomp(pca.clst,scale. = T, center = T)$x)
# histCLUST <- as.dendrogram(hclust(dist(scale(pca.clst),method = "euclidean",diag = F),method = "ward.D2"))
# histCLUST <- dendrapply(histCLUST,colLab)
#
# plot(histCLUST, type = "triangle", ylab = "")
# legend("topright",
#        legend = unique(names(clrs.clst)),
#        col = unique(clrs.clst),
#        pch = 15, bty = "n",  pt.cex = 1, cex = 0.5,
#        text.col = "black", horiz = FALSE)
# histPlot <- recordPlot()
#
# pcaPlot <- ggplot(as.data.frame(pca.clst.prc)) +
#   geom_point(aes(PC1,PC2,colour = as.factor(pData(gRS.flt)$Group),shape = as.factor(pData(gRS.flt)$Mutation_status))) +
#   theme_minimal() +
#   theme(legend.direction = "horizontal",legend.position = "bottom",axis.line = element_line()) +
#   scale_color_discrete(name = "Group") +
#   scale_shape_discrete(name = "Mutation status", labels = c("WT","MUT"))
# png(filename = paste("beta_PCAclustering_",paste(test_group,collapse = "_"),"_fltNorm.png",sep = ""),width = 8,height = 8,units = "in",res = 800)
# pcaPlot
# dev.off()
# png(filename = paste("beta_Histclustering_",paste(test_group,collapse = "_"),"_fltNorm.png",sep = ""),width = 12,height = 6,units = "in",res = 800)
# histPlot
# dev.off()

# ## Summary plot and stats
# sample_sheet[which(!sample_sheet$Group %in% c("NORMAL_GASTRIC","NORMAL_BOWEL","NORMAL_ADRENAL")),]
# sample_sheet2$split <- paste(sample_sheet2$Group,ifelse(sample_sheet2$Mutation_status == 1,"mutated","non-mutated"),sep = " ")
# splits <- split(x = sample_sheet2[,1],f = factor(sample_sheet2$split))
# meanBetaPerGroup <- lapply(X = splits,FUN = function(x){colMeans(getBeta(gRS.flt[,colnames(gRS.flt) %in% x]),na.rm = T)})

## Parametric
# summary(aov(formula = value ~ L1,data = melt(meanBetaPerGroup)))
#
# png("parametric_testing_normality_homogeneity.png",width = 8,height = 4,res = 800,units = "in")
# layout(cbind(1,2))
# plot(aov(formula = value ~ L1,data = melt(meanBetaPerGroup)),1)
# mtext(text = paste("Levene test: ",signif(leveneTest(y = value ~ as.factor(L1),data = melt(meanBetaPerGroup))$`Pr(>F)`[1],digits = 2),sep = ""),side = 4)
# plot(aov(formula = value ~ L1,data = melt(meanBetaPerGroup)),2)
# mtext(text = paste("Shapiro-wilk test: ",signif(shapiro.test(residuals(aov(formula = value ~ L1,data = melt(meanBetaPerGroup))))$p.value,digits = 2),sep = ""),side = 4)
# dev.off()
#
# lapply(split(melt(meanBetaPerGroup)$value,f = as.factor(melt(meanBetaPerGroup)$L1)),FUN = function(x) ks.test(x = x,y = pnorm,mean = mean(x),sd = sd(x))$p.value)
# leveneTest(y = value ~ as.factor(L1),data = melt(meanBetaPerGroup))
# shapiro.test(residuals(aov(formula = value ~ L1,data = melt(meanBetaPerGroup))))
# pairwise.t.test(melt(meanBetaPerGroup)$value,melt(meanBetaPerGroup)$L1,p.adjust = "none")

# ## Non-parametric
# kwt <- kruskal.test(value ~ L1,data = melt(meanBetaPerGroup))
# wilcox <- pairwise.wilcox.test(melt(meanBetaPerGroup)$value,melt(meanBetaPerGroup)$L1,p.adjust = "none")
#
# png("GroupBeta_boxplot_wilcox.png",width = 11,height = 6,units = "in",res = 800)
# ggplot(melt(meanBetaPerGroup)) +
#   geom_point(aes(y=value,x=L1,color = L1),position="jitter") +
#   geom_boxplot(aes(x=L1,y=value,fill= L1)) +
#   geom_signif(aes(x=L1,y=value),
#               comparisons = list(c("GLIOMA mutated","GLIOMA non-mutated"),
#                                  c("GIST mutated","GIST non-mutated"),
#                                  c("PPGL mutated","PPGL non-mutated")),
#               test = "wilcox.test",y_position = 0.7,map_signif_level = F) +
#   #geom_signif(aes(x=L1,y=value),comparisons = list(c("GLIOMA mutated","GIST non-mutated")),test = "wilcox.test",y_position = 0.725) +
#   #geom_signif(aes(x=L1,y=value),comparisons = list(c("GLIOMA mutated","GIST mutated")),test = "wilcox.test",y_position = 0.75) +
#   #geom_signif(aes(x=L1,y=value),comparisons = list(c("GLIOMA mutated","PPGL mutated")),test = "wilcox.test",y_position = 0.775) +
#   #geom_signif(aes(x=L1,y=value),comparisons = list(c("GLIOMA mutated","PPGL non-mutated")),test = "wilcox.test",y_position = 0.8) +
#   #geom_signif(aes(x=L1,y=value),comparisons = list(c("GLIOMA non-mutated","PPGL non-mutated")),test = "wilcox.test",y_position = 0.825) +
#   scale_y_continuous(limits = c(NA,0.85)) +
#   scale_fill_manual(name="Group",
#                     values = alpha(c("cadetblue2","cadetblue4","darkolivegreen1","darkolivegreen3","gold1","gold3"),0.7),
#                     aesthetics = c("color","fill")) +
#   ylab("Mean β value") +
#   xlab(NULL) +
#   labs(title = "Mean methylation per group") +
#   theme(legend.position = "bottom",
#         legend.direction = "horizontal")
# dev.off()
# png("GroupBeta_density.png",width = 10,height = 6,units = "in",res = 800)
# ggplot(melt(meanBetaPerGroup)) +
#   geom_density(aes(value, fill=L1)) +
#   #geom_histogram(aes(value, fill=L1),color ="grey15",binwidth = 0.05) +
#   scale_fill_manual(name="Group",
#                     values = alpha(c("cadetblue2","cadetblue4","darkolivegreen1","darkolivegreen3","gold1","gold3"),0.7)) +
#   scale_x_continuous(limits =  c(0,1)) +
#   ylab("Density") +
#   xlab(NULL) +
#   theme(legend.position = "bottom",
#         legend.direction = "horizontal") +
#   facet_wrap(.~L1,nrow = 3,ncol=2)
# dev.off()

# #Save fitlered methyl data for downstream analysis
# save(gRS.flt,file = "filtered.methyl.data.RData")
#
# ## Subset data by samples of interest defined by test_group var
sgRS <- gRS.flt[,pData(gRS.flt)$Group %in% test_group]

## MODEL
pheno <- pData(sgRS)[,colnames(pData(sgRS)) == test_pheno]
design <- model.matrix(~ pheno)

## DMP
beta <- getBeta(sgRS)
dmp <- dmpFinder(beta, pheno = pheno,type = "categorical")
write.table(x = dmp[dmp$pval < 0.05,],file = paste("DMP_table_",paste(test_group,collapse = "_"),".tsv",sep = ""),quote = F,sep = "\t",row.names = T,col.names = T)

## SVA implementation

if(SVA == TRUE){
  mval <- getM(sgRS)
  mval <- mval[-which(!is.finite(rowSums(mval))),]
  pheno.sva <- pData(sgRS)

  mod <- model.matrix(~as.factor(Mutation_status), data=pheno.sva)
  mod0 <- model.matrix(~1, data=pheno.sva)

  sva.results <- sva(mval, mod, mod0)
  modSv = cbind(mod,sva.results$sv)
  fit <- lmFit(mval,modSv)
  write.table(x = fit$design,file = paste("sva_bmpHunter_",paste(test_group,collapse = "_"),"_model.txt",sep = ""),quote = F,sep = "\t",row.names = T,col.names = T)

  registerDoParallel(cores = cores)
  dmrs <- bumphunter(sgRS, design = fit$design,cutoff = bmphunter_cutoff, B=bmphunter_permutation, type="Beta",nullMethod="bootstrap")
} else {
  ## CHECK FOR NUMBER OF BUMPS AND ADJUST CUTOFF according - No SVA implementation
  ## bumphunter
  registerDoParallel(cores = cores)
  #dmrs.len <- bumphunter(sgRS, design = design,cutoff = bmphunter_cutoff, B=0, type="Beta")
  #nrow(dmrs.len$table)
  dmrs <- bumphunter(sgRS, design = design,cutoff = bmphunter_cutoff, B=bmphunter_permutation, type="Beta")
}
## DMR annotation
dmrs.GRange <- annotate_dmrs(x = dmrs)
write.table(x = as.data.frame(dmrs.GRange),file = paste("dmrs_table_",paste(test_group,collapse = "_"),".tsv",sep = ""),append = F,quote = F,sep = "\t",row.names = F,col.names = T)

## GO analysis
GO.ids <- mapIds(x = org.Hs.eg.db,keys = as.character(dmrs.GRange$name),keytype = "SYMBOL",column = "GO",multiVals = "list")
geneL <- dmrs.GRange$fwer
names(geneL) <- as.character(dmrs.GRange$name)
GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneL,
              geneSel = geneSelction,
              annot = annFUN.gene2GO,
              gene2GO = GO.ids,
              nodeSize = 5)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
resultKS.weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "ks")

allRes <- GenTable(GOdata,Fisher = resultFisher,weight01KS = resultKS.weight01,KS = resultKS,KSelim = resultKS.elim,topNodes = 500,numChar = 500)
write.table(x = as.data.frame(allRes),file = paste("GO_table_",paste(test_group,collapse = "_"),".tsv",sep = ""),append = F,quote = F,sep = "\t",row.names = F,col.names = T)

## Sub plots
## Clustering subset
clrs.clst.sub <- brewer.pal(nlevels(as.factor(pData(sgRS)$Group)),"Dark2")[as.factor(pData(sgRS)$Group)]
names(clrs.clst.sub) <- as.factor(pData(sgRS)$Group)

pca.clst.sub <- t(getBeta(sgRS))
pca.clst.prc.sub <- as.data.frame(prcomp(pca.clst.sub,scale. = T, center = T)$x)
histCLUST.sub <- as.dendrogram(hclust(dist(scale(pca.clst.sub),method = "euclidean",diag = F),method = "ward.D2"))
histCLUST.sub <- dendrapply(histCLUST.sub,colLab.sub)

plot(histCLUST.sub, type = "triangle", ylab = "")
legend("topright",
       legend = unique(names(clrs.clst.sub)),
       col = unique(clrs.clst.sub),
       pch = 15, bty = "n",  pt.cex = 1, cex = 0.5,
       text.col = "black", horiz = FALSE)
histPlot.sub <- recordPlot()

pcaPlot.sub  <- ggplot(as.data.frame(pca.clst.prc.sub)) +
  geom_point(aes(PC1,PC2,colour = as.factor(pData(sgRS)$Mutation_status),shape = as.factor(pData(sgRS)$Group))) +
  geom_vline(xintercept = 348.6543*2) +
  geom_vline(xintercept = -348.6543*2) +
  geom_hline(yintercept = 298.865*2) +
  geom_hline(yintercept = -298.865*2) +
  theme_minimal() +
  theme(legend.direction = "horizontal",legend.position = "bottom",axis.line = element_line()) +
  scale_shape_discrete(name = "Group") +
  scale_color_discrete(name = "Mutation status", labels = c("WT","MUT"))

png(filename = paste("beta_PCAclustering_sub_",paste(test_group,collapse = "_"),"_fltNorm.png",sep = ""),width = 8,height = 8,units = "in",res = 800)
pcaPlot.sub
dev.off()
png(filename = paste("beta_Histclustering_sub_",paste(test_group,collapse = "_"),"_fltNorm.png",sep = ""),width = 12,height = 6,units = "in",res = 800)
histPlot.sub
dev.off()

##
hist_max <- max(hist((as.data.frame(dmrs.GRange)$value), breaks=seq(-1, 1, by=0.0075), plot=FALSE)$counts)
meth_hist <- ggplot(as.data.frame(values(dmrs.GRange))) +
                geom_rect(aes(xmin = -1, xmax = -(bmphunter_cutoff), ymin = 0,ymax = hist_max*1.1), fill = "lightblue") +
                geom_rect(aes(xmin = bmphunter_cutoff, xmax = 1, ymin = 0,ymax = hist_max*1.1), fill = "pink") +
                geom_histogram(data = base::subset(as.data.frame(values(dmrs.GRange)),value < 0),aes(value),binwidth = 0.0075,fill = "darkblue") +
                geom_histogram(data = base::subset(as.data.frame(values(dmrs.GRange)),value > 0),aes(value),binwidth = 0.0075,fill = "darkred") +
                geom_vline(xintercept = bmphunter_cutoff,linetype=2,color="grey15")  +
                geom_vline(xintercept = -(bmphunter_cutoff),linetype=2,color="grey15") +
                xlab("Relative methylation change per DMR") +
                ylab("Frequency") +
                scale_x_continuous(limits = c(-1,1),expand = c(0,0)) +
                scale_y_continuous(expand = c(0,0),limits = c(0,hist_max*1.1)) +
                annotate(geom = "text",x = 0,y = hist_max,label="Bumphunter cutoff value",fontface = "bold") +
                annotate(geom = "text",x = ((1 + bmphunter_cutoff) / 2),y = hist_max,label="Increased methylation",fontface = "bold") +
                annotate(geom = "text",x = ((-1 + -bmphunter_cutoff) / 2),y = hist_max,label="Decreased methylation",fontface = "bold") +
                annotate(geom = "line",y = hist_max*0.95,x = c(-0.1,0.1),linetype = 2,color="grey15")

png(filename = paste("meth_hist_dist_",paste(test_group,collapse = "_"),"_fltNorm.png",sep = ""),width = 8,height = 8,units = "in",res = 800)
meth_hist
dev.off()

##
dmr_density <- ggplot(as.data.frame(values(dmrs.GRange))) +
                geom_hex(aes(x=fwer,y=distance),bins=50,color="grey15") +
                scale_fill_gradient(low = "grey90",high = "firebrick1",trans = "log10",name="log10(count)") +
                ylab("Distance (bp)") +
                theme_minimal() +
                theme(axis.line = element_line(colour = "grey5"))

png(filename = paste("fwer_distribution_distance_",paste(test_group,collapse = "_"),"_fltNorm.png",sep = ""),width = 8,height = 8,units = "in",res = 600)
dmr_density
dev.off()

##
den1 <- ggplot(as.data.frame(values(dmrs.GRange))) +
          geom_density(aes(x = p.value)) +
          theme_minimal() +
          theme(axis.line = element_line(colour = "grey5"))
den2 <- ggplot(as.data.frame(values(dmrs.GRange))) +
          geom_density(aes(x = fwer)) +
          theme_minimal() +
          theme(axis.line = element_line(colour = "grey5"))
den3 <- ggplot(as.data.frame(values(dmrs.GRange))) +
          geom_density(aes(x = p.valueArea)) +
          theme_minimal() +
          theme(axis.line = element_line(colour = "grey5"))
den4 <- ggplot(as.data.frame(values(dmrs.GRange))) +
          geom_density(aes(x = fwerArea),colour="red") +
          theme_minimal() +
          theme(axis.line = element_line(colour = "grey5"))

png(filename = paste("signif_value_dists_",paste(test_group,collapse = "_"),".png",sep=""),width = 8,height = 8,units = "in",res = 800)
plot_grid(den1,den2,den3,den4)
dev.off()

## circos plot
bed_list = list(a=as.data.frame(dmrs.GRange[dmrs.GRange$value > 0 & dmrs.GRange$fwer < 0.05]),
                b=as.data.frame(dmrs.GRange[dmrs.GRange$value < 0 & dmrs.GRange$fwer < 0.05]))
circos.clear()
circos.par("gap.degree" = rep(2,24),start.degree = 90)
circos.initializeWithIdeogram(species = "hg19",
                              chromosome.index = c(paste(rep("chr",22),seq.int(1,22,1),sep = ""),"chrX","chrY"),
                              ideogram.height =  .075)
#circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col = c("#FF000080", "#0000FF80"))
circos.genomicDensity(bed_list$a, col = c("#FF000080"), track.height = 0.1)
circos.genomicDensity(bed_list$b, col = c("#0000FF80"), track.height = 0.1)
circosPlot <- recordPlot()

png(filename = paste("circos_dmr_",paste(test_group,collapse = "_"),".png",sep=""),width = 8,height = 8,units = "in",res = 800)
circosPlot
dev.off()

#### From nmf_methylation script

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
