sex.qc <- function(exp.data=Meth.data){
  set.raw <- mapToGenome(preprocessRaw(exp.data))
  set.raw <- addSex(set.raw)

  if(length(which(!colData(set.raw)$`predictedSex` == sample_sheet$Sex[!sample_sheet$Sample_Name %in% excluded])) > 0){
    conflicting_sex.raw <- as.data.frame(merge(x = colData(set.raw),y = sample_sheet[!sample_sheet$Sample_Name %in% excluded,]))[c("Sample_Name","Sex","predictedSex")]
    conflicting_sex.raw <- conflicting_sex.raw[which(conflicting_sex.raw$Sex != conflicting_sex.raw$predictedSex),]
    cat(paste("Conflicting Sex identified in ",nrow(conflicting_sex.raw)," cases (exluding NA values - preprocessRaw)\n",sep = ""))
    write.table(x = conflicting_sex.raw,file = paste("SampleSex_",paste(test_group,collapse = "_"),"_conflicts.txt",sep=""),sep = "\t",append = F,quote = F,row.names = F,col.names = T)
  }
  plot(x = set.raw$xMed, y = set.raw$yMed, type = "n", xlab = "X chr, median total intensity (log2)", ylab = "Y chr, median total intensity (log2)")
  title(main = "No normalisation - Actual (preprocessRaw)",cex.main = 0.75)
  text(x = set.raw$xMed, y = set.raw$yMed, labels = sample_sheet$Sample_Name,
       col = ifelse(sample_sheet$Sex[!sample_sheet$Sample_Name %in% excluded] == "M", "deepskyblue",ifelse(sample_sheet$Sex[!sample_sheet$Sample_Name %in% excluded] == "F", "deeppink3","grey50")))
  legend("bottomleft", c("M", "F","NA"), col = c("deepskyblue","deeppink3","grey50"), pch = 16)

  plot(x = set.raw$xMed, y = set.raw$yMed, type = "n", xlab = "X chr, median total intensity (log2)", ylab = "Y chr, median total intensity (log2)")
  title(main = "No normalisation - Predicted (preprocessRaw)",cex.main = 0.75)
  text(x = set.raw$xMed, y = set.raw$yMed, labels = set.raw$Sample_Name, col = ifelse(set.raw$predictedSex == "M", "deepskyblue", "deeppink3"))
  legend("bottomleft", c("M", "F","NA"), col = c("deepskyblue","deeppink3","grey50"), pch = 16)

}

## Normalisation function
normalise <- function(exp.data=Meth.data){
  if(preprocess.var == "Raw"){
    set <- preprocessRaw(exp.data)
    return(set)
  } else if(preprocess.var == "Quantile"){
    set <- preprocessQuantile(exp.data)
    return(set)
  } else if(preprocess.var == "Funnorm"){
    set <- preprocessFunnorm(exp.data,verbose = T)
    return(set)
  } else if(preprocess.var == "Noob"){
    set <- preprocessNoob(exp.data)
    return(set)
  } else if(preprocess.var == "Illumina"){
    set <- preprocessIllumina(exp.data)
    return(set)
  } else if(preprocess.var == "SWAN"){
    set <- preprocessSWAN(exp.data)
    return(set)
  } else {
    print("[preprocessing] Unknown preprocessing method - exiting now")
  }
}

## DMR annotate function
annotation <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene,annotation="org.Hs.eg.db")
annotate_dmrs <- function(x){
  dmrs.GRange <- makeGRangesFromDataFrame(as.data.frame(x$table),keep.extra.columns = T)
  genes <- matchGenes(x = dmrs.GRange, subject = annotation,skipExons = F,promoterDist = 2500,type = "fiveprime")
  values(dmrs.GRange) <- c(values(dmrs.GRange),genes)
  dmrs.GRange <- dmrs.GRange[dmrs.GRange$distance < dmr_distance,]
  return(dmrs.GRange)
}

## TOPGo gene selection function
geneSelction <- function(x){
  x <- x < 0.05
  return(x)
}

## Probe exclusion function
probe.QC <- function(x){
  probe.p <- detectionP(Meth.data)
  row.probe.p <- rowMeans(probe.p)
  row.keep <- row.probe.p < 0.05
  x <- x[row.keep,]
  if(length(names(row.probe.p[row.probe.p > 0.05])) > 0){
    write.table(x = paste(names(row.probe.p[row.probe.p > 0.05]),"dectionP",sep = "\t"),
                file = paste("excluded_probes_",paste(test_group,collapse = "_"),".txt",sep = ""),
                sep = "\t",
                append = T,
                col.names = F,
                row.names = F,
                quote = F)
  }
  return(x)
}

## Sample exclusion fucntion
sample.QC <- function(x){
  sample.p <- detectionP(Meth.data)
  col.probe.p <- colMeans(sample.p)
  col.keep <- col.probe.p < 0.05
  x <- x[,col.keep]
  if(length(names(col.probe.p[col.probe.p > 0.05])) > 0){
    write.table(x = paste(names(col.probe.p[col.probe.p > 0.05]),"dectionP",sep = "\t"),
                file = paste("excluded_samples_",paste(test_group,collapse = "_"),".txt",sep = ""),
                sep = "\t",
                append = T,
                col.names = F,
                row.names = F,
                quote = F)
  }
  return(x)
}

## dendrogram colours function all
colLab <- function(n){
  if(is.leaf(n)){
    a <- attributes(n)
    group <- pData(gRS.flt)$Group[which(rownames(pData(gRS.flt)) == attributes(n)$label)]
    point_type <- pData(gRS.flt)$Mutation_status[which(rownames(pData(gRS.flt)) == attributes(n)$label)]
    cols <- unique(clrs.clst[names(clrs.clst) %in% group])
    attr(n,"nodePar")<-c(a$nodePar,list(cex=1,lab.cex=1,pch=as.numeric(point_type)*15,col=cols,lab.font=1,lab.cex=1))
  }
  return(n)
}

## dendrogram colours function sub
colLab.sub <- function(n){
  if(is.leaf(n)){
    a <- attributes(n)
    group <- pData(sgRS)$Group[which(rownames(pData(sgRS)) == attributes(n)$label)]
    point_type <- pData(sgRS)$Mutation_status[which(rownames(pData(sgRS)) == attributes(n)$label)]
    cols <- unique(clrs.clst.sub[names(clrs.clst.sub) %in% group])
    attr(n,"nodePar")<-c(a$nodePar,list(cex=1,lab.cex=1,pch=as.numeric(point_type)*15,col=cols,lab.font=1,lab.cex=1))
  }
  return(n)
}
