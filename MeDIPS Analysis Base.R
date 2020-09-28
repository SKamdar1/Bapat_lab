source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
a
y
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
library(biomaRt)
BSgenome="BSgenome.Hsapiens.UCSC.hg19"
chr.select=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
uniq=1e-3
extend=300
shift=0
ws=100
#Parental LNCaP set
LNCAP = MEDIPS.createSet(file="C:/Users/Shivani/Desktop/hmC Project/bam/LNCAP_009_sorted_nodup.bam", BSgenome=BSgenome, extend=extend, shift=shift, window_size=ws, chr.select=chr.select, paired=FALSE)
#LNCaP_abl set
ABL = MEDIPS.createSet(file="C:/Users/Shivani/Desktop/hmC Project/bam/abl_010_sorted_nodup.bam", BSgenome=BSgenome, extend=extend, shift=shift, window_size=ws, chr.select=chr.select, paired=FALSE)

#CpG Density analysis
CS=MEDIPS.couplingVector(pattern="CG", refObj=LNCAP)

#MeDIPS diff meth Analysis
#LNCaP vs Abl
mr.edgeR = MEDIPS.meth(MSet1=LNCAP, MSet2=ABL, CSet=CS, p.adj="bonferroni", diff.method="edgeR", MeDIP=T, CNV=F, minRowSum=10)

#Identify regions with FC > 1.5 or < 0.5
mr.edgeR.s = MEDIPS.selectSig(results=mr.edgeR,p.value=1,adj=T, ratio=1.5,bg.counts=NULL,CNV=F)
mr.edgeR.s.gain = mr.edgeR.s[which(mr.edgeR.s[,(grep("edgeR.logFC",colnames(mr.edgeR.s)))] < 0),]
mr.edgeR.s.gain.m = MEDIPS.mergeFrames(frames=mr.edgeR.s.gain,distance=1)

mr.edgeR.s.loss = mr.edgeR.s[which(mr.edgeR.s[,grep("logFC",colnames(mr.edgeR.s))] > 0),]
mr.edgeR.s.loss.m = MEDIPS.mergeFrames(frames=mr.edgeR.s.loss,distance=1)

#Extract regions of interest (methylation gain/loss) into different table

mycolumns <- names(mr.edgeR)[grep("counts|rpkm|logFC|p.value",names(mr.edgeR))]
rois_gain <- MEDIPS.selectROIs(results=mr.edgeR.s,rois=mr.edgeR.s.gain.m,columns=mycolumns,summarize="avg")
rois_loss <- MEDIPS.selectROIs(results=mr.edgeR.s,rois=mr.edgeR.s.loss.m,columns=mycolumns,summarize="avg")


#Annotate the master file (ENSID) to hg19
anno.mart.gene=MEDIPS.getAnnotation(host="grch37.ensembl.org", dataset=c("hsapiens_gene_ensembl"),annotation=c("GENE"), chr=chr.select)
mr.edgeR.roi.gain = MEDIPS.setAnnotation(regions=rois_gain, annotation=anno.mart.gene)
mr.edgeR.roi.loss = MEDIPS.setAnnotation(regions=rois_loss, annotation=anno.mart.gene)
#Update to HGNC gene symbols
genes.gain <- mr.edgeR.roi.gain$'1_id'
genes.loss <- mr.edgeR.roi.loss$'1_id'
GeneList.gain <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "hgnc_symbol"),values=genes.gain,mart=useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl",host="grch37.ensembl.org"))
GeneList.loss <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "hgnc_symbol"),values=genes.loss,mart=useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl",host="grch37.ensembl.org"))
colnames(GeneList.gain)[1] <- c("1_id")
colnames(GeneList.loss)[1] <- c("1_id")
mr.edgeR.annot.gain <- merge(mr.edgeR.roi.gain, GeneList.gain, by="1_id", all=TRUE)
mr.edgeR.annot.loss <- merge(mr.edgeR.roi.loss, GeneList.loss, by="1_id", all=TRUE)

#Export data

write.table(mr.edgeR.annot.gain, file="LNCAP_abl_hMeDIP_gain.txt", sep="\t", col.names=TRUE)
write.table(mr.edgeR.annot.loss, file="LNCAP_abl_hMeDIP_loss.txt", sep="\t", col.names=TRUE)

########################################################################################
#Annotation

# Installing dependencies
source("https://bioconductor.org/biocLite.R")
if (require(ChIPpeakAnno)==FALSE) {
  biocLite("ChIPpeakAnno")
  a
}
library(ChIPpeakAnno)

if (require(biomaRt)==FALSE) {
  biocLite("biomaRt")
  a
}
library(biomaRt)

if (require(EnsDb.Hsapiens.v75)==FALSE) {
  biocLite("EnsDb.Hsapiens.v75")
  a
}
library(EnsDb.Hsapiens.v75)

#################################################################################################################################################################

#CpG island annotation (performed first so that we can adjust column names accordingly)
#ExtractCpGs
library(AnnotationHub)
hub <- AnnotationHub()
cpgs <- hub[["AH5086"]]

#Format the data to GRanges
CpGChip <- makeGRangesFromDataFrame(mr.edgeR.annot.loss,keep.extra.columns=TRUE,ignore.strand = TRUE, seqinfo=NULL, seqnames.field="chr",start.field="start",end.field="end")
annotatedCpGPeak <- annotatePeakInBatch(CpGChip,AnnotationData=cpgs)
AnnotCpGChip <- as.data.frame(annotatedCpGPeak)

#Format the data to prepare for non-CpG annotation
colnames(AnnotCpGChip)[c(40:length(colnames(AnnotCpGChip)))] <- paste("CpG_",colnames(AnnotCpGChip)[c(40:length(colnames(AnnotCpGChip)))],sep="")
colnames(AnnotCpGChip)[1] <- "chr"
AnnotCpGChip <- AnnotCpGChip[,-c(4:5)]

################################################################################################################################################################

#Regular TSS annotation

#Formatting your data
RangeChip <- makeGRangesFromDataFrame(AnnotCpGChip, keep.extra.columns=TRUE, ignore.strand=TRUE, seqinfo=NULL, seqnames.field="chr",start.field="start", end.field="end")
#The actual annotation
annoData <- toGRanges(EnsDb.Hsapiens.v75, feature="gene")
annotatedPeak <- annotatePeakInBatch(RangeChip, AnnotationData=annoData)
AnnotChip <- as.data.frame(annotatedPeak)

#Output the actual table

write.table(AnnotChip,file="LNCAP_abl_hMeDIP_loss_annot.txt",sep="\t")

########################################################################################
#QC CHECKS

#Saturation analysis

sr_LNCAP <- MEDIPS.saturation(file="C:/Users/Shivani/Desktop/hmC Project/bam/LNCAP_009_sorted_nodup.bam",BSgenome=BSgenome,uniq=uniq,extend=extend,shift=shift,window_size=ws,chr.select=chr.select,nit=10,nrit=1,empty_bins=TRUE,rank=TRUE)
sr_ABL <- MEDIPS.saturation(file="C:/Users/Shivani/Desktop/hmC Project/bam/abl_010_sorted_nodup.bam",BSgenome=BSgenome,uniq=uniq,extend=extend,shift=shift,window_size=ws,chr.select=chr.select,nit=10,nrit=1,empty_bins=TRUE,rank=TRUE)
#Plot saturation
png(filename="LNCAP_hMeDIP_SaturationPlot.png")
MEDIPS.plotSaturation(sr_LNCAP)
dev.off()
png(filename="ABL_hMeDIP_SaturationPlot.png")
MEDIPS.plotSaturation(sr_ABL)
dev.off()

#Coverage analysis
cr_LNCAP <- MEDIPS.seqCoverage(file="C:/Users/Shivani/Desktop/hmC Project/bam/LNCAP_009_sorted_nodup.bam",pattern="CG",BSgenome=BSgenome,uniq=uniq,extend=extend,shift=shift,chr.select=chr.select)
cr_ABL <- MEDIPS.seqCoverage(file="C:/Users/Shivani/Desktop/hmC Project/bam/abl_010_sorted_nodup.bam",pattern="CG",BSgenome=BSgenome,uniq=uniq,extend=extend,shift=shift,chr.select=chr.select)
#Plot coverage analysis
png(filename="LNCAP_hMeDIP_SequenceCoverage.png")
MEDIPS.plotSeqCoverage(seqCoverageObj=cr_LNCAP,type="pie",cov.level=c(0,1,2,3,4,5))
dev.off()
png(filename="ABL_hMeDIP_SequenceCoverage.png")
MEDIPS.plotSeqCoverage(seqCoverageObj=cr_ABL,type="pie",cov.level=c(0,1,2,3,4,5))
dev.off()


########################################################################################

#Plot the CNV alterations - preparation
if (require(ggplot2)==FALSE) {
  install.packages("ggplot2")
}
library(ggplot2)

mr.edgeR.s <- mr.edgeR.s[order(mr.edgeR.s$start),]
CNVPlot <- data.frame(mr.edgeR.s$chr, mr.edgeR.s$stop, mr.edgeR.s$CNV.log2.ratio)
colnames(CNVPlot) <- c("chr", "pos", "CNV")

#Iterate through chromosomes and ggplot/save

myplot <- ggplot(data=CNVPlot,aes(x=pos,y=CNV)) + geom_point(alpha=.1) + theme_bw() + xlab("Chr11 (Bases)") + ylab("log2 CNV Ratio")

pdf("Chromosome11.pdf")
print(myplot)
dev.off()

#Output the actual table
write.table(mr.edgeR.s,file="Chromosome11.txt",sep="\t")

#Then don't forget to do the same thing for RWPE-1 vs CR2!
