source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
library("MEDIPSData")
library(biomaRt)
BSgenome="BSgenome.Hsapiens.UCSC.hg19"
chr.select=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
uniq=1
extend=300
shift=0
ws=100
#RWPE-1 MEDIPS sets (3 replicates)
RW_MeDIP = MEDIPS.createSet(file="M:/shared/Hydroxymethylation Project/Raw BAM Files - MACS Peaks/RWPE-1/mC/RWPE1_rp1CoordSrt.bam", BSgenome=BSgenome, extend=extend, shift=shift, window_size=ws, chr.select=chr.select, paired=TRUE)
RW_MeDIP = c(RW_MeDIP, MEDIPS.createSet(file="M:/shared/Hydroxymethylation Project/Raw BAM Files - MACS Peaks/RWPE-1/mC/RWPE1_rp2CoordSrt.bam", BSgenome=BSgenome, extend=extend, shift=shift, window_size=ws, chr.select=chr.select, paired=TRUE))
RW_MeDIP = c(RW_MeDIP, MEDIPS.createSet(file="M:/shared/Hydroxymethylation Project/Raw BAM Files - MACS Peaks/RWPE-1/mC/RWPE1_rp3CoordSrt.bam", BSgenome=BSgenome, extend=extend, shift=shift, window_size=ws, chr.select=chr.select, paired=TRUE))
#CR1 MEDIPS sets (2 replicates)
Rv_MeDIP = MEDIPS.createSet(file="M:/shared/Hydroxymethylation Project/Raw BAM Files - MACS Peaks/22Rv1/mC/22RV1_rp1CoordSrt.bam", BSgenome=BSgenome, extend=extend, shift=shift, window_size=ws, chr.select=chr.select, paired=TRUE)
Rv_MeDIP = c(Rv_MeDIP, MEDIPS.createSet(file="M:/shared/Hydroxymethylation Project/Raw BAM Files - MACS Peaks/22Rv1/mC/22RV1_rp2CoordSrt.bam", BSgenome=BSgenome, extend=extend, shift=shift, window_size=ws, chr.select=chr.select, paired=TRUE))
Rv_MeDIP = c(Rv_MeDIP, MEDIPS.createSet(file="M:/shared/Hydroxymethylation Project/Raw BAM Files - MACS Peaks/22Rv1/mC/22RV1_rp3CoordSrt.bam", BSgenome=BSgenome, extend=extend, shift=shift, window_size=ws, chr.select=chr.select, paired=TRUE))
#RWPE-1 Input
RW_Input = MEDIPS.createSet(file="M:/shared/Hydroxymethylation Project/Raw BAM Files - MACS Peaks/RWPE-1/mC/RWPE1_InputCoordSrt.bam", BSgenome=BSgenome, extend=extend, shift=shift, window_size=ws, chr.select=chr.select, paired=TRUE)
#22Rv1 Input
Rv_Input = MEDIPS.createSet(file="M:/shared/Hydroxymethylation Project/Raw BAM Files - MACS Peaks/22Rv1/mC/22RV1_InputCoordSrt.bam", BSgenome=BSgenome, extend=extend, shift=shift, window_size=ws, chr.select=chr.select, paired=TRUE)

#CpG Density analysis
CS=MEDIPS.couplingVector(pattern="CG", refObj=RW_MeDIP[[1]])

#MeDIPS CNV Analysis
#RW v CR1
mr.edgeR = MEDIPS.meth(MSet1=RW_MeDIP, MSet2=Rv_MeDIP, CSet=CS, ISet1 = RW_Input, ISet2=Rv_Input, p.adj="bonferroni", diff.method="edgeR", MeDIP=T, CNV=T, minRowSum=10)

#Annotate the master file (ENSID) to hg19
anno.mart.gene=MEDIPS.getAnnotation(host="grch37.ensembl.org", dataset=c("hsapiens_gene_ensembl"),annotation=c("GENE"), chr=chr.select)
mr.edgeR2 = MEDIPS.setAnnotation(regions=mr.edgeR, annotation=anno.mart.gene)
#Update to HGNC gene symbols
genes <- mr.edgeR2$'1_id'
GeneList <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id", "hgnc_symbol"),values=genes,mart=useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl",host="grch37.ensembl.org"))
colnames(GeneList)[1] <- c("1_id")
CNVAnnot <- merge(mr.edgeR2, GeneList, by="1_id", all=TRUE)
mr.edgeR.s = MEDIPS.selectSig(results=CNVAnnot, p.value=0.1, adj=T, ratio=1, bg.counts=NULL,CNV=T)

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
