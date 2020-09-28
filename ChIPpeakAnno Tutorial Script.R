# ChIPpeakAnno Tutorial v1.0
# Annotation of genomic loci (same method used by TCAG)
#
# Purpose: Annotate a list of genomic regions with a
#          given chromosome, start, and end region.
#          Annotations will include gene names
#          (Ensembl gene and HGNC IDs) and will
#          specify whether the location is intragenic,
#          intergenic, or promoter (and the distance
#          from each feature).
#
# Version: v1.0
# Date: 2017-05-26
# Author: Shivani Kamdar
#
# Input: Tab-delimited .txt files - can have
#        any data in it, as long as it includes 
#        three key columns:
#        one with the chromosome ID (format: chr1,
#        chr2, etc.), one with the start co-ordinate,
#        and one with the end co-ordinate.
#        FOR EASE OF USE: name these columns "chr",
#        "start", and "end". If you don't, you'll
#        need to change column names manually (see
#        below)
# Output: Same file, but with gene IDs and genomic
#         locations added.
#
# Dependencies: ChIPpeakAnno, biomaRt,
#               EnsDb.Hsapiens.v75
#               (v75 is for hg19-aligned data). If
#               you have data aligned to a different
#               genome build, change the version to
#               the appropriate one in all the code
#               before running.
#
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

#Importing your data
## Instructions: replace "mydata.txt" with the name
## of your file. If it is in a different folder,
## indicate the folder path in the filename (ie:
## C:/my/file/path/filename.txt)

Basetable <- read.table("mydata.txt", header=TRUE, sep="\t")
## If you need to change your column names in R, see
## below for an example.
## colnames(Basetable) <- c("chr","start","end")

#################################################################################################################################################################

#CpG island annotation (performed first so that we can adjust column names accordingly)
#ExtractCpGs
library(AnnotationHub)
hub <- AnnotationHub()
cpgs <- hub[["AH5086"]]

#Format the data to GRanges
CpGChip <- makeGRangesFromDataFrame(Basetable,keep.extra.columns=TRUE,ignore.strand = TRUE, seqinfo=NULL, seqnames.field="chr",start.field="start",end.field="end")
annotatedCpGPeak <- annotatePeakInBatch(CpGChip,AnnotationData=cpgs)
AnnotCpGChip <- as.data.frame(annotatedCpGPeak)

#Format the data to prepare for non-CpG annotation
colnames(AnnotCpGChip)[c(38:length(colnames(AnnotCpGChip)))] <- paste("CpG_",colnames(AnnotCpGChip)[c(38:length(colnames(AnnotCpGChip)))],sep="")
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


#Formatting the output table and merging
##Uncomment if you are having issues with formatting
##AnnotChip$EnsGene <- NA
##AnnotChip$EnsGene <- rownames(AnnotChip)
##AnnotChip$EnsGene <- substring(AnnotChip$EnsGene,8,22)



## Annotate with HGNC IDs
mart <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
HGNCAnnot <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),filters="ensembl_gene_id",values=AnnotChip$feature,mart=mart)
colnames(HGNCAnnot) <- c("feature", "HGNC_ID")
## Merge
FinalAnnot <- merge(AnnotChip, HGNCAnnot,by="feature", all=TRUE)

#Export the file - and we're done!
write.table(FinalAnnot, "Annotated Peak Data.txt", sep="\t")
