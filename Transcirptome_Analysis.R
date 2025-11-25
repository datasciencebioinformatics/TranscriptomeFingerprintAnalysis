#!/usr/bin/env Rscript
## 1) Set the project folder either on Linux or in Window
### a) If in Linux
##### project_folder="/home/felipe/Documents/TranscriptomeFingerprintAnalysis/"

### b) Or in Windows
project_folder="C:/Users/User/Documents/GitHub/TranscriptomeFingerprintAnalysis/"

### c) set outputfolder
output_dir=project_folder
############################################################################################################################################
# Load library 
options(java.parameters = "-Xmx16000m")
library("RSQLite")
library("systemPipeR")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("gridExtra")
library("grid")
library("lattice")
library("org.Hs.eg.db")
library("xlsx")
library("topGO")
library("biomaRt")
library("matrixStats")
library("reshape")
#############################################################################################################################################
# Define functions
# Function to convert IDs
convertIDs <- function( ids, from, to, db) 
{
	# Return ID
	return(getBM(attributes = c(from, to), filters = "ensembl_transcript_id", values = ids, mart = db)[,2])   

}
############################################################################################################################################
# Read the inputdatafile argument
file_tableReadsCount=paste(project_folder,"GSE123658_read_counts.gene_level.tsv.txt",sep="")

# Read the inputdatafile2 argument
file_dataInfo=paste(project_folder,"Metadata.txt",sep="")
############################################################################################################################################																																																																																																			
# Set the analysis level
analysis_level="ENSEMBL"

# Initialize mart object
mart = useMart("ensembl")

############################################################################################################################################
# Setting the annotation
# If species is human

# Get the annotation
anno = org.Hs.eg.db
mapping = "org.Hs.eg.db"
dataset="hsapiens_gene_ensembl"
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset=dataset)

############################################################################################################################################
# Load maxtrix as table
tableReadsCount<-read.table(file_tableReadsCount,header = TRUE,row.names=1, check.names = FALSE)

# Sort the columns of the table
tableReadsCount<-tableReadsCount[ , order(names(tableReadsCount))]

# Load maxtrix as table
dataInfo<-data.frame(read.table(file_dataInfo,header = FALSE))

# Set colnames
colnames(dataInfo)<-c("sequence_file","diagnosis","id")

# Set colnames
rownames(dataInfo)<-dataInfo$id

# Set the order
dataInfo<-dataInfo[colnames(tableReadsCount),]

# Create DESeq object from data, grouped by Sample_Group
DESeqData_group_phenotype <- DESeqDataSetFromMatrix(countData = tableReadsCount,colData = dataInfo,design = as.formula(paste("~ diagnosis")))
############################################################################################################################################
# Differential expression analysis
# Comparison will be the last level of this variable over the first level - Sample_Group
DESeqData_group_phenotype <- DESeq(DESeqData_group_phenotype)
############################################################################################################################################
# Extract results for the default comparison
res <- results(DESeqData_group_phenotype,contrast=c("diagnosis","healthy", "T1D"))

# Print a summary of the results
summary(res)

# Get significant genes (e.g., with padj < 0.05)
significant_genes <- subset(res, padj < 0.05)

significant_genes$log2FoldChange
############################################################################################################################################
# Take the average expression in healthy +- sd, average expression in patients +- sd, number of interactions,
significant_genes<-significant_genes[which(significant_genes$log2FoldChange>1.5),]

# Add collumns for mean expression
significant_genes<-cbind(significant_genes,mean_healthy=0,sd_healthy=0,mean_T1D=0,sd_T1D=0)

# for each genes, calculate the averages
for (genes in rownames(significant_genes))
{
	#  Take the mean of the healthy samples
	significant_genes[genes,"mean_healthy"]<-mean(counts(DESeqData_group_phenotype,norm=TRUE)[genes,dataInfo[which(dataInfo$diagnosis == "healthy"),"id"]])
	significant_genes[genes,"sd_healthy"]<-sd(counts(DESeqData_group_phenotype,norm=TRUE)[genes,dataInfo[which(dataInfo$diagnosis == "healthy"),"id"]])
	significant_genes[genes,"mean_T1D"]<-mean(counts(DESeqData_group_phenotype,norm=TRUE)[genes,dataInfo[which(dataInfo$diagnosis == "T1D"),"id"]])
	significant_genes[genes,"sd_T1D"]<-sd(counts(DESeqData_group_phenotype,norm=TRUE)[genes,dataInfo[which(dataInfo$diagnosis == "T1D"),"id"]])	
}
############################################################################################################################################
