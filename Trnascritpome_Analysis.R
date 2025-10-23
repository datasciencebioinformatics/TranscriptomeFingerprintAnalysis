#!/usr/bin/env Rscript
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
dataInfo<-data.frame(read.table(file_dataInfo,header = TRUE,row.names=4))

# Sort the rows of the table
dataInfo<-dataInfo[ order(row.names(dataInfo)), ]

# Transform Run coluln into factor
dataInfo[,batch_run]<-as.factor(dataInfo[,batch_run])

# Create DESeq object from data, grouped by Sample_Group
DESeqData_group_phenotype <- DESeqDataSetFromMatrix(countData = tableReadsCount,colData = dataInfo,design = as.formula(paste("+",group_phenotype)))

############################################################################################################################################
# Pre-filtering - removing rows in which there are no reads or nearly no reads
# Filter DESeq object grouped by Sample_Group
DESeqData_group_phenotype <- DESeqData_group_phenotype[ rowSums(counts(DESeqData_group_phenotype)) > 1, ]

############################################################################################################################################

# Differential expression analysis
# Comparison will be the last level of this variable over the first level - Sample_Group
DESeqData_group_phenotype <- DESeq(DESeqData_group_phenotype)
############################################################################################################################################
# Get TPM count
