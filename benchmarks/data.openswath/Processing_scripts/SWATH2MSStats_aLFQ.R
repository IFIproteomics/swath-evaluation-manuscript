### Install aLFQ (only once)
source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
biocLite("graph")
biocLite("RBGL")
install.packages("devtools")
#install Rtools 3.1 from http://cran.r-project.org/bin/windows/Rtools/
find_rtools()
library(devtools)
install_github("grosenberger/aLFQ")
library(reshape)

### 1.1 Loading the data
rm(list=ls())
library(SWATH2MSstats)

# set working directory
setwd("~/IMSB/Experiments/QTOF_Panning/1503_Pedro_comparisons/")
# Input data file (openSWATH output with or without requant)
file.name <- "TTOF5600_1ug_32w_feature_alignment.tsv"

# File to assign file names to condition, biol. replicates and run number
#will be created later on.
annotation.file <- "TTOF5600_1ug_32w_sample_annotation.txt"

# load data
data <- fread(file.name, sep="\t", header=TRUE)
data <- as.data.frame(data)

### 1.2 Annotate the data
# list number and different Files present
nlevels(factor(data$align_origfilename)) # prints the number of runs
levels(factor(data$align_origfilename)) # prints the names of the runs

# create a tab-delimited file (e.g. in Excel) with the columns
#Filename, Condition, BioReplicate, Run
# for the filename column it is enough to provide: userid_K140613_001_SW_xxx
sample_annotation <- read.delim2(file.path(getwd(), annotation.file), dec=".", sep="\t", header=TRUE)

# annotate data
data <- sample_annotation(data, sample_annotation)
# OPTIONAL: for human, shorten Protein Name to remove non-unique information
#(sp|Q9GZL7|WDR12_HUMAN --> Q9GZL7)
head(unique(data$ProteinName))

### 1.3 Filter the data
data.filtered.decoy <- assess_decoy_rate(data) # count and remove decoys
#data.filtered.mscore <- filter_mscore(data.filtered.decoy, 2)

data <- convert4aLFQ(data.filtered.decoy)


### actual aLFQ

library(aLFQ)

prots <- ProteinInference(data, peptide_method = "top", peptide_topx = 3, peptide_strictness = "loose", peptide_summary = "mean", transition_topx = 3, transition_strictness = "loose", transition_summary = "sum", fasta = NA, apex_model = NA, combine_precursors = FALSE, combine_peptide_sequences = FALSE, consensus_proteins = FALSE, consensus_peptides = TRUE, consensus_transitions = TRUE)

#write.csv(prots, file="5600_64w_alfq_proteininference_output.csv", row.names=FALSE, quote=FALSE)

prots_c <- cast(prots,protein_id~run_id)

write.csv(prots_c, file="TTOF5600_1ug_32w_alfq_proteininference1_output.csv", row.names=FALSE, quote=FALSE)

