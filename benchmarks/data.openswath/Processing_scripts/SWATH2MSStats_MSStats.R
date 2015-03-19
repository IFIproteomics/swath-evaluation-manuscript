### Install MSStats (only once)
# step 1: install dependency packages
# install.packages(c("gplots","lme4","ggplot2","reshape","data.table","Rcpp"))
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("limma","marray","preprocessCore","MSnbase"))
# biocLite("Biobase")

# step 2: select 'Install package(s) from local zip files...' under 'Packages' in Menu bar. Then select 'MSstats.daily_2.1.6.zip' from your local directory

# step 3: load the library of MSstats.daily
# library(MSstats.daily)

# step 4: getting started


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

# select necessary columns
data <- filter_openSWATH_output(data)
# remove the iRT peptides (usually not in there anyway)
data <- data[grep("iRT", data$ProteinName, invert=TRUE),]

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
#data.filtered.mscore <- filter_mscore_requant(data.filtered.decoy, 0.01, 0.8)
#data.filtered.mscore <- filter_mscore_condition(data.filtered.decoy, 0.01, 4)

#data <- filter_proteotypic_peptides(data.filtered.mscore)
#data <- filter_all_peptides(data.filtered.mscore)
#data.filtered.max <- filter_on_max_peptides(data.filtered.mscore, 5)
data <- data.filtered.decoy


### 2.1 Results on peptide level
#data.peptide <- data
#data.peptide$aggr_Fragment_Annotation <- NULL
#data.peptide$aggr_Peak_Area <- NULL
#write.csv(data.peptide, file="TTOF5600_1ug_32w_peptide_level_output.csv", row.names=FALSE, quote=FALSE)
#write.csv(data.peptide, file="5600_64w_peptide_level_output.csv", row.names=FALSE, quote=FALSE)

### 2.2.1 MSstats within R
#nrow(data)
#raw <- disaggregate(data)
#MSstats.input <- convert4MSstats(raw, bash.script = FALSE)
#write.csv(MSstats.input, file="5600_32w_MSstats_input_trs.csv", row.names=FALSE, quote=FALSE)


### 2.2.2 MSstats using the bash script
data <- convert4MSstats(data, bash.script = TRUE)

write.csv(data, file='TTOF5600_1ug_32w_feature_alignment_mod.csv', row.names=FALSE, quote=FALSE)
head(data)

## in brutus: sh ./featurealigner2msstats_gene.bash TTOF5600_1ug_32w_feature_alignment_mod.csv TTOF5600_1ug_32w_feature_alignment_mod_trs.csv

## load data back to R
raw <- fread("TTOF5600_1ug_32w_feature_alignment_mod_trs.csv", sep=",", header=TRUE)
raw <- as.data.frame(raw)
#head(raw)
dim(raw)
unique(raw$ProteinName)
### actual MSStats analysis
library(MSstats.daily)
#head(raw)
sum(raw$Intensity < 0, na.rm=T)
sum(raw$Intensity == 0, na.rm=T)
raw[raw$Intensity %in% 0,"Intensity"]<-NA
dim(raw)
unique(raw$ProteinName)
#msstats.data <- dataProcess(raw, logTrans=2, normalization=TRUE,betweenRunInterferenceScore=FALSE)
# Data preprocessing and QC in MSstats
data1 <- dataProcess(raw, logTrans=2, normalization="quantile",nameStandards=NULL, betweenRunInterferenceScore=FALSE,fillIncompleteRows=TRUE)
dim(data1)
unique(data1$PROTEIN)
#data<-dataProcess(raw, logTrans=2, normalization="quantile",nameStandards=NULL,
#            betweenRunInterferenceScore=FALSE,fillIncompleteRows=FALSE,
#            FeatureSelection=TRUE,lambda=seq(0,1.4,0.1),eta=0.05, address="")

# Profile plot of peptides for every protein
# dataProcessPlots(data=data1.trans,type="ProfilePlot",featureName="Peptide",ylimUp=FALSE,ylimDown=FALSE,scale=FALSE,x.axis.size=10,y.axis.size=10,text.size=4,text.angle=0,legend.size=7,dot.size=3,width=10, height=10, which.Protein="all")

# QC plot for every protein
# dataProcessPlots(data=data1,type="QCPlot", ylimUp=FALSE,ylimDown=FALSE,scale=FALSE, x.axis.size=10,y.axis.size=10,text.size=4,text.angle=0,legend.size=7,dot.size=3, width=10, height=10, which.Protein="all")

#Condition plot of peptides for every protein
# dataProcessPlots(data=data1,type="ConditionPlot",ylimUp=FALSE,ylimDown=FALSE,scale=TRUE,interval="SD",x.axis.size=10,y.axis.size=10,text.size=4,text.angle=0,legend.size=7,dot.size=3,width=10, height=10, which.Protein="all")

#To show the order of the conditions
#levels(data1$GROUP_ORIGINAL)

#Comparison Matrix: creating row per comparison
#comparison1 <- matrix(c(-1,1), nrow=1)

#bind all the rows together into a matrix
#comparison <- rbind(comparison1)
#row.names(comparison) <- c("A.vs.B")

#Testing/ group comparison
#data2 <- groupComparison(contrast.matrix=comparison, data=data1, labeled=FALSE,scopeOfBioReplication="expanded", scopeOfTechReplication="expanded",interference=TRUE,missing.action="nointeraction")
#write.csv(data2$ComparisonResult, file="TTOF5600_1ug_32w_groupcomparison.csv",row.names=FALSE)

#To save the results in wide formats
#data2_wide_FC <- dcast(data2$ComparisonResult, Protein ~ Label, value.var="log2FC")
#data2_wide_pval <- dcast(data2$ComparisonResult, Protein ~ Label, value.var="adj.pvalue")
#data2_wide <- merge(data2_wide_FC, data2_wide_pval, by.x = "Protein", by.y="Protein")

#write.csv(data2_wide, file ="TTOF5600_1ug_32w_groupComparison_wide.csv", row.names=FALSE)

#Plots
#?groupComparisonPlots
#Volcano Plot
#groupComparisonPlots(data=data2$ComparisonResult,type="VolcanoPlot",sig=0.05,FCcutoff=1.5,ylimUp=40,ylimDown=FALSE,xlimUp=FALSE,x.axis.size=10,y.axis.size=10,dot.size=3,text.size=3.4,legend.size=7,ProteinName=TRUE,ProteinNameLoc=1, clustering="both", width=10, height=10, which.Comparison="all")
#Heatmap
#groupComparisonPlots(data=data2$ComparisonResult,type="Heatmap",sig=0.05,FCcutoff=1.5,ylimUp=FALSE,ylimDown=FALSE,xlimUp=FALSE,x.axis.size=10,y.axis.size=10,text.size=3,legend.size=7,ProteinName=TRUE, numProtein=60, clustering="both", width=12, height=12, which.Comparison="all")
#Comparison Plot
#groupComparisonPlots(data=data2$ComparisonResult,type="ComparisonPlot",sig=0.05,ylimUp=FALSE,ylimDown=FALSE,xlimUp=FALSE,x.axis.size=10,y.axis.size=10,dot.size=3,text.size=4,legend.size=7,width=10, height=10, which.Comparison="all")
#?modelBasedQCPlots
#To test for normal distribution of errors (residuals against predicted values)
#modelBasedQCPlots(data2$ModelQC,type="ResidualPlots",axis.size=10,dot.size=3,text.size=7,legend.size=7,width=10, height=10,featureName=TRUE,which.Protein="all")

#To test for constant variance of errors (QQ plots)
#modelBasedQCPlots(data2$ModelQC,type="QQPlots", axis.size=10,dot.size=3,text.size=7,legend.size=7,width=10, height=10,feature.QQPlot="all",which.Protein="all")

### Quantification
#?quantification
#Model-based quantification for each condition/ for each biological samples per protein. 
#Quantification takes the processed data by dataProcess (data1) as input and automatically generate the
#quantification results (data.frame) wth long/ matrix format.

#groupquant <- quantification(data=data1,type="Group",format="matrix")

#groupquant <- quantification(data=data1,type="Group",format="matrix", scopeOfTechReplication="restricted", scopeOfBioReplication="restricted", interference=TRUE, missing.action="nointeraction",equalFeatureVar=TRUE)

#write.csv(groupquant,file="group_Quant_exp.csv",row.names=TRUE, quote= FALSE)

samplequant <- quantification(data1,type="Sample",format="matrix", scopeOfTechReplication="restricted", scopeOfBioReplication="restricted", interference=TRUE, missing.action="nointeraction", equalFeatureVar=TRUE)
#samplequant2 <- quantification(data1,type="Sample",format="matrix", scopeOfTechReplication="expanded", scopeOfBioReplication="expanded", interference=TRUE, missing.action="nointeraction", equalFeatureVar=TRUE)

write.csv(samplequant,file="TTOF5600_1ug_32w_Sample_Quant_exp_again.csv",row.names=TRUE)
#write.csv(samplequant2,file="TTOF5600_1ug_32w_Sample_Quant2_exp.csv",row.names=TRUE)

