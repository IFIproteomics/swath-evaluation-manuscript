## Select Software-depending variable names #########
quantitative.var <- NA
quantitative.var.tag <- NA
protein.var <- NA
filename.var <- NA
sequence.mod.var <- NA
charge.var <- NA
q_filter_threshold <- NA
decoy.var <- NA
input.extension <- "*.csv$"
nastrings = "NA"

if(software_source == "Spectronaut"){
    quantitative.var <- "FG.NormalizedTotalPeakArea"
    protein.var <- "EG.ProteinId"
    filename.var <- "R.FileName"
    sequence.mod.var <- "EG.ModifiedSequence"
    charge.var <- "FG.Charge" 
    nastrings <- "NaN"
    input.extension <- "*.tsv$"
}

if(software_source == "test"){
    quantitative.var <- "sumArea"
    protein.var <- "ProteinName"
    filename.var <- "ReplicateName"
    sequence.mod.var <- "PeptideSequence"
    charge.var <- "PrecursorCharge" 
    nastrings <- "NA"
    input.extension <- "*.tsv$"
}

if(software_source == "DIAumpire"){
    quantitative.var <- "area"
    protein.var <- "Protein"
    filename.var <- "ReplicateName"
    sequence.mod.var <- "ModSeq"
    charge.var <- "Charge" 
    nastrings <- "NA"
    input.extension <- "*.csv$"
}

if(software_source == "Skyline"){
    q_filter_threshold <- 0.15
    qvalue.var <- "annotation_QValue"
    quantitative.var <- "TotalArea"
    protein.var <- "ProteinName"
    filename.var <- "FileName"
    sequence.mod.var <- "ModifiedSequence"
    charge.var <- "PrecursorCharge"  
    decoy.var <- "IsDecoy"
    nastrings <- "#N/A"
}
if(software_source == "PeakView"){
    quantitative.var <- "TotalAreaFragment"
    protein.var <- "Protein"
    filename.var <- "R.FileName"
    sequence.mod.var <- "Peptide"
    charge.var <- "Precursor Charge"
    input.extension <- "*.tsv$"
}
if(software_source == "openSWATH"){
    #q_filter_threshold <- 0.05
    #qvalue.var <- "score"
    quantitative.var.tag <- "Intensity_"
    quantitative.var <- "Intensity"
    protein.var <- "Protein"
    filename.var <- "FileName"
    sequence.mod.var <- "Peptide"
    charge.var <- "Charge"  
    decoy.var <- "IsDecoy"
    input.extension <- "*.tsv$"
}

quantitative.var <- gsub(" ", ".", quantitative.var)
protein.var <- gsub(" ", ".", protein.var)
filename.var <- gsub(" ", ".", filename.var)
sequence.mod.var <- gsub(" ", ".", sequence.mod.var)
charge.var <- gsub(" ", ".", charge.var)
if(!is.na(decoy.var)){
    decoy.var <- gsub(" ", ".", decoy.var)
}

#####################################################

