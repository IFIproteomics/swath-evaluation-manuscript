rm(list=ls())

if(!require(data.table)) { install.packages("data.table") ;library(data.table) }
if(!require(reshape2)) { install.packages("reshape2") ;library(reshape2) }
if(!require(dplyr)) { install.packages("dplyr") ;library(dplyr)}
if(!require(tidyr)) { install.packages("tidyr") ;library(tidyr)}
if(!require(tools)) { install.packages("tools") ;library(tools)}

working_dir <- "/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/Spectronaut7"
software_source <- "Spectronaut"
input_format <- "long"  # Options: "long", "wide"
results_dir <- "formatted_files"
setwd(working_dir)

if(!file.exists(results_dir)) { dir.create(file.path(working_dir, results_dir)) }

AllInputFiles = list.files( path=working_dir, pattern="*.tsv", full.names= FALSE )
#
species <- vector(mode="list", length=3)
names(species) <- c("HUMAN", "YEAST", "ECOLI")
species[[1]] <- "_HUMAN"
species[[2]] <- "_YEAS"
species[[3]] <- "_ECOLI"

experiments <- vector(mode="list", length=4)
names(experiments) <- c("5600-32w", "5600-64w", "6600-32w", "6600-64w")
# for each experiment:
#   three first ones are sample A, other three ones are sample B
experiments[[1]] <- c("lgillet_L150206_001", "lgillet_L150206_003", "lgillet_L150206_005",   # A
                      "lgillet_L150206_002", "lgillet_L150206_013", "lgillet_L150206_014")   # B 

experiments[[2]] <- c("lgillet_L150206_007", "lgillet_L150206_009", "lgillet_L150206_011",   # A
                      "lgillet_L150206_008", "lgillet_L150206_010", "lgillet_L150206_012")   # B

experiments[[3]] <- c("lgillet_I150211_002", "lgillet_I150211_004", "lgillet_I150211_006",   # A
                      "lgillet_I150211_003", "lgillet_I150211_005", "lgillet_I150211_007")   # B

experiments[[4]] <- c("lgillet_I150211_008", "lgillet_I150211_010", "lgillet_I150211_012",   # A
                      "lgillet_I150211_009", "lgillet_I150211_011", "lgillet_I150211_013")   # B



## Select Software-depending variable names #########
if(software_source == "Spectronaut"){
    quantitative.var <- "FG.NormalizedTotalPeakArea"
    protein.var <- "EG.ProteinId"
    filename.var <- "R.FileName"
    sequence.mod.var <- "EG.ModifiedSequence"
    charge.var <- "FG.Charge"    
}

#####################################################

## Ops related to the quantitation variable ##
sumquant <- paste0("sum(", quantitative.var, ")")
medianquant <- paste0("median(", quantitative.var,")")
####



substrRight <- function(x, n) { substr(x, nchar(x)-n+1, nchar(x)) }
rmlastchars <- function(x, n) { substr(x, 1, nchar(x) - n) }

take1stentry <- function(entries){
  first_entry <- strsplit(as.character(entries), "\\||\\/", fixed=F, perl=T)
  first_entry <- first_entry[[1]][4]
  if(substrRight(first_entry, 3)=="/sp") first_entry <- rmlastchars(first_entry, 3)
  return(first_entry)
}

guessExperiment <- function(exp, injections){
    # exp: one of the experiments: i.e. experiments[[1]]
    # injections: list of injections (preferably a unique list!) of the experiment
    all(sapply(injections, is.element, exp ))
}

guessSpecie <- function(proteinid){
    sp <- names(which(sapply(species, grepl, proteinid)))
    if(length(sp) == 0) sp <- "NA"
    if(length(sp) > 1) sp <- "multiple"
    sp
}

sum_top_n <- function(values, n, minimum = 1){
    # This top N approach is INDIVIDUAL, that is, there is no consensus among replicates to 
    # choose the top N peptides.
    if(length(which(!is.na(values))) < minimum) {return (NA)}
    if (n > length(values)) n <- length(values)
    sum(sort(values, decreasing=T)[1:n], na.rm=T)
}

generateReports <- function(experiment_file){

    # Read file
    #experiment_file <- "20150317_182627_TTOF5600_1ug_64w_Spectronaut7_Report.tsv"
    print(paste0("Generating peptide report for ", experiment_file))
    df <- read.table(experiment_file, na.strings="NaN",
                             header=T, sep="\t",
                             comment.char="", quote = "", stringsAsFactors =F)
    
    # Attach specie and remove peptides belonging to multiple species, and not considered species ("NA"s)
    df <- df %>% rowwise()
    df <- eval( substitute(mutate(df, "specie" = guessSpecie(var)), list(var = as.name(protein.var)) ) ) 
    df <- filter(df, specie != "NA", specie != "multiple")
    
    # Guess the experiment by the filename.var column
    injections <- distinct(select_(df, filename.var))  
    experiment <- which(sapply(experiments, guessExperiment, injections))
    
    # Remove any duplicate: (based on: injectionfile, ModifiedSequence, ChargeState)
    # (These cases are likely to be due to several entries in the library)
    data <- df %>% distinct_(filename.var, sequence.mod.var, charge.var)  # I am not sure I removed all duplicates, or there is still one value per duplicate
    
    # For each peptide: Sum quantitative values of charge states
    data <- data %>% group_by_(filename.var, sequence.mod.var) %>% 
        summarise_( proteinID = protein.var, specie = "specie" ,  quant_value = sumquant )  

    # Case of long formats: pivot filename values (and assign already samples?)
    if(input_format == "long"){
        peptides_wide <- spread_(data, filename.var, "quant_value") 
        experiment.order <- match(experiments[[experiment]], names(peptides_wide[-c(1:3)])) + 3
        peptides_wide <- peptides_wide[, c(c(1:3), experiment.order)]
    }
    else if(input_format == "wide"){
        #Not yet implemented!
    }
    
    names(peptides_wide) <- c("sequenceID", "proteinID", "specie", "A1", "A2", "A3", "B1", "B2", "B3") #Rename the samples is unnecessary, but...

    expfile_noext <- file_path_sans_ext(basename(experiment_file))
    
    write.table(peptides_wide, file=file.path(working_dir, results_dir ,paste0(expfile_noext, "_peptides.tsv")), 
                sep="\t", row.names=F, col.names=T)


    
    ## PROTEIN REPORT
    print(paste0("Generating protein report for ", experiment_file))
    
    proteins_wide <- peptides_wide %>% 
                        select(-sequenceID) %>% 
                        arrange(proteinID, specie) %>%
                        group_by(proteinID, specie) %>%  
                        summarise_each(funs(sum_top_n(., 3, 2)))  # , n_distinct(.) (if you want to report the num of peptides)

    # Remove "empty" proteins (all values are NAs). I wish I could find a more elegant way to do it. I am tired.
    proteins_wide <- filter(proteins_wide, !is.na(A1) | !is.na(A2) | !is.na(A3) | !is.na(B1) | !is.na(B2) | !is.na(B3))
    
    write.table(proteins_wide, file=file.path(working_dir, results_dir ,paste0(expfile_noext, "_proteins.tsv")), 
                sep="\t", row.names=F, col.names=T)
    
    
}


nix <- sapply(AllInputFiles, generateReports)


