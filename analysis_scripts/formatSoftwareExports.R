rm(list=ls())

loadLibrary <- function(lib) { 
    if(!require(lib, character.only = T)) { 
        install.packages(lib)
        library(lib, character.only = T)
    } 
}

loadLibrary("data.table")
loadLibrary("reshape2")
loadLibrary("dplyr")
loadLibrary("tidyr")
loadLibrary("tools")



#working_dir <- "/Users/napedro/Dropbox/PAPER_SWATHbenchmark/hye.r/data.peakview/RAW.PeakView.output"
#working_dir <- "/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/Skyline"
#working_dir <- "/Users/napedro/Dropbox/PAPER_SWATHbenchmark/hye.r/data.new.openswath/Raw_OpenSWATH_Output"
#working_dir <-"/Users/napedro/Dropbox/PAPER_SWATHbenchmark/hye.r/data.new.openswath/Raw_OpenSWATH_Output"
#working_dir <-"/Users/napedro/Dropbox/PAPER_SWATHbenchmark_prv/Skyline/transitions/Qvalue0.02_top_ranked_transitions"
working_dir <-"/Users/napedro/Dropbox/tmp_wrk_home2/SWATHbenchmark/DIA-umpire/round1/peptides"


software_source <- "DIAumpire"    # Options: "Spectronaut", "PeakView", "Skyline", "openSWATH", "DIAumpire"
input_format <- "wide"          # Options: "long", "wide"

results_dir <- "formatted_files"

source("fswe.variables.R")
#q_filter_threshold <- 0.005
source("fswe.functions.R")
source("fswe.datasets.R")

if(!file.exists(file.path(working_dir, results_dir))) { dir.create(file.path(working_dir, results_dir)) }
AllInputFiles = list.files( path=working_dir, pattern=input.extension, full.names= FALSE )


## Ops related to the quantitation variable ##
sumquant <- paste0("sum(", quantitative.var, ")")
medianquant <- paste0("median(", quantitative.var,")")
####

generateReports <- function(experiment_file){

    # Read file
    #  experiment_file <- AllInputFiles[1]
    cat(paste0("Generating peptide report for ", experiment_file, "\n"))
    experiment_file <- file.path(working_dir, experiment_file)
    df <- read.table(experiment_file, na.strings= nastrings,
                             header=T, sep=guessSep(experiment_file),
                             , stringsAsFactors =F, fill = T)
    
    # Attach specie and remove peptides belonging to multiple species, and not considered species ("NA"s)
    if(!is.na(q_filter_threshold)){
        df <- eval( substitute(filter(df, var < q_filter_threshold), list( var = as.name(qvalue.var)) ) )
    }
    df <- df %>% rowwise()
    df <- eval( substitute(mutate(df, "specie" = guessSpecie(var)), list(var = as.name(protein.var)) ) ) 
    df <- filter(df, specie != "NA", specie != "multiple")
    
    experiment <- NA
    if(input_format == "wide"){
        # At this point, it is better to transform the data frame into a long-format data frame, and
        # continue the analysis commonly for wide and long inputs.
        
        # find the columns containing the quantitative values (pivoted values). 
        # They will be use to gather the key-value pairs
        experiment <- which(sapply(experiments, guessExperiment_wide, colnames(df) ))
        
        df <- df %>% gather_(filename.var, quantitative.var, 4:9) %>%  # TODO: change this hard-coded 6:11
                    arrange_(protein.var, sequence.mod.var)
        
    }else if(input_format == "long"){
        # Guess the experiment by the filename.var column
        injections <- distinct(select_(df, filename.var))  
        experiment <- which(sapply(experiments, guessExperiment, injections))
    }

    # Remove any duplicate: (based on: injectionfile, ModifiedSequence, ChargeState)
    # (These cases are likely to be due to several entries in the library)
    data <- df %>% distinct_(filename.var, sequence.mod.var, charge.var)  # I am not sure I removed all duplicates, or there is still one value per duplicate
    
    #Remove NAs of the quant variable
    data <- data %>% na.omit() 
    
    # For each peptide: Sum quantitative values of charge states
    data <- data %>% group_by_(filename.var, sequence.mod.var, protein.var , "specie") %>% 
        summarise_(  quant_value = sumquant ) # proteinID = protein.var, specie = "specie" , 
    
    peptides_wide <- spread_(data, filename.var, "quant_value") 
    
    # Change variable names of injections by their right name (removing additions from software to the name, etc)
    common_names <- names(peptides_wide)[c(1:3)]
    inj_names <- names(peptides_wide)[-c(1:3)]
    inj_names <- as.character(sapply(inj_names, guessInjection, experiment))
    names(peptides_wide) <- c( common_names, inj_names )
    # names(peptides_wide) <- gsub("\\.[^\\.]*$","", names(peptides_wide))
    
    experiment.order <- match(experiments[[experiment]], names(peptides_wide[-c(1:3)])) + 3
    peptides_wide <- peptides_wide[, c(c(1:3), experiment.order)]

    
    names(peptides_wide) <- c("sequenceID", "proteinID", "specie", "A1", "A2", "A3", "B1", "B2", "B3") #Rename the samples is unnecessary, but...

    expfile_noext <- file_path_sans_ext(basename(experiment_file))
    
    write.table(peptides_wide, file=file.path(working_dir, results_dir ,paste0(expfile_noext, "_peptides.tsv")), 
                sep="\t", row.names=F, col.names=T)


    
    ## PROTEIN REPORT
    cat(paste0("Generating protein report for ", experiment_file,"\n"))
    
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


