### Common functions for formatting Software inputs 


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
    
    #remove extensions
    injections <- as.vector(sapply(injections, file_path_sans_ext))
    all(sapply(injections, is.element, exp ))
}

guessInjection <- function(varname, exp){
    injections <- experiments[[exp]]
    as.character(names(which(sapply(injections,  grepl, varname))))
}

guessExperiment_wide <- function(exp, varnames){
    # exp: one of the experiments: i.e. experiments[[1]]
    # varnames: column names of the data frame
    any(sapply(exp, grepl, varnames))
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

guessSep <- function(filename){
    extension <- file_ext(filename)
    if(extension == "tsv") return ("\t")
    else if(extension == "csv") return (",")
    else return(NA)
}

