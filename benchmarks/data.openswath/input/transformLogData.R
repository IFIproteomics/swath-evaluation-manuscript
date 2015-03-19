library(stringr)

transformlogdata <- function(filename){
    # data in datalog is actually log2(intensity) --> we need to transform it to 2^(log2(intensity))
    datalog <- read.csv(filename,header=T)
    datalog[datalog == 0] <- NA
    data <- 2 ** datalog[,2:7]
    data <- as.data.frame(cbind(as.character(datalog[,1]), data))
    names(data) <- c("protein", "A1", "A2", "A3", "B1", "B2", "B3")
    
    outputfile <- paste0(str_sub(filename, 0, -5),"_noLog.csv")
    
    write.csv(data, outputfile, row.names=F )
}



workdir <- "/Users/napedro/Dropbox/PAPER_SWATHbenchmark/hye.r/data.openswath/input"

setwd(workdir)

filelist <- list.files(pattern = "Sample_Quant_exp\\.csv$")


xx <- lapply(filelist, transformlogdata)
