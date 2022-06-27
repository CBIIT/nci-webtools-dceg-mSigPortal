library(dplyr, warn.conflicts = FALSE)
library(vroom)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript combineDatasets.R 'inputFileGlobPattern' processorFunction outputFile.csv")
}

inputFiles <- Sys.glob(args[1])
processorFunction <- args[2]
outputFile <- args[3]

extract <- function(str, separator, index) {
    unlist(strsplit(str, separator))[index]
}

combineAssociationFiles <- function(x) {
    x %>% 
        mutate(
            Study=extract(basename(filepath), '_' , 1),
            Dataset=extract(basename(filepath), '_' , 2),
            .before=Cancer_Type
        )
}

datasets <- sapply(inputFiles, function(f) get(load(f)), simplify=F)
combinedDatasets <- bind_rows(datasets, .id="filepath")
processedDatasets <- do.call(processorFunction, list(x=combinedDatasets)) %>% select(-filepath)
vroom_write(processedDatasets, file = outputFile, delim=",", na="")

