library(dplyr, warn.conflicts = FALSE)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript exportJsonAsCsv.R inputFile.json columnName=columnValue outputFile.csv")
}

inputFile <- args[1]
additionalColumn <- unlist(strsplit(args[2], '='))
outputFile <- args[3]

loadJson <- function(filepath, additionalColumn) {
  dataset <- read_json(filepath, simplifyVector = T) %>% 
    rowwise() %>% 
    mutate(across(where(~ class(.x) == "list"), ~ paste(., collapse = " ")))
  columnName = additionalColumn[1]
  columnValue = additionalColumn[2]
  dataset[[columnName]] <- columnValue
  dataset
}

dataset <- loadJson(inputFile, additionalColumn)
write.csv(dataset, file = outputFile, na="", row.names = F)

