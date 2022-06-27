library(dplyr, warn.conflicts = FALSE)
library(vroom)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript exportPartitionedTables.R inputFile group1,group2,group3 outputFolder")
}

inputFile <- args[1]
groups <- unlist(strsplit(args[2], ','))
outputFolder <- args[3]

inputData <- vroom(inputFile)

dir.create(outputFolder, recursive = T, showWarnings = F)
setwd(outputFolder)

writeJsonGroup <- function(x, group) {
  folderPath <- paste0(group, collapse = '/')
  filePath <- paste(folderPath, 'data.json', sep='/')
  dir.create(folderPath, recursive = T, showWarnings = F)
  write_json(x, filePath)
}

writeJsonGroups <- function(datatable, groups) {
  datatable %>% 
    group_by(across(all_of(groups))) %>%
    group_walk(writeJsonGroup, .keep=T)
}

writeJsonGroups(inputData, groups)
