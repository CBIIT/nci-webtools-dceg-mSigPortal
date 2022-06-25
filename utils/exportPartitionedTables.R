library(dplyr, warn.conflicts = FALSE)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript exportPartitionedTables.R inputFile inputDataKey group1,group2,group3 outputFolder/")
}

inputFile <- args[1]
inputDataKey <- args[2]
groups <- unlist(strsplit(args[3], ','))
outputFolder <- args[4]

load(inputFile)

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

writeJsonGroups(
  get(inputDataKey), 
  groups
)
