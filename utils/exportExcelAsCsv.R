library(dplyr, warn.conflicts = FALSE)
library(readxl)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript exportExcelAsCsv.R outputFile.csv")
}

inputFile <- args[1]
outputFile <- args[2]

loadExcel <- function(filepath) {
  sheets <- excel_sheets(filepath)
  datasets <- sapply(
    sheets,
    function(s) read_excel(filepath, s) %>% mutate_all(as.character),
    simplify = FALSE
  )
  bind_rows(datasets, .id = "Sheet")
}

dataset <- loadExcel(inputFile)
write.csv(dataset, file = outputFile, na = "", row.names = FALSE)