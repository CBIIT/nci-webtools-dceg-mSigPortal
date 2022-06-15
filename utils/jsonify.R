# convert Database/ .RData files to gzip JSON
library(jsonlite)
library(R.utils)


# Signature (Seqmatrix) Public Data
paths <- c("Seqmatrix")
files <- list.files(path = paths, pattern = "*refdata.RData", full.names = TRUE)
lapply(files, function(e) {
    print(paste0("Processing ", e))
    rData <- get(load(e))
    jsonFile <- sub(".RData", ".json", e)
    write_json(rData, path = jsonFile)
    gzip(jsonFile, destname = sprintf("%s.gz", jsonFile), overwrite = TRUE, remove = TRUE)
    print(paste0("Done"))
})