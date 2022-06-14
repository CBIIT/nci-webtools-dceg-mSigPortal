# convert Database/ .RData files to gzip JSON
library(jsonlite)
library(R.utils)

paths <- c("Signature", "Exposure")
files <- list.files(path = paths, pattern = "*.RData")

lapply(files, function(e) {
    rData <- get(load(x))
    jsonFile <- sub(".RData", ".json", e)
    write_json(rData, auto_unbox = TRUE, path = jsonFile)
    gzip(jsonFile, destname = sprintf("%s.gz", jsonFile), overwrite = TRUE, remove = TRUE)
})