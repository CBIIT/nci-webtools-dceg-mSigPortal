# convert Database/ .RData files to gzip CSV

# Signature (Seqmatrix) Public Data
paths <- c("Database/Seqmatrix")
files <- list.files(path = paths, pattern = "*seqmatrix_refdata.RData", full.names = TRUE)
lapply(files, function(e) {
    print(paste0("Processing ", e))
    rData <- get(load(e))
    csvFile <- sub(".RData", ".csv.gz", e)
    write.csv(rData, file = gzfile(csvFile), row.names = FALSE)
    print(paste0("Done"))
})