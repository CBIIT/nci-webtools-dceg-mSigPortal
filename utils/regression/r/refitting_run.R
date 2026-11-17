#!/usr/bin/env Rscript
# Separate service preflight: real SATS/source loading and input validation.
# Scientific fitting needs the service's MAF, panel, clinical and reference files.
argv <- commandArgs(TRUE)
option <- function(key, default = NULL) {
  i <- match(key, argv)
  if (is.na(i)) return(default)
  if (i == length(argv)) stop("Missing value for ", key)
  argv[[i + 1L]]
}
root <- normalizePath(option("--server", "/refitting-service"), mustWork = TRUE)
output <- normalizePath(option("--output", "/tmp/msigportal-refitting-regression"), mustWork = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
setwd(root)
if (!is.null(option("--baseline"))) stop("Refitting preflight has no numerical baseline; fitting integration requires real service fixtures")
cases <- list(source = list(target = "refitting.R", scope = "deployment-load"))
for (fn in c("run_sbs_refitting", "run_dbs_refitting")) {
  cases[[paste0(fn, ".missing-files")]] <- list(target = fn, scope = "validation", pattern = "Missing user files", genome = "hg19")
  cases[[paste0(fn, ".invalid-genome")]] <- list(target = fn, scope = "validation", pattern = "arg.*should be one of", genome = "invalid")
}
selected <- names(cases)
if (!is.null(option("--case"))) selected <- selected[grepl(option("--case"), selected)]
if (!length(selected)) stop("No refitting cases selected")
if ("--list" %in% argv) { writeLines(selected); quit(status = 0L) }
writeLines(capture.output(sessionInfo()), file.path(output, "session-info.txt"))
write.csv(as.data.frame(installed.packages()[, c("Package", "Version", "Built")]), file.path(output, "packages.csv"), row.names = FALSE)
write.csv(data.frame(file = "refitting.R", md5 = unname(tools::md5sum("refitting.R"))), file.path(output, "sources.csv"), row.names = FALSE)
app <- new.env(parent = globalenv())
load_error <- NULL
tryCatch(source("refitting.R", local = app), error = function(e) load_error <<- conditionMessage(e))
results <- lapply(selected, function(name) {
  spec <- cases[[name]]
  log <- file.path(output, paste0(name, ".log"))
  message <- ""
  status <- "PASS"
  tryCatch({
    if (!is.null(load_error)) stop(load_error)
    if (name == "source") {
      stopifnot(is.function(app$run_sbs_refitting), is.function(app$run_dbs_refitting))
    } else {
      missing <- file.path(output, "intentionally-missing-input")
      caught <- NULL
      text <- capture.output(tryCatch(do.call(get(spec$target, app), list(maf_file = missing,
        genomic_file = missing, clinical_file = missing, output_dir = output,
        common_files_dir = output, genome = spec$genome)), error = function(e) caught <<- e))
      writeLines(text, log)
      if (!inherits(caught, "error") || !grepl(spec$pattern, conditionMessage(caught))) stop("Expected validation error was not raised")
    }
  }, error = function(e) { status <<- "FAIL"; message <<- conditionMessage(e) })
  cat(name, status, message, "\n")
  data.frame(case = name, target = spec$target, scope = spec$scope, status = status, message = message)
})
summary <- do.call(rbind, results)
write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
cat(sum(summary$status == "PASS"), "passed,", sum(summary$status == "FAIL"), "failed\n")
quit(status = as.integer(any(summary$status == "FAIL")))
