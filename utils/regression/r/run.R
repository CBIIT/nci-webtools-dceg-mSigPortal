#!/usr/bin/env Rscript
# Base-R runner. Every case gets a fresh R process; no package installation occurs.
argv <- commandArgs(trailingOnly = TRUE)
option <- function(key, default = NULL) {
  i <- match(key, argv)
  if (is.na(i)) return(default)
  if (i == length(argv)) stop("Missing value for ", key)
  argv[[i + 1L]]
}
script <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
script <- normalizePath(script)
suite_dir <- dirname(script)
server_dir <- normalizePath(option("--server", getwd()), mustWork = TRUE)
output_dir <- normalizePath(option("--output", "/tmp/msigportal-r-regression"), mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
baseline_dir <- option("--baseline")
if (!is.null(baseline_dir)) baseline_dir <- normalizePath(baseline_dir, mustWork = TRUE)
setwd(server_dir)

definitions <- function(file) {
  Filter(function(x) is.call(x) && as.character(x[[1]]) %in% c("<-", "=") &&
    length(x) == 3L && is.call(x[[3]]) && identical(x[[3]][[1]], as.name("function")), parse(file))
}
files <- list.files("services/R", pattern = "[.]R$", full.names = TRUE)
# Refitting runs in its own service and suite, not in the web backend.
files <- files[basename(files) != "refittingWrapper.R"]
if (!length(files)) stop("No services/R/*.R files found under --server: ", server_dir)
inventory <- do.call(rbind, lapply(files, function(file) {
  defs <- definitions(file)
  data.frame(file = basename(file), function_name = vapply(defs, function(x) as.character(x[[2]]), ""))
}))
inventory$key <- paste(inventory$file, inventory$function_name, sep = "::")

# Test helpers are sourced into a separate environment from application globals.
harness <- new.env(parent = globalenv())
sys.source(file.path(suite_dir, "helpers.R"), harness)
sys.source(file.path(suite_dir, "cases.R"), harness)
cases <- harness$make_cases()
unknown_targets <- setdiff(vapply(cases, function(x) x$target, ""), inventory$key)
if (length(unknown_targets)) stop("Stale case targets: ", paste(unknown_targets, collapse = ", "))
unregistered <- setdiff(inventory$key, vapply(cases, function(x) x$target, ""))
if (length(unregistered)) {
  writeLines(unregistered, file.path(output_dir, "unregistered-functions.txt"))
  stop("Functions need explicit test cases: ", paste(unregistered, collapse = ", "))
}
if ("--list" %in% argv) {
  print(data.frame(case = names(cases), target = vapply(cases, function(x) x$target, ""),
    scope = vapply(cases, function(x) x$scope, "")), row.names = FALSE)
  quit(status = 0L)
}

worker <- option("--worker")
if (!is.null(worker)) {
  case <- cases[[worker]]
  if (is.null(case)) stop("Unknown case: ", worker)
  result <- list(case = worker, target = case$target, scope = case$scope,
    status = "FAIL", message = "", warnings = character(), hits = character())
  hits <- new.env(parent = emptyenv())
  instrument <- function(file) {
    if (!file.exists(file)) return(invisible(NULL))
    if (!basename(file) %in% inventory$file) return(invisible(NULL))
    for (expr in definitions(file)) {
      name <- as.character(expr[[2]])
      if (!exists(name, globalenv(), inherits = FALSE)) next
      fun <- get(name, globalenv())
      if (inherits(fun, "functionWithTrace")) next
      key <- paste(basename(file), name, sep = "::")
      tracer <- substitute(assign(KEY, TRUE, envir = .regression_hits), list(KEY = key))
      invisible(capture.output(trace(name, tracer = tracer, print = FALSE, where = globalenv())))
    }
  }
  assign(".regression_hits", hits, globalenv())
  # Preserve production source() behavior, then instrument newly defined helpers.
  assign("source", function(file, ...) {
    value <- base::source(file, ...)
    instrument(file)
    invisible(value)
  }, globalenv())
  started <- proc.time()[["elapsed"]]
  tryCatch(withCallingHandlers({
    set.seed(5177184)
    options(warn = 1, width = 120, lifecycle_verbosity = "warning")
    fixture <- harness$make_fixture(output_dir, server_dir)
    options(device = function(...) grDevices::pdf(file.path(fixture$work, "implicit-device.pdf"), ...))
    harness$install_transport(fixture)
    file <- file.path("services/R", strsplit(case$target, "::", fixed = TRUE)[[1]][1])
    # Shared plotting helpers normally execute after the visualization wrapper
    # attaches its graphics dependencies and sources utils.R.
    if (basename(file) == "Sigvisualfunc.R") source("services/R/visualizeWrapper.R")
    source(file)
    if (case$load_helpers) source("services/R/Sigvisualfunc.R")
    value <- case$run(fixture)
    harness$assert_success(value)
    while (dev.cur() > 1L) dev.off()
    snapshot <- harness$canonical(value, fixture)
    artifacts <- harness$collect_artifacts(fixture)
    snapshot <- list(value = snapshot, tables = artifacts$tables)
    saveRDS(snapshot, file.path(output_dir, "snapshot.rds"), version = 2)
    if (!is.null(baseline_dir)) {
      baseline_file <- file.path(baseline_dir, worker, "result.rds")
      harness$expect(file.exists(baseline_file), "Baseline case result missing")
      previous <- readRDS(baseline_file)
      harness$expect(identical(previous$status, "PASS"), "Baseline case did not pass")
      snapshot_file <- file.path(baseline_dir, worker, "snapshot.rds")
      harness$expect(file.exists(snapshot_file), "Baseline snapshot missing")
      difference <- all.equal(readRDS(snapshot_file), snapshot, tolerance = 1e-6, check.attributes = TRUE)
      harness$expect(isTRUE(difference), paste("Baseline differs:", paste(difference, collapse = "; ")))
      old_visuals <- file.path(baseline_dir, worker, "artifacts")
      result$visual_changes <- harness$visual_changes(old_visuals, fixture$artifacts)
    }
    result$status <- "PASS"
  }, warning = function(w) {
    result$warnings <<- c(result$warnings, conditionMessage(w))
  }), error = function(e) {
    result$message <<- conditionMessage(e)
    cat("\nFAIL: ", conditionMessage(e), "\n", sep = "")
    if (!is.null(conditionCall(e))) print(conditionCall(e))
  })
  while (dev.cur() > 1L) try(dev.off(), silent = TRUE)
  result$hits <- ls(hits)
  result$seconds <- proc.time()[["elapsed"]] - started
  if (result$status == "PASS" && !case$target %in% result$hits) {
    result$status <- "FAIL"
    result$message <- "Case never invoked its declared target"
  }
  saveRDS(result, file.path(output_dir, "result.rds"), version = 2)
  quit(status = if (result$status == "PASS") 0L else 1L)
}

selected <- names(cases)
pattern <- option("--case")
if (!is.null(pattern)) selected <- selected[grepl(pattern, selected)]
if (!length(selected)) stop("No cases match --case")
if (any(file.exists(file.path(output_dir, selected)))) stop("Output already contains selected cases; use a fresh output directory")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session-info.txt"))
write.csv(as.data.frame(installed.packages()[, c("Package", "Version", "Built")]),
  file.path(output_dir, "packages.csv"), row.names = FALSE)
write.csv(data.frame(file = files, md5 = unname(tools::md5sum(files))),
  file.path(output_dir, "sources.csv"), row.names = FALSE)
results <- list()
for (name in selected) {
  cat(sprintf("[%d/%d] %s ... ", length(results) + 1L, length(selected), name))
  flush.console()
  case_dir <- file.path(output_dir, name)
  dir.create(case_dir)
  args <- c("--vanilla", script, "--server", server_dir, "--output", case_dir, "--worker", name)
  if (!is.null(baseline_dir)) args <- c(args, "--baseline", baseline_dir)
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(args),
    stdout = file.path(case_dir, "run.log"), stderr = file.path(case_dir, "run.log"),
    timeout = as.numeric(option("--timeout", "600")))
  path <- file.path(case_dir, "result.rds")
  result <- if (file.exists(path)) readRDS(path) else list(case = name, target = cases[[name]]$target,
    scope = cases[[name]]$scope, status = "FAIL", message = paste("Worker exited without report; exit", status),
    warnings = character(), hits = character(), seconds = NA_real_)
  if (status != 0L && result$status == "PASS") {
    result$status <- "FAIL"
    result$message <- paste("Worker exit", status)
  }
  results[[name]] <- result
  cat(result$status, if (nzchar(result$message)) paste0(": ", result$message) else "", "\n", sep = "")
}
summary <- do.call(rbind, lapply(results, function(x) data.frame(case = x$case, target = x$target,
  scope = x$scope, status = x$status, seconds = x$seconds, warnings = length(x$warnings), message = x$message)))
write.csv(summary, file.path(output_dir, "summary.csv"), row.names = FALSE)
inventory$status <- vapply(inventory$key, function(key) {
  relevant <- Filter(function(x) key %in% x$hits, results)
  if (any(vapply(relevant, function(x) x$status == "PASS", TRUE))) "EXERCISED"
  else if (length(relevant)) "INVOKED_IN_FAILED_CASE" else "NOT_EXERCISED"
}, "")
inventory$cases <- vapply(inventory$key, function(key) paste(names(Filter(function(x) key %in% x$hits, results)), collapse = ";"), "")
write.csv(inventory, file.path(output_dir, "coverage.csv"), row.names = FALSE)
writeLines(unlist(lapply(results, function(x) c(paste0(x$case, ":"), unique(x$warnings)))),
  file.path(output_dir, "warnings.txt"))
writeLines(as.character(unlist(lapply(results, function(x) if (length(x$visual_changes)) paste(x$case, x$visual_changes)))),
  file.path(output_dir, "visual-changes.txt"))
failed <- sum(summary$status != "PASS")
cat(sprintf("\n%d passed, %d failed. %d/%d functions exercised successfully.\nReports: %s\n",
  sum(summary$status == "PASS"), failed, sum(inventory$status == "EXERCISED"), nrow(inventory), output_dir))
quit(status = if (failed) 1L else 0L)
