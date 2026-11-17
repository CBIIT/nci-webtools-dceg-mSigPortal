#!/usr/bin/env Rscript
# Checks the harness assertions and snapshot behavior without plotting packages.
script <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
env <- new.env(parent = globalenv())
sys.source(file.path(dirname(normalizePath(script)), "helpers.R"), env)
with(env, {
  equal(c(1, 2), c(1, 2 + 1e-10))
  expect_error(equal(c(1, 2), c(1, 3)), "difference")
  expect_error(assert_success(list(output = list(uncaughtError = "hidden error"))), "hidden error")
  expect_error(assert_success(list(landscapeError = "caught by wrapper")), "caught by wrapper")
  assert_success(list(error = NULL, data = data.frame(error = character())))
  fixture <- list(work = "/tmp/test-work", output = "/tmp/test-output")
  x <- list(path = "/tmp/test-work/output.svg", timestamp = "yesterday", uid = "random-a",
    data = data.frame(value = c(1, 2), group = factor(c("A", "B"))))
  y <- x; y$timestamp <- "today"; y$uid <- "random-b"
  equal(canonical(x, fixture), canonical(y, fixture))
  expect(identical(canonical(x, fixture)$path, "<WORK>/output.svg"))
  y$data$value[1] <- 3
  expect(!isTRUE(all.equal(canonical(x, fixture), canonical(y, fixture))))
  y <- x; levels(y$data$group) <- c("B", "A")
  expect(!isTRUE(all.equal(canonical(x, fixture), canonical(y, fixture))))
  widget_a <- list(attrs = list(abc123 = list(value = 42)), visdat = list(abc123 = "source"), cur_data = "abc123")
  widget_b <- list(attrs = list(def456 = list(value = 42)), visdat = list(def456 = "source"), cur_data = "def456")
  equal(canonical(widget_a, fixture), canonical(widget_b, fixture))
  widget_b$attrs$def456$value <- 43
  expect(!isTRUE(all.equal(canonical(widget_a, fixture), canonical(widget_b, fixture))))
  root <- tempfile("r-harness-artifacts-")
  artifact_fixture <- list(output = root, work = file.path(root, "work"), artifacts = file.path(root, "artifacts"))
  dir.create(artifact_fixture$work, recursive = TRUE)
  dir.create(artifact_fixture$artifacts)
  con <- gzfile(file.path(artifact_fixture$work, "table.txt.gz"), "wt")
  write.table(data.frame(Sample = c("A", "B"), Mutations = c(2.5, 3.5)), con, sep = "\t", row.names = FALSE)
  close(con)
  artifacts <- collect_artifacts(artifact_fixture)
  equal(artifacts$tables[["table.txt.gz"]]$data$Mutations, c(2.5, 3.5))
  cat("Harness assertion and snapshot checks passed\n")
})
