expect <- function(ok, message = "Assertion failed") {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  invisible(TRUE)
}
equal <- function(actual, expected, tolerance = 1e-8) {
  diff <- all.equal(actual, expected, tolerance = tolerance, check.attributes = FALSE)
  expect(isTRUE(diff), paste(diff, collapse = "; "))
  invisible(actual)
}
expect_error <- function(expr, pattern, class = "error") {
  caught <- tryCatch({ force(expr); NULL }, error = identity)
  expect(inherits(caught, class), paste("Expected error class", class))
  expect(grepl(pattern, conditionMessage(caught)), paste("Unexpected error:", conditionMessage(caught)))
  TRUE
}
assert_success <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    error_fields <- grep("(^error$|^uncaughtError$|Error$)", names(x), value = TRUE)
    for (key in error_fields) expect(is.null(x[[key]]) || length(x[[key]]) == 0L,
      paste(key, paste(unlist(x[[key]]), collapse = "; ")))
    if (identical(x$success, FALSE)) stop(paste(x$error), call. = FALSE)
    for (item in x) assert_success(item)
  }
  invisible(x)
}

make_fixture <- function(output, server) {
  work <- file.path(output, "work")
  artifacts <- file.path(output, "artifacts")
  dir.create(file.path(work, "output", "plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(artifacts, showWarnings = FALSE)
  expect(file.symlink(file.path(server, "services"), file.path(work, "services")), "Cannot link application services into fixture work directory")
  sbs <- unlist(lapply(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"), function(type)
    as.vector(outer(c("A", "C", "G", "T"), c("A", "C", "G", "T"), function(a, b) paste0(a, "[", type, "]", b)))))
  dbs <- readLines(file.path(suite_dir, "dbs78.txt"))
  id_groups <- c("1:Del:C", "1:Del:T", "1:Ins:C", "1:Ins:T", paste0(2:5, ":Del:R"),
    paste0(2:5, ":Ins:R"), paste0(2:5, ":Del:M"))
  ids <- unlist(Map(function(group, n) paste0(group, ":", seq_len(n) - as.integer(!grepl(":M$", group))),
    id_groups, c(rep(6L, 12), 1L, 2L, 3L, 5L)))
  profile_types <- list(SBS96 = sbs, DBS78 = dbs, ID83 = ids,
    SBS192 = c(paste0("T:", sbs), paste0("U:", sbs)))
  samples <- sprintf("Sample%02d", 1:24)
  sig_names <- c("SBS1", "SBS2", "SBS3")
  W <- matrix(runif(96 * 3, 0.1, 3), ncol = 3, dimnames = list(NULL, sig_names))
  W <- sweep(W, 2L, colSums(W), "/")
  H <- matrix(round(runif(24 * 3, 200, 2500)), ncol = 3, dimnames = list(samples, sig_names))
  G <- W %*% t(H)
  signature <- data.frame(MutationType = sbs, W, check.names = FALSE)
  activities <- data.frame(Sample = samples, H, check.names = FALSE)
  matrix_data <- data.frame(MutationType = sbs, G, check.names = FALSE)
  cancers <- rep(c("CancerA", "CancerB"), each = 12)
  exposure <- expand.grid(Sample = samples, Signature_name = sig_names, stringsAsFactors = FALSE)
  exposure$Exposure <- as.vector(H)
  exposure$Cancer_Type <- cancers[match(exposure$Sample, samples)]
  exposure$Study <- "PCAWG"
  exposure$Dataset <- "WGS"
  exposure$Signature_set_name <- "COSMIC v3 Signatures (SBS)"
  ref <- expand.grid(MutationType = sbs, Signature_name = sig_names, stringsAsFactors = FALSE)
  ref$Contribution <- as.vector(W)
  ref$Profile <- "SBS96"
  ref$Signature_set_name <- "COSMIC v3 Signatures (SBS)"
  ref$Dataset <- "WGS"
  ref$Source <- "Reference"
  seqdata <- expand.grid(MutationType = sbs, Sample = samples, stringsAsFactors = FALSE)
  seqdata$Mutations <- as.vector(G)
  seqdata$Cancer_Type <- cancers[match(seqdata$Sample, samples)]
  seqdata$Study <- "PCAWG"
  seqdata$Dataset <- "WGS"
  seqdata$Profile <- "SBS96"
  age <- data.frame(Sample = samples, variable_name = "Age", variable_value = as.character(25 + seq_along(samples) * 1.7), variable_value_type = "numeric")
  group <- data.frame(Sample = samples, variable_name = "Group", variable_value = rep(c("A", "B"), 12), variable_value_type = "character")
  vardata <- rbind(age, group)
  vardata$data_source <- "clinical"
  vardata$data_type <- "demographic"
  vardata$Cancer_Type <- cancers[match(vardata$Sample, samples)]
  assoc_data <- data.frame(x = 1:48, y = (1:48) * 0.7 + sin(1:48) * 4,
    category = rep(c("A", "B"), 24), category3 = rep(c("A", "B", "C"), 16),
    group = rep(c("SBS1", "SBS2"), each = 24))
  config <- list(wd = work, savePath = paste0(work, "/output/plots/"),
    prefix = "fixture/", bucket = "regression-offline", id = "", localData = TRUE)
  write.table(matrix_data, file.path(work, "matrix.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(signature, file.path(work, "signature.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(activities, file.path(work, "exposure.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.csv(data.frame(SAMPLE_ID = samples, Signature = "SBS1", H_Burden = H[, 1]),
    file.path(work, "job_H_Burden_est.csv"), row.names = FALSE)
  list(work = work, artifacts = artifacts, output = output, server = server, config = config,
    samples = samples, signature = signature, activities = activities, matrix = matrix_data,
    exposure = exposure, ref = ref, seqdata = seqdata, vardata = vardata, association = assoc_data,
    profile_types = profile_types, W = W, H = H, G = G, set = unique(ref$Signature_set_name),
    matrix_file = file.path(work, "matrix.tsv"))
}

install_transport <- function(F) {
  # Only the external storage boundary is replaced; analysis code remains real.
  loader <- function(object, bucket, envir = parent.frame(), ...) {
    expect(identical(bucket, "regression-offline"), "Unexpected S3 bucket")
    name <- if (grepl("subset_files", object)) "seqmatrix_refdata_subset_files"
      else if (grepl("seqmatrix_refdata_info", object)) "seqmatrix_refdata_info"
      else if (grepl("signature_refsets", object)) "signature_refsets"
      else if (grepl("exposure_refdata", object)) "exposure_refdata"
      else if (grepl("vardata", object)) "vardata_refdata"
      else if (grepl("content_data_all", object)) "content_data_all"
      else if (grepl("seqmatrix", object)) "seqmatrix_refdata"
      else stop("Unregistered S3 fixture: ", object)
    value <- switch(name, signature_refsets = F$ref, exposure_refdata = F$exposure,
      seqmatrix_refdata = if (grepl("cancerA", object)) F$seqdata[F$seqdata$Cancer_Type == "CancerA", ]
        else if (grepl("cancerB", object)) F$seqdata[F$seqdata$Cancer_Type == "CancerB", ] else F$seqdata,
      vardata_refdata = F$vardata,
      seqmatrix_refdata_subset_files = data.frame(Study = "PCAWG", Dataset = "WGS",
        Cancer_Type = c("CancerA", "CancerB"), file = c("cancerA_seqmatrix.RData", "cancerB_seqmatrix.RData")),
      seqmatrix_refdata_info = data.frame(Study = "PCAWG", Dataset = "WGS", Cancer_Type = "CancerA",
        Sample = F$samples[1:12], Profile = "SBS96", Profile_Type = "SBS", Matrix_Size = 96L, Path = "sample.svg"),
      content_data_all = data.frame(Study = "PCAWG@CancerA", Pattern = c("TCW>TTW", "NCG>NTG"), N1 = c(0.2, 0.3)))
    assign(name, value, envir)
    invisible(name)
  }
  assign("s3load", loader, globalenv())
  # Qualified calls cannot be shadowed by a global binding. Change only this
  # read-only existence probe inside the disposable worker process.
  if (requireNamespace("aws.s3", quietly = TRUE)) {
    ns <- asNamespace("aws.s3")
    unlockBinding("object_exists", ns)
    assign("object_exists", function(...) FALSE, ns)
    lockBinding("object_exists", ns)
  }
  assign("config", F$config, globalenv())
  assign("study", "PCAWG", globalenv())
  assign("cancer_type", "CancerA", globalenv())
}

profile <- function(F, kind = "SBS96") {
  types <- F$profile_types[[kind]]
  data.frame(MutationType = types, Value = seq_along(types) %% 17 + 1)
}
api_args <- function(F) {
  list(exposureData = dplyr::transmute(F$exposure, sample = Sample, signatureName = Signature_name,
      exposure = Exposure, cancer = Cancer_Type),
    signatureData = dplyr::transmute(F$ref, mutationType = MutationType, signatureName = Signature_name,
      contribution = Contribution, profile = "SBS", matrix = 96L),
    seqmatrixData = dplyr::transmute(F$seqdata, mutationType = MutationType, sample = Sample,
      mutations = Mutations, cancer = Cancer_Type, profile = "SBS", matrix = 96L))
}
public_args <- function(F, relative = FALSE) {
  list(study = "PCAWG", cancerType = "CancerA", experimentalStrategy = "WGS", strategy = "WGS",
    profileType = "SBS", matrixSize = "96", profileName = "SBS96", signatureSet = F$set,
    signatureSetName = F$set, rsSet = F$set, rsSet1 = F$set, rsSet2 = F$set,
    signatureSource = "Reference", signatureName = "SBS1", signatureName1 = "SBS1", signatureName2 = "SBS2",
    sampleName = F$samples[1], sampleName1 = F$samples[1], sampleName2 = F$samples[2],
    userSample = F$samples[1], publicSample = F$samples[2],
    compare = "SBS1", pattern = "TCW>TTW", proportion = 0.02,
    matrixFile = if (relative) "matrix.tsv" else F$matrix_file)
}
association_args <- function(F) {
  variable <- list(source = "clinical", type = "demographic", name = "Age", log2 = FALSE)
  list(study = "PCAWG", strategy = "WGS", cancer = "CancerA", rsSet = F$set,
    source = "clinical", type = "demographic", assocName = "Group", expName = "Signature_exposure",
    associationVar = variable, exposureVar = list(name = "Signature_exposure", log2 = FALSE),
    associationVars = data.frame(source = "clinical", type = "demographic", name = "Age", log2 = "FALSE"),
    signature = "SBS1", testType = "parametric", xlab = "Age", ylab = "Exposure")
}
render <- function(plot, F, name = "returned-plot") {
  expect(inherits(plot, c("ggplot", "grob", "gtable", "patchwork")), "Function did not return a plot")
  path <- file.path(F$artifacts, paste0(name, ".svg"))
  ggplot2::ggsave(path, plot = plot, device = svglite::svglite, width = 12, height = 8, limitsize = FALSE)
  expect(file.exists(path) && file.info(path)$size > 100, "Plot rendering produced an empty file")
  if (inherits(plot, "ggplot") && !inherits(plot, "patchwork")) {
    return(list(layers = ggplot2::ggplot_build(plot)$data))
  }
  list(rendered = TRUE)
}
capture_plot_data <- function(run) {
  # Some helpers return a gtable-backed plot. Observe their real input to the
  # graphics build so numerical assertions can precede the rendering check.
  had_binding <- exists("ggplot_build", globalenv(), inherits = FALSE)
  previous <- if (had_binding) get("ggplot_build", globalenv()) else NULL
  on.exit({
    if (had_binding) assign("ggplot_build", previous, globalenv()) else
      rm(list = "ggplot_build", envir = globalenv())
  })
  captured <- list()
  assign("ggplot_build", function(plot, ...) {
    captured[[length(captured) + 1L]] <<- plot$data
    ggplot2::ggplot_build(plot, ...)
  }, globalenv())
  plot <- run()
  expect(length(captured) > 0L, "No plot data reached the graphics build")
  list(plot = plot, data = captured[[1]])
}
checked_paths <- function(value, F) {
  assert_success(value)
  if (!is.list(value) || is.data.frame(value)) return(value)
  paths <- unlist(Filter(is.character,
    value[grepl("(Path$|Data$|^pca[123]$|^heatmap$|^path$)", names(value))]))
  for (path in paths) expect(file.exists(path) && file.info(path)$size > 0, paste("Missing/empty output:", path))
  value
}

canonical <- function(x, F, depth = 0L) {
  if (depth > 30L) stop("Snapshot nesting too deep")
  if (is.factor(x)) return(list(values = as.character(x), levels = levels(x)))
  if (is.data.frame(x)) return(list(columns = names(x), rows = nrow(x),
    data = lapply(x, function(item) canonical(item, F, depth + 1L))))
  if (is.matrix(x) || is.array(x)) return(list(dim = dim(x), dimnames = dimnames(x), values = as.vector(x)))
  if (is.list(x)) {
    # plotly_json() keys its unevaluated source-data metadata by process-specific
    # hashes. Preserve the metadata and its references, with stable local IDs.
    if (is.list(x$attrs) && is.list(x$visdat)) {
      ids <- unique(c(names(x$attrs), names(x$visdat)))
      ids <- ids[grepl("^[0-9a-f]+$", ids)]
      stable <- setNames(paste0("plotly-data-", seq_along(ids)), ids)
      for (slot in c("attrs", "visdat")) {
        original <- names(x[[slot]])
        found <- original %in% ids
        original[found] <- unname(stable[original[found]])
        names(x[[slot]]) <- original
      }
      if (is.character(x$cur_data) && all(x$cur_data %in% ids)) x$cur_data <- unname(stable[x$cur_data])
    }
    # Widget element IDs and trace UIDs are random rendering metadata.
    if (!is.null(names(x))) x <- x[!names(x) %in% c("uid", "elementId", "lastModified", "timestamp", "fileSize")]
    return(lapply(x, function(item) canonical(item, F, depth + 1L)))
  }
  if (is.character(x)) {
    x <- gsub(F$work, "<WORK>", x, fixed = TRUE)
    x <- gsub(F$output, "<OUTPUT>", x, fixed = TRUE)
  }
  if (is.language(x)) return(paste(deparse(x), collapse = " "))
  if (is.environment(x) || is.function(x) || isS4(x)) return(class(x))
  x
}
collect_artifacts <- function(F) {
  paths <- list.files(F$work, recursive = TRUE, full.names = TRUE)
  paths <- paths[!grepl("/services/", paths, fixed = TRUE)]
  # Fixture input files are intentionally excluded from outputs/baselines.
  paths <- paths[!basename(paths) %in% c("matrix.tsv", "signature.tsv", "exposure.tsv", "job_H_Burden_est.csv")]
  tables <- list()
  for (path in paths) {
    rel <- substring(path, nchar(F$work) + 2L)
    dst <- file.path(F$artifacts, rel)
    dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
    file.copy(path, dst, overwrite = TRUE)
    if (grepl("[.](txt|tsv|csv)([.]gz)?$", path)) {
      table <- read.table(path, header = TRUE, sep = if (grepl("[.]csv([.]gz)?$", path)) "," else "\t",
        quote = "\"", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
      tables[[rel]] <- canonical(table, F)
    }
    if (grepl("[.]svg$", path)) expect(any(grepl("<svg", readLines(path, warn = FALSE))), paste("Invalid SVG:", path))
  }
  list(tables = tables)
}
visual_changes <- function(before, after) {
  paths <- sort(unique(c(list.files(before, pattern = "[.](svg|pdf|png)$", recursive = TRUE),
    list.files(after, pattern = "[.](svg|pdf|png)$", recursive = TRUE))))
  Filter(function(p) !file.exists(file.path(before, p)) || !file.exists(file.path(after, p)) ||
    !identical(unname(tools::md5sum(file.path(before, p))), unname(tools::md5sum(file.path(after, p)))), paths)
}
