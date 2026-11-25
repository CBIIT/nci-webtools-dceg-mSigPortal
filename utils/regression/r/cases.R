# Each top-level backend function must have an explicit target below. Additional
# cases cover important branches. A new function without a case fails discovery.
make_cases <- function() {
  cases <- list()
  add <- function(file, name, run, suffix = NULL, helpers = FALSE, scope = "functional") {
    id <- paste(sub("[.]R$", "", file), name, suffix, sep = ".")
    id <- sub("[.]$", "", id)
    expect(is.null(cases[[id]]), paste("Duplicate case", id))
    cases[[id]] <<- list(target = paste(file, name, sep = "::"), run = run,
      load_helpers = helpers, scope = scope)
  }
  sig <- function(name, run, suffix = NULL, scope = "functional")
    add("Sigvisualfunc.R", name, run, suffix, scope = scope)
  invoke <- function(name, ...) do.call(get(name, envir = globalenv()), list(...))
  plot_case <- function(name, args, suffix = NULL) {
    force(name); force(args)
    sig(name, function(F) render(do.call(get(name, globalenv()), args(F)), F), suffix)
  }
  wide_exposure <- function(F) { x <- F$activities; names(x)[1] <- "Samples"; x }
  motif <- function(F) dplyr::mutate(F$seqdata, Type = substr(MutationType, 3, 5),
    SubType1 = substr(MutationType, 1, 1), SubType2 = substr(MutationType, 7, 7))
  similarities <- function(F) calculate_similarities(F$matrix, F$signature, F$activities)

  # Independent numerical oracles, including perfect reconstruction.
  sig("cos_sim", function(F) {
    equal(cos_sim(c(1, 0), c(0, 1)), 0)
    equal(cos_sim(c(1, 2), c(1, 2)), 1)
    equal(cos_sim(c(1, 0), c(1, 1)), 1 / sqrt(2))
  })
  sig("calculate_similarities", function(F) {
    x <- similarities(F)
    equal(x$Cosine_similarity, rep(1, 24)); equal(x$L1_Norm, rep(0, 24))
    equal(x$L2_Norm, rep(0, 24)); equal(x$Total_Mutations, colSums(F$G), tolerance = 1e-5)
    shifted <- F$matrix
    shifted[1, -1] <- shifted[1, -1] + 17
    perturbed <- calculate_similarities(shifted, F$signature, F$activities)
    equal(perturbed$L1_Norm, rep(17, 24)); equal(perturbed$L2_Norm, rep(17, 24))
    oracle <- vapply(seq_len(24), function(j) {
      actual <- shifted[[j + 1L]]; reconstructed <- F$G[, j]
      round(sum(actual * reconstructed) / sqrt(sum(actual^2) * sum(reconstructed^2)), 3)
    }, 0)
    equal(perturbed$Cosine_similarity, oracle)
    list(perfect = x, perturbed = perturbed)
  })
  sig("probabilities", function(F) {
    H <- F$activities; names(H)[1] <- "Samples"
    x <- probabilities(F$signature, H)
    equal(rowSums(as.matrix(x[, -(1:2)])), rep(1, 24 * 96))
    x
  })
  for (name in c("cos_sim_df", "cos_sim_df_old")) local({
    target <- name
    sig(target, function(F) {
      x <- invoke(target, F$signature, F$signature, output_matrix = TRUE)
      equal(diag(as.matrix(x)), rep(1, 3)); x
    })
  })
  sig("cos_sim_df", function(F) {
    x <- cos_sim_df(F$signature[, 1:2], F$signature[, 1:2])
    expect(nrow(x) == 1 && ncol(x) == 2); equal(x[[2]], 1); x
  }, "single-column")
  sig("signature_sum_operation", function(F) {
    x <- signature_sum_operation(F$ref, F$set, "0.25*SBS1;0.75*SBS2")
    expected <- 0.25 * F$W[, 1] + 0.75 * F$W[, 2]
    equal(x$Contribution[match(F$signature$MutationType, x$MutationType)], expected); x
  })
  sig("iupac", function(F) {
    equal(sort(iupac("N")), c("A", "C", "G", "T")); equal(iupac("Y"), c("C", "T")); iupac("R")
  })
  sig("genome2size", function(F) {
    equal(genome2size("GRCh37"), 3101976562 / 1e6)
    equal(genome2size("GRCh38"), 3217346917 / 1e6); TRUE
  })
  sig("rowAny", function(F) equal(rowAny(matrix(c(0, 0, 1, 0), nrow = 2, byrow = TRUE)), c(FALSE, TRUE)))
  sig("change_data_type", function(F) {
    equal(change_data_type(c("1", "2")), c(1, 2)); equal(change_data_type(c(1, 2)), c("1", "2")); TRUE
  })
  sig("restore_var_names", function(F) equal(restore_var_names(c("X7qOther", "Age", "(Intercept)"),
    c(X7q = "7q")), c("7qOther", "Age", "(Intercept)")))
  sig("finalize_group_result", function(F) {
    x <- finalize_group_result(data.frame(p.value = c(0.04, NA, 0.01)), "parametric")
    equal(x$fdr, c(0.02, 0.04)); x
  })
  sig("finalize_group_result", function(F) {
    x <- finalize_group_result(data.frame(bf10 = c(2, NA, 10)), "bayes")
    equal(x$bf10, c(10, 2)); expect(!"fdr" %in% names(x)); x
  }, "bayes")
  sig("validate_vardf", function(F) {
    x <- validate_vardf(data.frame(Age = 1:12, Group = rep(c("A", "B"), 6)))
    expect(is.numeric(x$Age)); expect(is.factor(x$Group)); x
  })
  sig("multivariable_inputs", function(F) {
    x <- multivariable_inputs(F$vardata, list(list(source = "clinical", type = "demographic", name = "Age", log2 = "FALSE")))
    equal(x$Age, as.numeric(F$vardata$variable_value[F$vardata$variable_name == "Age"])); x
  })
  sig("resolve_conflicts", function(F) { resolve_conflicts(); TRUE })

  for (kind in c("SBS96", "SBS192", "DBS78", "ID83")) local({
    type <- kind
    sig("profile_format_df", function(F) {
      x <- profile_format_df(profile(F, type), factortype = TRUE)
      expect(nrow(x) == nrow(profile(F, type))); expect(!anyNA(x$Type));
      equal(sum(x$Value), sum(profile(F, type)$Value)); x
    }, type)
    plot_case("plot_profile_logo", function(F) list(profile(F, type)), type)
    plot_case("plot_compare_profiles_diff", function(F) {
      x <- profile(F, type); y <- x; y$Value <- rev(y$Value)
      list(profile1 = x, profile2 = y)
    }, type)
    plot_case("plot_compare_profiles_diff2", function(F) {
      x <- profile(F, type); x$Second <- rev(x$Value); list(profile = x)
    }, type)
    sig("plot_compare_profiles_diff2", function(F) {
      x <- profile(F, type)
      x$Scaled <- 100 * x$Value
      x$Reversed <- rev(x$Value)
      captured <- capture_plot_data(function() plot_compare_profiles_diff2(x,
        profile_names = c("Original", "Scaled", "Reversed"), condensed = TRUE))
      data <- captured$data
      equal(as.numeric(tapply(data$value, data$name, sum)), rep(1, 3))
      equal(data$value[data$name == "Original"], data$value[data$name == "Scaled"])
      expect(nrow(data) == 3L * nrow(x))
      if (type == "SBS192") equal(sort(unique(as.character(data$Strand))), c("T", "U"))
      else expect(!"Strand" %in% names(data))
      list(data = data, rendered = render(captured$plot, F))
    }, paste0(type, ".normalization"))
  })
  plot_case("plot_sbs_96_profile", function(F) list(data = profile(F), samplename = "Fixture"))
  plot_case("plot_dbs_78_profile", function(F) list(data = profile(F, "DBS78"), samplename = "Fixture"))
  plot_case("plot_dbs_78_profile_old", function(F) list(data = profile(F, "DBS78"), samplename = "Fixture"))
  for (name in c("plot_profile_id29", "plot_profile_id29_logo")) local({
    target <- name
    plot_case(target, function(F) list(profile = data.frame(MutationType = id29_col1, Value = seq_along(id29_col1))))
  })
  plot_case("plot_profile_logo_cn48", function(F) list(profile = data.frame(
    MutationType = cn48_mutationtype, Value = seq_along(cn48_mutationtype))))
  for (name in c("plot_profile_dbs78_signal", "plot_profile_dbs78_signal_logo")) local({
    target <- name; plot_case(target, function(F) list(profile = profile(F, "DBS78")))
  })
  plot_case("plot_compare_profiles_diff_strand", function(F) {
    x <- profile_format_df(profile(F, "SBS192"), factortype = TRUE)
    list(profile1 = x, profile2 = x)
  })
  sig("plot_compare_profiles_diff_strand", function(F) {
    first <- profile_format_df(profile(F, "SBS192"), factortype = TRUE)
    second <- first
    second$Value <- rev(second$Value) * 7
    oracle <- merge(
      data.frame(MutationType = first$MutationType, Strand = first$Strand, first = first$Value / sum(first$Value)),
      data.frame(MutationType = second$MutationType, Strand = second$Strand, second = second$Value / sum(second$Value)),
      by = c("MutationType", "Strand"))
    first <- first[, c("Value", "Type", "MutationType", "SubType", "Strand")]
    second <- second[sample(nrow(second)), c("Strand", "SubType", "Value", "MutationType", "Type")]
    names(second)[names(second) == "Value"] <- "Counts"
    plot <- plot_compare_profiles_diff_strand(first, second,
      profile_names = c("Original", "Expected"), condensed = TRUE)
    data <- plot$data
    difference <- data[data$name == "Difference", ]
    keys <- function(x) paste(x$MutationType, x$Strand, sep = ":")
    expected <- oracle$first - oracle$second
    equal(difference$value, expected[match(keys(difference), keys(oracle))])
    equal(sum(data$value[data$name == "Original"]), 1)
    equal(sum(data$value[data$name == "Expected"]), 1)
    equal(sort(unique(as.character(data$Strand))), c("T", "U"))
    expect(nrow(data) == 576L)
    list(data = data, rendered = render(plot, F))
  }, "reordered-channels")
  sig("plot_compare_profiles_diff_strand", function(F) {
    x <- profile_format_df(profile(F, "SBS192"), factortype = TRUE)
    expect_error(plot_compare_profiles_diff_strand(x, x[-1, ]), "same mutation types and strands")
    expect_error(plot_compare_profiles_diff_strand(x, rbind(x, x[1, ])), "uniquely identify")
    x$Value <- 0
    expect_error(plot_compare_profiles_diff_strand(x, x), "positive total")
  }, "invalid-channels", scope = "validation")
  sig("plot_compare_profiles_diff2", function(F) {
    captured <- capture_plot_data(function() plot_compare_profiles_diff2(profile(F, "SBS192"), profile_names = "Only"))
    expect(nrow(captured$data) == 192L)
    equal(sum(captured$data$value), 1)
    list(data = captured$data, rendered = render(captured$plot, F))
  }, "single-profile")
  sig("plot_compare_profiles_diff2", function(F) {
    x <- profile(F); x$Value <- 0
    expect_error(plot_compare_profiles_diff2(x), "positive total")
  }, "zero-total", scope = "validation")
  sig("position_stack_and_nudge", function(F) {
    p <- ggplot2::ggplot(data.frame(x = "A", y = c(2, 3)), ggplot2::aes(x, y)) +
      ggplot2::geom_col(position = position_stack_and_nudge(x = 0.2, vjust = 0.5))
    render(p, F)
  })
  plot_case("plot_individual_samples", function(F) list(
    exposure_refdata_input = F$exposure[F$exposure$Sample == F$samples[1], c("Signature_name", "Exposure")],
    signature_refsets_input = F$ref[, c("Signature_name", "MutationType", "Contribution")],
    seqmatrix_refdata_input = F$matrix[, 1:2]))
  plot_case("plot_cosine_heatmap_df", function(F) list(cos_sim_df(F$signature, F$signature)))
  plot_case("profile_heatmap_plot", function(F) list(data = data.frame(Sample = F$samples,
    Profile = "SBS96", Mutations = colSums(F$G))))
  plot_case("L1_normal_displot", function(F) list(similarities_dataframe = similarities(F)))
  plot_case("decompsite_distribution", function(F) list(decompsite = dplyr::mutate(
    dplyr::rename(similarities(F), Sample = Sample_Names), Cancer_Type = "CancerA")))
  plot_case("TMBplot", function(F) list(data = data.frame(Cancer_Type = F$exposure$Cancer_Type[1:24],
    Sample = F$samples, Burden = log10(rowSums(F$H) / 3100))))
  sig("prevalence_plot", function(F) {
    x <- prevalence_plot(wide_exposure(F), nmutation = 100)
    equal(x$freq_data$Value[x$freq_data$Type == "Prevalence by mutations"], colSums(F$H) / sum(F$H))
    equal(x$freq_data$Value[x$freq_data$Type == "Prevalence by samples"], rep(1, 3))
    list(data = x$freq_data, pie = render(x$pie_plot, F, "prevalence-pie"),
      bar = render(x$bar_plot, F, "prevalence-bar"))
  })
  for (name in c("Exposure_Clustering", "Exposure_Clustering_old")) local({
    target <- name
    sig(target, function(F) {
      path <- file.path(F$artifacts, "clustering.svg")
      invoke(target, sigdata = wide_exposure(F), clustern = 3, output_plot = path,
        cosinedata = data.frame(Samples = F$samples, Similarity = seq(0.8, 0.99, length.out = 24)))
      expect(file.exists(path) && file.info(path)$size > 100); TRUE
    })
  })
  sig("SBS96_Clustering", function(F) {
    path <- file.path(F$artifacts, "legacy-clustering.pdf")
    SBS96_Clustering(wide_exposure(F), cosinedata = data.frame(Samples = F$samples, Similarity = 0.99),
      clustern = 3, sigcolor = c(SBS1 = "red", SBS2 = "blue", SBS3 = "green"),
      studydata = data.frame(Samples = F$samples, Study = rep(c("A", "B"), 12)),
      studycolor = c(A = "#336699", B = "#993366"), filename = path)
    expect(file.exists(path) && file.info(path)$size > 100); TRUE
  })
  plot_case("signature_piechart", function(F) list(data = data.frame(Profile = "SBS96",
    Signature_set_name = F$set, N = 3, N2 = "3"), colset = setNames("#336699", F$set)))
  plot_case("piechart_plot", function(F) list(data = data.frame(Type = "Fixture",
    Category = c("A", "B"), Frequency = c(3, 7), Label = c("3", "7"))))
  for (name in c("barchart_plot", "barchart_plot2")) local({
    target <- name
    plot_case(target, function(F) {
      args <- list(data = data.frame(Type = "Fixture", Category = c("A", "B"), Frequency = c(0.3, 0.7)))
      if (target == "barchart_plot") args$legend_name <- "Category"
      args
    })
  })
  sig("context_plot", function(F) {
    x <- context_plot(motif(F), "TCW>TTW", data_return = TRUE)
    expect(nrow(x) == 24); expect(all(x$N1 >= 0 & x$N1 <= 1)); x
  })
  plot_case("context_plot", function(F) list(data = motif(F), pattern = "TCW>TTW"), "render")
  sig("content_extraction", function(F) {
    x <- content_extraction(motif(F)[motif(F)$Sample %in% F$samples[1:2], ])
    expect(nrow(x) > 1000); expect(all(x$N1 >= 0 & x$N1 <= 1, na.rm = TRUE)); x
  })
  plot_case("signature_association", function(F) list(data = data.frame(
    Exposure1 = F$association$x * 100, Exposure2 = F$association$y * 100)))
  sig("katPoint", function(F) {
    positions <- c(seq(10000, by = 10, length.out = 10), seq(30000, by = 10000, length.out = 10))
    x <- data.frame(build = "hg19", chr = "1", pos = positions,
      dis = c(10000, diff(positions)), ref = "C")
    value <- katPoint(x, min.mut = 5)
    expect(is.data.frame(value)); expect(all(c("start", "end", "number.mut") %in% names(value))); value
  })
  sig("mutSNP.input", function(F) {
    library(BSgenome.Hsapiens.UCSC.hg19)
    x <- mutSNP.input(data.frame(chr = "chr1", pos = seq(100000, by = 2000, length.out = 12),
      ref = "C", alt = "T"), build = "hg19")
    expect(nrow(x) > 0); expect(all(c("seq", "pos", "dis") %in% names(x))); x
  }, scope = "genome-backed")
  sig("kataegis_rainfall_plot", function(F) {
    hg19 <- data.frame(chr = "chr1", len = 249250621, start2 = 1, end2 = 249250621)
    save(hg19, file = file.path(F$work, "hg19_ref.RData"))
    path <- file.path(F$artifacts, "rainfall.svg")
    positions <- c(seq(100000, by = 20, length.out = 10), seq(150000, by = 3000, length.out = 10))
    value <- kataegis_rainfall_plot(data.frame(chr = "1", pos = positions,
      ref = "C", alt = "T"), sample_name = "Fixture", genome_build = "hg19",
      reference_data_folder = F$work, filename = path)
    expect(file.exists(path) && file.info(path)$size > 100)
    expect(is.data.frame(value) && nrow(value) > 0 && "geneID" %in% names(value),
      "Kataegis annotation returned no data; inspect run.log for a swallowed annotation error")
    value
  }, scope = "genome-backed")
  sig("seqmatrix_public_download", function(F) {
    path <- file.path(F$work, "download")
    seqmatrix_public_download(F$seqdata, path)
    archive <- paste0(path, ".tar.gz")
    expect(file.exists(archive))
    extracted <- file.path(F$work, "extracted")
    dir.create(extracted)
    utils::untar(archive, exdir = extracted)
    tables <- list.files(extracted, pattern = "[.]txt$", recursive = TRUE, full.names = TRUE)
    expect(length(tables) == 1L, "Expected one exported SBS96 table")
    value <- read.delim(tables[1], check.names = FALSE)
    equal(value, F$matrix); value
  })
  for (variant in c("numeric", "two-groups", "three-groups", "categorical", "regression")) local({
    type <- variant
    sig("mSigPortal_associaiton", function(F) {
      if (type == "regression") p <- mSigPortal_associaiton(F$association, regression = TRUE,
        formula = "y ~ x + category", type = "lm")
      else p <- mSigPortal_associaiton(F$association,
        Var1 = switch(type, numeric = "x", `two-groups` = "category", `three-groups` = "category3", categorical = "category"),
        Var2 = if (type == "categorical") "category3" else "y", type = "parametric")
      render(p, F)
    }, type)
    sig("mSigPortal_associaiton_group", function(F) {
      if (type == "regression") library(broom)
      if (type == "regression") x <- mSigPortal_associaiton_group(F$association,
        Group_Var = "group", regression = TRUE, formula = "y ~ x + category", type = "lm")
      else x <- mSigPortal_associaiton_group(F$association,
        Var1 = switch(type, numeric = "x", `two-groups` = "category", `three-groups` = "category3", categorical = "category"),
        Var2 = if (type == "categorical") "category3" else "y", Group_Var = "group", type = "parametric")
      expect(nrow(x) > 0)
      if (type == "regression") {
        expect(all(c("term", "estimate", "p.value") %in% names(x)))
        expect(nrow(x) == 4 && !any(x$term == "(Intercept)"))
      } else expect("fdr" %in% names(x))
      x
    }, type)
  })

  # Utilities exercise real RData loading semantics with the local S3 transport.
  add("utils.R", "known_error", function(F) expect_error(known_error("fixture failure"), "fixture failure", "known_error"))
  add("utils.R", "assert_signature_set", function(F) {
    assert_signature_set(F$ref, F$set)
    expect_error(assert_signature_set(F$ref[FALSE, ], "Missing", "SBS96", "PCAWG"), "Missing.*SBS96.*PCAWG", "known_error")
  })
  add("utils.R", "s3load_as", function(F) equal(s3load_as("fixture/signature_refsets.RData", F$config$bucket), F$ref))
  add("utils.R", "read_signature_refsets", function(F) equal(read_signature_refsets("fixture/signature_refsets.RData", F$config$bucket), F$ref))
  add("utils.R", "load_signature_refsets", function(F) equal(load_signature_refsets(F$config), F$ref))

  # Dispatcher serialization and error handling, in each independent wrapper.
  for (file in c("explorationWrapper.R", "visualizeWrapper.R", "associationWrapper.R")) local({
    script <- file
    add(script, "wrapper", function(F) {
      assign("msigportal.regressionSuccess", function(args, config) list(value = args$value), globalenv())
      assign("msigportal.regressionKnownError", function(args, config) known_error("expected fixture error"), globalenv())
      x <- wrapper("regressionSuccess", list(value = 42), F$config)
      if (is.character(x)) {
        decoded <- jsonlite::fromJSON(x); equal(decoded$output$value, 42)
        err <- jsonlite::fromJSON(wrapper("regressionKnownError", list(), F$config))
        equal(err$output$error, "expected fixture error")
      } else {
        equal(x$value, 42)
        err <- wrapper("regressionKnownError", list(), F$config)
        equal(err$error, "expected fixture error")
      }
      TRUE
    }, scope = "dispatcher-contract")
  })

  # Current JSON-driven Exploration APIs.
  add("explorationWrapper.R", "msigportal.msLandscape", function(F) {
    x <- msigportal.msLandscape(api_args(F))
    equal(sort(x$cosineData$sample), sort(F$samples))
    equal(x$cosineData$similarity, rep(1, 24)); expect(length(x$dendrogram$data) > 0)
    expect(nrow(x$exposureData) == nrow(F$exposure)); x
  })
  add("explorationWrapper.R", "msigportal.msLandscape", function(F) {
    args <- api_args(F)
    args$exposureData <- dplyr::filter(args$exposureData, signatureName == "SBS1")
    x <- msigportal.msLandscape(args)
    expect(length(x$dendrogram) == 0); expect(nrow(x$exposureData) == 24); x
  }, "single-signature")
  add("explorationWrapper.R", "msigportal.msLandscape", function(F) {
    args <- api_args(F)
    args$exposureData <- dplyr::filter(args$exposureData, sample %in% F$samples[1:2])
    expect_error(msigportal.msLandscape(args), "Too few samples", "known_error")
  }, "too-few-samples", scope = "validation")
  add("explorationWrapper.R", "msigportal.msDecomposition", function(F) {
    x <- msigportal.msDecomposition(api_args(F)); expect(nrow(x$download) == 24)
    equal(x$download$Cosine_similarity, rep(1, 24)); x
  })
  add("explorationWrapper.R", "msigportal.cosineSimilarity", function(F) {
    args <- api_args(F)
    x <- msigportal.cosineSimilarity(list(signatureData1 = args$signatureData, signatureData2 = args$signatureData), F$config)
    equal(x$data$similarity[x$data$sample == x$data$signature], rep(1, 3)); x
  })

  for (name in c("msigportal.getSignatureNames", "msigportal.getSampleNames", "msigportal.exposureDownload",
      "msigportal.referenceSignatures", "msigportal.mutationalProfiles", "msigportal.rsCosineSimilarity",
      "msigportal.mutationalSignatureComparison")) local({
    target <- name
    add("explorationWrapper.R", target, function(F) {
      x <- invoke(target, public_args(F), F$config)
      expect(length(x) > 0)
      if (target == "msigportal.getSignatureNames") equal(sort(x), c("SBS1", "SBS2", "SBS3"))
      if (target == "msigportal.getSampleNames") equal(x, F$samples[1:12])
      if (target != "msigportal.mutationalProfiles") checked_paths(x, F) else x
    })
  })
  legacy <- list(
    msigportal.tumorMutationalBurden = function(F, path) list(3100, path, F$exposure),
    msigportal.mutationalSignatureBurdenSeparated = function(F, path) list(3100, "CancerA", path, F$exposure),
    msigportal.mutationalSignatureBurdenAcrossCancer = function(F, path) list("SBS1", 3100, path, F$exposure),
    msigportal.mutationalSignatureAssociation = function(F, path) list(FALSE, "CancerA", FALSE, "SBS1", "SBS2", path, F$exposure),
    msigportal.mutationalSignatureDecomposition = function(F, path) list(path, file.path(F$config$savePath, "decomposition.txt"), F$exposure, F$ref, F$seqdata),
    msigportal.mutationalSignatureLandscape = function(F, path) list("CancerA", "", path, F$exposure, F$ref, F$seqdata),
    msigportal.mutationalSignaturePrevalence = function(F, path) list(100, "CancerA", path, F$exposure),
    msigportal.mutationalSignatureIndividual = function(F, path) list(F$samples[1], "CancerA", path, F$exposure, F$ref, F$seqdata))
  for (name in names(legacy)) local({
    target <- name; arguments <- legacy[[name]]
    add("explorationWrapper.R", target, function(F) {
      path <- file.path(F$config$savePath, "plot.svg")
      do.call(get(target, globalenv()), arguments(F, path))
      expect(file.exists(path) && file.info(path)$size > 100); TRUE
    }, helpers = TRUE)
  })
  for (name in c("msigportal.exposurePublic", "msigportal.exposureUser")) local({
    target <- name
    add("explorationWrapper.R", target, function(F) {
      args <- list(fn = "all", common = jsonlite::toJSON(list(study = "PCAWG", strategy = "WGS",
        rsSet = F$set, useCancerType = TRUE, cancerType = "CancerA", genome = "GRCh37"), auto_unbox = TRUE),
        files = jsonlite::toJSON(list(exposureFile = "/exposure.tsv", matrixFile = "/matrix.tsv", signatureFile = "/signature.tsv"), auto_unbox = TRUE),
        burden = jsonlite::toJSON(list(signatureName = "SBS1"), auto_unbox = TRUE),
        association = jsonlite::toJSON(list(useCancer = FALSE, both = FALSE, signatureName1 = "SBS1", signatureName2 = "SBS2"), auto_unbox = TRUE),
        landscape = '{"variableFile":""}', prevalence = '{"mutation":100}',
        individual = jsonlite::toJSON(list(sample = F$samples[1]), auto_unbox = TRUE))
      x <- invoke(target, args, F$config); checked_paths(x, F)
    })
  })

  # Visualization wrappers run with local uploads or synthetic public tables.
  viz <- c("msigportal.getReferenceSignatureSets", "msigportal.getSignaturesR", "msigportal.getPublicData",
    "msigportal.profilerSummary", "msigportal.profilerSummaryPublic", "msigportal.cosineSimilarityWithin",
    "msigportal.cosineSimilarityWithinPublic", "msigportal.cosineSimilarityRefSig", "msigportal.cosineSimilarityRefSigPublic",
    "msigportal.cosineSimilarityPublic", "msigportal.profileComparisonWithin", "msigportal.profileComparisonWithinPublic",
    "msigportal.profileComparisonRefSig", "msigportal.profileComparisonRefSigPublic", "msigportal.profileComparisonPublic",
    "msigportal.mpeaUser", "msigportal.mutationalPattern", "msigportal.mutationalPatternPublic",
    "msigportal.pca", "msigportal.pcaPublic", "msigportal.pcaWithPublic", "msigportal.downloadPublicData",
    "msigportal.getHierarchy", "msigportal.getTL", "msigportal.getTreeLeaf", "msigportal.getTreeLeafUser")
  for (name in viz) local({
    target <- name
    add("visualizeWrapper.R", target, function(F) {
      relative <- grepl("profileComparison|mutationalPattern$", target)
      args <- public_args(F, relative)
      args$matrixList <- jsonlite::toJSON(data.frame(Path = "matrix.tsv", Profile_Type = "SBS", Matrix_Size = 96))
      args$seqmatrixData <- api_args(F)$seqmatrixData
      args$seqmatrixData$profileMatrix <- "SBS96"
      args$seqmatrixData$filter <- "PASS"
      args$exposureData <- api_args(F)$exposureData
      if (target %in% c("msigportal.mpeaUser", "msigportal.mutationalPattern")) {
        write.table(F$matrix[, 1:5], F$matrix_file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
      if (target == "msigportal.pcaWithPublic") {
        upload <- F$matrix
        names(upload)[-1] <- paste0("User", names(upload)[-1])
        write.table(upload, F$matrix_file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
      if (target == "msigportal.getTreeLeaf") {
        args$seqmatrixData$profileMatrix <- "SBS96"
      }
      x <- invoke(target, args, F$config)
      expect(length(x) > 0, "Empty wrapper response")
      checked_paths(x, F)
    })
  })
  add("visualizeWrapper.R", "load_public_seqmatrix_data", function(F) {
    subset <- data.frame(Study = "PCAWG", Cancer_Type = "CancerA", Dataset = "WGS", file = "cancerA_seqmatrix.RData")
    x <- load_public_seqmatrix_data("PCAWG", "CancerA", "WGS", F$config, subset)
    equal(x$data, F$seqdata[F$seqdata$Cancer_Type == "CancerA", ]); x
  })
  add("visualizeWrapper.R", "load_public_seqmatrix_data", function(F) {
    x <- load_public_seqmatrix_data("PCAWG", "*ALL", NULL, F$config, data.frame())
    expect(grepl("Experimental strategy required", x$error)); TRUE
  }, "pancancer-validation", scope = "validation")
  for (name in c("msigportal.getAssocVarData", "msigportal.getExpVarData", "msigportal.loadCollapse",
      "msigportal.univariable", "msigportal.loadCollapseMulti", "msigportal.multivariable")) local({
    target <- name
    add("associationWrapper.R", target, function(F) {
      args <- association_args(F)
      if (target == "msigportal.multivariable") args$testType <- "lm"
      x <- invoke(target, args, F$config)
      expect(length(x) > 0); checked_paths(x, F)
    })
  })

  cases
}
