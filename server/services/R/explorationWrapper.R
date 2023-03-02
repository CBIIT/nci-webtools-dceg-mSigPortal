library(tidyverse)
library(cowplot)
library(ggsci)
library(ggrepel)
library(ggforce)
library(ggtext)
library(ggpubr)
library(jsonlite)
library(stringr)
library(aws.s3)

# capture console output for all functions called in wrapper
wrapper <- function(fn, args, config) {
  stdout <- vector("character")
  con <- textConnection("stdout", "wr", local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  output <- list()

  tryCatch(
    {
      output <- get(fn)(args, config)
    },
    error = function(e) {
      print(e)
      output <<- append(output, list(uncaughtError = e$message))
    },
    finally = {
      sink(con)
      sink(con)
      return(toJSON(list(stdout = stdout, output = output), pretty = TRUE, auto_unbox = TRUE))
    }
  )
}

# Util Functions for retrieving data
# get dataframe with column and filter arguments
# getReferenceSignatureData <- function(args, config) {
#   library(stringr)
#   s3load(paste0(config$s3Data, 'Signature/signature_refsets.RData'), config$bucket)

#   columns = unlist(args$columns, use.names = FALSE)
#   filters = args$filters
#   # get selected columns
#   data = signature_refsets %>% select(columns) %>% unique()
#   # apply filters
#   if (length(filters) > 0) {
#     for (i in 1:length(names(filters))) {
#       data = str_sort(data %>% filter(get(names(filters)[[i]]) == filters[[i]]), numeric = TRUE)
#     }
#   }

#   return(data)
# }

# retrieve signature names filtered by cancer type
getSignatureNames <- function(args, config) {
  library(stringr)
  s3load(paste0(config$s3Data, "Exposure/exposure_refdata.RData"), config$bucket)

  exposure_refdata_selected <- exposure_refdata %>% filter(Study == args$study, Dataset == args$strategy, Signature_set_name == args$rsSet)
  if (!is.null(args$cancerType)) exposure_refdata_selected <- exposure_refdata_selected %>% filter(Cancer_Type == args$cancerType)

  # available siganture name, Dropdown list for all the signature name
  signature_name_avail <- str_sort(exposure_refdata_selected %>% filter(Exposure > 0) %>% pull(Signature_name) %>% unique(), numeric = TRUE)

  return(signature_name_avail)
}

# retrieve sample names filtered by cancer type
getSampleNames <- function(args, config) {
  library(stringr)
  s3load(paste0(config$s3Data, "Exposure/exposure_refdata.RData"), config$bucket)

  exposure_refdata_selected <- exposure_refdata %>% filter(Study == args$study, Dataset == args$strategy, Signature_set_name == args$rsSet)

  # available siganture name, Dropdown list for all the signature name
  sampleNames <- str_sort(exposure_refdata_selected %>% filter(Cancer_Type == args$cancerType) %>% pull(Sample) %>% unique(), numeric = TRUE)

  return(sampleNames)
}

exposureDownload <- function(args, config) {
  s3load(paste0(config$s3Data, "Exposure/exposure_refdata.RData"), config$bucket)
  setwd(config$wd)

  exposure_refdata_selected <- exposure_refdata %>% filter(Study == args$study, Dataset == args$strategy, Signature_set_name == args$rsSet, Cancer_Type == args$cancerType)

  dfile_name <- paste(args$study, args$strategy, args$cancerType, str_remove_all(str_remove_all(args$rsSet, "\\("), "\\)"), "exposure_data.txt.gz", sep = "_")
  dfile_name <- str_replace_all(dfile_name, " +", "_")
  filepath <- paste0(config$savePath, "/", dfile_name)
  exposure_refdata_input <- exposure_refdata_selected
  exposure_refdata_input %>%
    select(Sample, Signature_name, Exposure) %>%
    pivot_wider(names_from = Signature_name, values_from = Exposure) %>%
    write_delim(file = filepath, delim = "\t", col_names = T)

  return(list(path = filepath, filename = dfile_name))
}

# Catalog - Signature  -------------------------------------------------------
# section 1: Current reference signatures in mSigPortal -------------------
referenceSignatures <- function(args, config) {
  s3load(paste0(config$s3Data, "Signature/signature_refsets.RData"), config$bucket)
  source("services/R/Sigvisualfunc.R")
  setwd(config$wd)

  plotPath <- paste0(config$savePath, "/reference_signatures.svg")

  nsig_data <- signature_refsets %>%
    group_by(Profile, Signature_set_name) %>%
    summarise(N = n_distinct(Signature_name)) %>%
    ungroup() %>%
    mutate(Profile = factor(Profile, levels = c("SBS96", "SBS192", "SBS1536", "ID83", "DBS78", "RS32")))

  nsig_data <- nsig_data %>%
    left_join(
      nsig_data %>% group_by(Profile) %>% summarise(Total = sum(N))
    ) %>%
    mutate(freq = N / Total) %>%
    mutate(N2 = if_else(freq > 0.02, as.character(N), ""))

  # put the follow pie-chart on website
  signature_piechart(nsig_data, sigsetcolor, output_plot = plotPath)

  return(list(plotPath = plotPath))
}

# section 2: Mutational signature profile  --------------------------------------------------------------
mutationalProfiles <- function(args, config) {
  s3load(paste0(config$s3Data, "Signature/signature_refsets.RData"), config$bucket)
  source("services/R/Sigvisualfunc.R")
  setwd(config$wd)

  path_profile <- paste0(config$s3Data, "Signature/Reference_Signature_Profiles_SVG/")
  signature_profile_files <- signature_refsets %>%
    select(Source, Profile, Signature_set_name, Dataset, Signature_name) %>%
    unique() %>%
    mutate(Path = str_replace_all(Signature_set_name, " ", "_"), Path = str_remove_all(Path, "[()]"), Path = paste0(path_profile, Path, "/", Signature_name, ".svg"))
  svgfile_selected <- signature_profile_files %>%
    filter(Source == args$signatureSource, Profile == args$profileName, Signature_set_name == args$rsSet, Dataset == args$experimentalStrategy, Signature_name == args$signatureName) %>%
    pull(Path)

  # fix filename
  splitPath <- strsplit(svgfile_selected, "/")[[1]]
  profileType <- str_extract_all(args$profileName, "[aA-zZ]+")
  matrixSize <- str_extract_all(args$profileName, "[0-9]+")
  filenamePrefix <- paste0(profileType, "_", matrixSize, "_plots_mSigPortal_")
  splitPath[length(splitPath)] <- paste0(filenamePrefix, splitPath[length(splitPath)])
  newPath <- paste0(splitPath, collapse = "/")

  return(list(plotPath = newPath))
}

# section3: Cosine similarities among mutational signatures -------------------------
rsCosineSimilarity <- function(args, config) {
  # The parameters will be “Matrix Size”, “Reference Signature Set1” and “Reference Signature Set2”.
  source("services/R/Sigvisualfunc.R")
  s3load(paste0(config$s3Data, "Signature/signature_refsets.RData"), config$bucket)
  setwd(config$wd)

  plotPath <- paste0(config$savePath, "signature_cos_sim_res.svg")
  txtPath <- paste0(config$savePath, "signature_cos_sim_res.txt")

  signature_refsets %>%
    filter(Profile == args$profileName) %>%
    pull(Signature_set_name) %>%
    unique()
  sigrefset1_data <- signature_refsets %>%
    filter(Profile == args$profileName, Signature_set_name == args$rsSet1) %>%
    select(Signature_name, MutationType, Contribution) %>%
    pivot_wider(names_from = Signature_name, values_from = Contribution)

  sigrefset2_data <- signature_refsets %>%
    filter(Profile == args$profileName, Signature_set_name == args$rsSet2) %>%
    select(Signature_name, MutationType, Contribution) %>%
    pivot_wider(names_from = Signature_name, values_from = Contribution)

  cos_sim_res <- cos_sim_df(sigrefset1_data, sigrefset2_data)

  # put this heatmap on the web
  plot_cosine_heatmap_df(cos_sim_res, cluster_rows = TRUE, plot_values = FALSE, output_plot = plotPath)
  # add a link bellow the heatmap
  cos_sim_res %>% write_delim(txtPath, delim = "\t", col_names = T)

  return(list(plotPath = plotPath, txtPath = txtPath))
}

# section4: Mutational signatures comparisons
## A comparison of two reference signatures
# There will be five parameters: “Profile Type”,  “Reference Signature Set1”, “Signature Name1”, “Reference Signature Set2”, “Signature Name2”;
mutationalSignatureComparison <- function(args, config) {
  # The parameters will be “Matrix Size”, “Reference Signature Set1” and “Reference Signature Set2”.
  s3load(paste0(config$s3Data, "Signature/signature_refsets.RData"), config$bucket)
  source("services/R/Sigvisualfunc.R")
  setwd(config$wd)

  plotPath <- paste0(config$savePath, "mutationalSignatureComparison.svg")
  txtPath <- paste0(config$savePath, "mutationalSignatureComparison.txt")

  signature_refsets %>%
    filter(Profile == args$profileName) %>%
    pull(Signature_set_name) %>%
    unique()
  profile1 <- signature_refsets %>%
    filter(Profile == args$profileName, Signature_set_name == args$rsSet1) %>%
    select(Signature_name, MutationType, Contribution) %>%
    pivot_wider(names_from = Signature_name, values_from = Contribution) %>%
    select(MutationType, one_of(args$signatureName1))

  profile2 <- signature_refsets %>%
    filter(Profile == args$profileName, Signature_set_name == args$rsSet2) %>%
    select(Signature_name, MutationType, Contribution) %>%
    pivot_wider(names_from = Signature_name, values_from = Contribution) %>%
    select(MutationType, one_of(args$signatureName2))

  plot_compare_profiles_diff(profile1, profile2, condensed = FALSE, output_plot = plotPath, output_data = txtPath)

  return(list(plotPath = plotPath, txtPath = txtPath))
}

# Exploration - Exposure  -------------------------------------------------------
tumorMutationalBurden <- function(genomesize, plotPath, exposure_refdata) {
  data_input <- exposure_refdata %>%
    group_by(Cancer_Type, Sample) %>%
    summarise(Burden = log10(sum(Exposure) / genomesize)) %>%
    ungroup()
  # put this barplot on the web
  TMBplot(data_input, output_plot = plotPath)
}

mutationalSignatureBurdenSeparated <- function(genomesize, cancerType, plotPath, exposure_refdata) {
  data_input <- exposure_refdata %>%
    filter(Cancer_Type == cancerType) %>%
    mutate(Burden = log10((Exposure) / genomesize)) %>%
    select(-Cancer_Type) %>%
    rename(Cancer_Type = Signature_name)
  # put this barplot on the web

  TMBplot(data_input, output_plot = plotPath)
}

mutationalSignatureBurdenAcrossCancer <- function(signatureName, genomesize, plotPath, exposure_refdata) {
  data_input <- exposure_refdata %>%
    filter(Signature_name == signatureName) %>%
    group_by(Cancer_Type, Sample) %>%
    summarise(Burden = log10(sum(Exposure) / genomesize)) %>%
    ungroup()
  # put this barplot on the web
  TMBplot(data_input, output_plot = plotPath, addnote = signatureName)
}

mutationalSignatureAssociation <- function(useCancer, cancerType, both, signatureName1, signatureName2, plotPath, exposure_refdata) {
  ## cancerType: toggle to select specific cancer type or combine all cancer type data (default)
  ## both: toggle to choose samples with both signature detected
  data_input <- left_join(
    exposure_refdata %>%
      filter(Signature_name == signatureName1) %>%
      rename(Exposure1 = Exposure) %>%
      select(-Signature_name),
    exposure_refdata %>%
      filter(Signature_name == signatureName2) %>%
      rename(Exposure2 = Exposure) %>%
      select(-Signature_name)
  )

  if (useCancer) {
    data_input <- data_input %>% filter(Cancer_Type == cancerType)
  }

  signature_association(data = data_input, signature_name_input1 = signatureName1, signature_name_input2 = signatureName2, signature_both = both, output_plot = plotPath)
}

mutationalSignatureDecomposition <- function(plotPath, txtPath, exposure_refdata_selected, signature_refsets_selected, seqmatrix_refdata_selected) {
  exposure_refdata_input <- exposure_refdata_selected %>%
    mutate(Sample = paste0(Cancer_Type, "@", Sample)) %>%
    select(Sample, Signature_name, Exposure) %>%
    pivot_wider(id_cols = Sample, names_from = Signature_name, values_from = Exposure, values_fill = 0)

  signature_refsets_input <- signature_refsets_selected %>%
    select(MutationType, Signature_name, Contribution) %>%
    pivot_wider(id_cols = MutationType, names_from = Signature_name, values_from = Contribution) %>%
    arrange(MutationType) # have to sort the mutationtype

  seqmatrix_refdata_input <- seqmatrix_refdata_selected %>%
    mutate(Sample = paste0(Cancer_Type, "@", Sample)) %>%
    select(MutationType, Sample, Mutations) %>%
    pivot_wider(id_cols = MutationType, names_from = Sample, values_from = Mutations) %>%
    arrange(MutationType) ## have to sort the mutation type

  decompsite_input <- drop_na(calculate_similarities(orignal_genomes = seqmatrix_refdata_input, signature = signature_refsets_input, signature_activaties = exposure_refdata_input))

  if (!is.data.frame(decompsite_input)) {
    stop("ERROR: Evaluating step failed due to missing the data")
  } else {
    decompsite_input <- decompsite_input %>% separate(col = Sample_Names, into = c("Cancer_Type", "Sample"), sep = "@")
    decompsite_distribution(decompsite = decompsite_input, output_plot = plotPath) # put the distribution plot online.
    decompsite_input %>% write_delim(txtPath, delim = "\t", col_names = T) ## put the link to download this table
  }
}

mutationalSignatureLandscape <- function(cancerType, varDataPath, plotPath, exposure_refdata, signature_refsets, seqmatrix_refdata) {
  exposure_refdata_input <- exposure_refdata %>%
    filter(Cancer_Type == cancerType) %>%
    select(Sample, Signature_name, Exposure) %>%
    pivot_wider(id_cols = Sample, names_from = Signature_name, values_from = Exposure)

  signature_refsets_input <- signature_refsets %>%
    select(MutationType, Signature_name, Contribution) %>%
    pivot_wider(id_cols = MutationType, names_from = Signature_name, values_from = Contribution) %>%
    arrange(MutationType) # have to sort the mutationtype

  seqmatrix_refdata_input <- seqmatrix_refdata %>%
    filter(Cancer_Type == cancerType) %>%
    select(MutationType, Sample, Mutations) %>%
    pivot_wider(id_cols = MutationType, names_from = Sample, values_from = Mutations) %>%
    arrange(MutationType) ## have to sort the mutationtype

  decompsite_input <- calculate_similarities(orignal_genomes = seqmatrix_refdata_input, signature = signature_refsets_input, signature_activaties = exposure_refdata_input)


  data_input <- exposure_refdata %>%
    filter(Cancer_Type == cancerType) %>%
    select(Sample, Signature_name, Exposure) %>%
    pivot_wider(id_cols = Sample, names_from = Signature_name, values_from = Exposure) %>%
    rename(Samples = Sample)

  data_input <- data_input %>% select_if(~ !is.numeric(.) || sum(.) > 0)
  data_input <- data_input %>% filter(rowAny(across(where(is.numeric), ~ .x > 0)))

  sigdata <- data_input

  if (!is.data.frame(decompsite_input)) {
    cosinedata <- sigdata %>%
      select(Samples) %>%
      mutate(Similarity = NA_real_)
  } else {
    cosinedata <- decompsite_input %>% select(Samples = Sample_Names, Similarity = Cosine_similarity)
  }

  ## two parameters to add the two bars: vardata1, vardata1_cat, vardata2, vardata2_cat
  # studydata <- data_input %>% select(Samples) %>% mutate(Study=if_else((seq_along(Samples) %% 2 ==0), "A","B"))
  # puritydata <-  data_input %>% select(Samples) %>% mutate(Purity=0)
  # puritydata$Purity <- runif(n=length(puritydata$Purity), min=1e-12, max=.9999999999)
  # highlight <-  c('SP124389','SP124273')
  # Exposure_Clustering(sigdata = sigdata,studydata = studydata,studyname = "VAR1",puritydata = puritydata,purityname = 'VAR2',cosinedata = cosinedata,clustern=5,output_plot = 'tmp.svg' )
  # Exposure_Clustering(sigdata = sigdata,cosinedata = cosinedata,clustern=5,output_plot = 'tmp.svg' )

  # parameter: Cancer Type, Vardata_input_file
  if (stringi::stri_length(varDataPath) > 0) {
    vardata_input <- read_delim(varDataPath, delim = "\t", col_names = T)

    vardata1_input <- vardata_input %>% select(1:2)
    colnames(vardata1_input) <- c("Samples", "Study")
    vardata1_cat_input <- if_else(is.character(vardata1_input$Study), TRUE, FALSE)

    if (dim(vardata_input)[2] > 2) {
      vardata2_input <- vardata_input %>% select(1, 3)
      colnames(vardata2_input) <- c("Samples", "Purity")
      vardata2_cat_input <- if_else(is.character(vardata2_input$Purity), TRUE, FALSE)
    } else {
      vardata2_input <- NULL
      vardata2_cat_input <- FALSE
    }
  } else {
    vardata_input <- NULL
    vardata1_input <- NULL
    vardata1_cat_input <- NULL
    vardata2_input <- NULL
    vardata2_cat_input <- NULL
  }

  Exposure_Clustering(sigdata = sigdata, studydata = vardata1_input, studydata_cat = vardata1_cat_input, puritydata = vardata2_input, puritydata_cat = vardata2_cat_input, cosinedata = cosinedata, clustern = 5, output_plot = plotPath)
}

mutationalSignaturePrevalence <- function(mutation, cancerType, plotPath, exposure_refdata) {
  require(janitor)
  require(scales)
  sigdata <- exposure_refdata %>%
    filter(Cancer_Type == cancerType) %>%
    select(Sample, Signature_name, Exposure) %>%
    pivot_wider(id_cols = Sample, names_from = Signature_name, values_from = Exposure) %>%
    rename(Samples = Sample)
  sigdata <- sigdata %>% select_if(~ !is.numeric(.) || sum(.) > 0)

  nsams <- sigdata %>%
    pivot_longer(cols = -Samples) %>%
    filter(value > 100) %>%
    dim()

  if (nsams[1] > 0) {
    prevalence_plot(sigdata = sigdata, nmutation = mutation, output_plot = plotPath)
  } else {
    stop(paste0("No signature in any sample with number of mutation larger than ", mutation))
  }
}

mutationalSignatureIndividual <- function(sample, cancerType, plotPath, exposure_refdata, signature_refsets, seqmatrix_refdata) {
  require(scales)
  exposure_refdata_input <- exposure_refdata %>%
    filter(Sample == sample) %>%
    select(Signature_name, Exposure)
  signature_refsets_input <- signature_refsets %>% select(Signature_name, MutationType, Contribution)
  seqmatrix_refdata_input <- seqmatrix_refdata %>%
    filter(Sample == sample) %>%
    select(MutationType, Mutations)

  plot_individual_samples(exposure_refdata_input = exposure_refdata_input, signature_refsets_input = signature_refsets_input, seqmatrix_refdata_input = seqmatrix_refdata_input, condensed = FALSE, output_plot = plotPath)
}

# exposurePublic <- function(fn, common, burden = '{}', association = '{}', landscape = '{}', prevalence = '{}', individual = '{}', id, pythonOutput, rootDir, config$savePath, config$s3Data, config$localData, config$bucket) {
exposurePublic <- function(args, config) {
  source("services/R/Sigvisualfunc.R")
  setwd(config$wd)

  totalTime <- proc.time()

  output <- list()
  errors <- list()
  tmbPath <- paste0(config$savePath, "tumorMutationalBurden.svg")
  signaturePath <- paste0(config$savePath, "tmbSeparated.svg")
  burdenPath <- paste0(config$savePath, "burdenAcrossCancer.svg")
  associationPath <- paste0(config$savePath, "mutationalSignatureAssociation.svg")
  decompositionPath <- paste0(config$savePath, "mutationalSignatureDecomposition.svg")
  decompositionData <- paste0(config$savePath, "mutationalSignatureDecomposition.txt")
  landscapePath <- paste0(config$savePath, "landscapeMutationalSignature.svg")
  prevalencePath <- paste0(config$savePath, "prevalenceMutationalSignature.svg")
  individualPath <- paste0(config$savePath, "msIndividual.svg")

  # parse arguments
  common <- fromJSON(args$common)
  burden <- if (is.null(args$burden)) list() else fromJSON(args$burden)
  association <- if (is.null(args$association)) list() else fromJSON(args$association)
  landscape <- if (is.null(args$landscape)) list() else fromJSON(args$landscape)
  prevalence <- if (is.null(args$prevalence)) list() else fromJSON(args$prevalence)
  individual <- if (is.null(args$individual)) list() else fromJSON(args$individual)

  genome <- case_when(
    common$study == "PCAWG" ~ "GRCh37",
    common$study == "TCGA" ~ "GRCh37",
    TRUE ~ "GRCh37"
  )
  genomesize <- genome2size(genome)

  study_signature_file <- paste0(config$s3Data, "Exposure/Study_Signatures/", common$study, "_", common$strategy, "_signature_refsets.RData")
  if (aws.s3::object_exists(study_signature_file, config$bucket)) {
    s3load(study_signature_file, config$bucket)
  } else {
    s3load(paste0(config$s3Data, "Signature/signature_refsets.RData"), config$bucket)
  }
  s3load(paste0(config$s3Data, "Exposure/", common$study, "_", common$strategy, "_exposure_refdata.RData"), config$bucket)
  s3load(paste0(config$s3Data, "Seqmatrix/seqmatrix_refdata_subset_files.RData"), config$bucket)

  # filter data
  exposure_refdata_selected <- exposure_refdata %>% filter(Study == common$study, Dataset == common$strategy, Signature_set_name == common$rsSet)
  signature_refsets_selected <- signature_refsets %>% filter(Signature_set_name == common$rsSet)

  if (common$useCancerType) {
    seqmatrixFile <- seqmatrix_refdata_subset_files %>%
      filter(Study == common$study, Dataset == common$strategy, Cancer_Type == common$cancerType) %>%
      pull(file)
  } else {
    seqmatrixFile <- paste0(common$study, "_", common$strategy, "_seqmatrix_refdata.RData")
  }
  file <- get_object(paste0(config$s3Data, "Seqmatrix/", seqmatrixFile), config$bucket)
  seqmatrix_refdata_selected <- get(load(rawConnection(file)))
  seqmatrix_refdata_selected <- seqmatrix_refdata_selected %>% filter(Profile == signature_refsets_selected$Profile[1], Cancer_Type == common$cancerType)

  # Tumor Overall Mutational Burden
  if ("all" %in% args$fn || "tmb" %in% args$fn) {
    fnTime <- proc.time()
    tryCatch(
      {
        print("Tumor Mutational Burden")
        tumorMutationalBurden(genomesize, tmbPath, exposure_refdata_selected)
        output[["tmbPath"]] <- tmbPath
      },
      error = function(e) {
        errors[["tmbError"]] <<- e$message
        print(e)
      }
    )
    print(paste0("Tumor Mutational Burden Runtime: ", (proc.time() - fnTime)[["elapsed"]]))
  }

  # Tumor Mutational Burden separated by signatures
  if ("all" %in% args$fn || "tmbSig" %in% args$fn) {
    fnTime <- proc.time()
    tryCatch(
      {
        print("Tumor Mutational Burden Separated by Signatures")
        mutationalSignatureBurdenSeparated(genomesize, common$cancerType, signaturePath, exposure_refdata_selected)
        output[["signaturePath"]] <- signaturePath
      },
      error = function(e) {
        errors[["signatureError"]] <<- e$message
        print(e)
      }
    )
    print(paste0("Tumor Mutational Burden Separated by Signatures Runtime: ", (proc.time() - fnTime)[["elapsed"]]))
  }

  # Mutational signature burden across cancer types
  if (length(burden)) {
    fnTime <- proc.time()
    tryCatch(
      {
        print("Mutational Signature Burden Across Cancer Types")
        mutationalSignatureBurdenAcrossCancer(burden$signatureName, genomesize, burdenPath, exposure_refdata_selected)
        output[["burdenPath"]] <- burdenPath
      },
      error = function(e) {
        errors[["burdenError"]] <<- e$message
        print(e)
      }
    )
    print(paste0("Mutational Signature Burden Across Cancer Types Runtime: ", (proc.time() - fnTime)[["elapsed"]]))
  }

  # Mutational Signature Association
  print(association)
  if (length(association)) {
    fnTime <- proc.time()
    tryCatch(
      {
        print("Mutational Signature Association")
        mutationalSignatureAssociation(common$useCancerType, common$cancerType, association$both, association$signatureName1, association$signatureName2, associationPath, exposure_refdata_selected)
        output[["associationPath"]] <- associationPath
      },
      error = function(e) {
        errors[["associationError"]] <<- e$message
        print(e)
      }
    )
    print(paste0("Mutational Signature Association Runtime: ", (proc.time() - fnTime)[["elapsed"]]))
  }

  # Evaluating the Performance of Mutational Signature Decomposition --------
  if ("all" %in% args$fn || "decomposition" %in% args$fn) {
    fnTime <- proc.time()
    tryCatch(
      {
        print("Evaluating the Performance of Mutational Signature Decomposition")
        mutationalSignatureDecomposition(decompositionPath, decompositionData, exposure_refdata_selected, signature_refsets_selected, seqmatrix_refdata_selected)
        output[["decompositionPath"]] <- decompositionPath
        output[["decompositionData"]] <- decompositionData
      },
      error = function(e) {
        errors[["decompositionError"]] <<- e$message
        print(e)
      }
    )
    print(paste0("Evaluating the Performance of Mutational Signature Decomposition Runtime: ", (proc.time() - fnTime)[["elapsed"]]))
  }

  # Landscape of Mutational Signature Activity
  if ("all" %in% args$fn | "landscape" %in% args$fn) {
    fnTime <- proc.time()
    tryCatch(
      {
        print("Landscape of Mutational Signature Activity")
        varDataPath <- ""
        if (stringi::stri_length(landscape$variableFile) > 0) {
          varDataPath <- file.path(paste0(config$wd, "/", config$id), landscape$variableFile)
        }
        mutationalSignatureLandscape(common$cancerType, varDataPath, landscapePath, exposure_refdata_selected, signature_refsets_selected, seqmatrix_refdata_selected)
        output[["landscapePath"]] <- landscapePath
      },
      error = function(e) {
        errors[["landscapeError"]] <<- e$message
        print(e)
      }
    )
    print(paste0("Landscape of Mutational Signature Activity Runtime: ", (proc.time() - fnTime)[["elapsed"]]))
  }

  # Prevalence plot
  if ("all" %in% args$fn | "prevalence" %in% args$fn) {
    fnTime <- proc.time()
    tryCatch(
      {
        print("Prevalence of Mutational Signature")
        mutationalSignaturePrevalence(prevalence$mutation, common$cancerType, prevalencePath, exposure_refdata_selected)
        output[["prevalencePath"]] <- prevalencePath
      },
      error = function(e) {
        errors[["prevalenceError"]] <<- e$message
        print(e)
      }
    )
    print(paste0("Prevalence of Mutational Signature Runtime: ", (proc.time() - fnTime)[["elapsed"]]))
  }

  # Individual plot
  if (length(individual)) {
    fnTime <- proc.time()
    tryCatch(
      {
        print("Mutational Signature in Individual Sample")
        mutationalSignatureIndividual(individual$sample, common$cancerType, individualPath, exposure_refdata_selected, signature_refsets_selected, seqmatrix_refdata_selected)
        output[["individualPath"]] <- individualPath
      },
      error = function(e) {
        errors[["individualError"]] <<- e$message
        print(e)
      }
    )
    print(paste0("Prevalence of Mutational Signature Runtime: ", (proc.time() - fnTime)[["elapsed"]]))
  }


  print(paste0("TotalRuntime: ", (proc.time() - totalTime)[["elapsed"]]))
  return(c(output, errors))
}

# exposureUser <- function(args$fn, files, common, burden = '{}', association = '{}', landscape = '{}', prevalence = '{}', individual = '{}', id, pythonOutput, rootDir, config$savePath, config$s3Data, config$localData, config$bucket) {
exposureUser <- function(args, config) {
  source("services/R/Sigvisualfunc.R")
  setwd(config$wd)

  output <- list()
  errors <- list()
  totalTime <- proc.time()
  tmbPath <- paste0(config$savePath, "tumorMutationalBurden.svg")
  signaturePath <- paste0(config$savePath, "tmbSeparated.svg")
  burdenPath <- paste0(config$savePath, "burdenAcrossCancer.svg")
  associationPath <- paste0(config$savePath, "mutationalSignatureAssociation.svg")
  decompositionPath <- paste0(config$savePath, "mutationalSignatureDecomposition.svg")
  decompositionData <- paste0(config$savePath, "mutationalSignatureDecomposition.txt")
  landscapePath <- paste0(config$savePath, "landscapeMutationalSignature.svg")
  prevalencePath <- paste0(config$savePath, "prevalenceMutationalSignature.svg")
  individualPath <- paste0(config$savePath, "msIndividual.svg")

  # parse arguments
  common <- fromJSON(args$common)
  files <- fromJSON(args$files)
  burden <- if (is.null(args$burden)) list() else fromJSON(args$burden)
  association <- if (is.null(args$association)) list() else fromJSON(args$association)
  landscape <- if (is.null(args$landscape)) list() else fromJSON(args$landscape)
  prevalence <- if (is.null(args$prevalence)) list() else fromJSON(args$prevalence)
  individual <- if (is.null(args$individual)) list() else fromJSON(args$individual)

  exposure_refdata_selected <- read_delim(file.path(paste0(config$wd, "/", config$id), files$exposureFile), delim = "\t", col_names = T)
  seqmatrix_refdata_selected <- read_delim(file.path(paste0(config$wd, "/", config$id), files$matrixFile), delim = "\t", col_names = T)

  if (stringi::stri_length(files$signatureFile) > 0) {
    # if using user uploaded signature file
    signature_refsets_selected <- read_delim(file.path(paste0(config$wd, "/", config$id), files$signatureFile), delim = "\t", col_names = T)
  } else {
    # else use public signature data
    s3load(paste0(config$s3Data, "Signature/signature_refsets.RData"), config$bucket)

    signature_refsets_selected <- signature_refsets %>%
      filter(Signature_set_name == common$rsSet) %>%
      select(MutationType, Signature_name, Contribution) %>%
      pivot_wider(names_from = Signature_name, values_from = Contribution)
  }

  genomesize <- genome2size(common$genome)
  cancer_type_user <- "Input"

  ## format the input as suggest by the public dataset ##
  colnames(exposure_refdata_selected)[1] <- "Sample"
  colnames(seqmatrix_refdata_selected)[1] <- "MutationType"
  colnames(signature_refsets_selected)[1] <- "MutationType"
  seqmatrix_refdata_selected <- seqmatrix_refdata_selected %>%
    select(1, any_of(exposure_refdata_selected$Sample)) %>%
    profile_format_df()
  signature_refsets_selected <- signature_refsets_selected %>%
    select(MutationType, any_of(colnames(exposure_refdata_selected))) %>%
    profile_format_df()

  exposure_refdata_selected <- exposure_refdata_selected %>%
    pivot_longer(cols = -Sample, names_to = "Signature_name", values_to = "Exposure") %>%
    mutate(Cancer_Type = cancer_type_user)
  signature_refsets_selected <- signature_refsets_selected %>%
    select(-Type, -SubType) %>%
    pivot_longer(cols = -MutationType, names_to = "Signature_name", values_to = "Contribution")
  seqmatrix_refdata_selected <- seqmatrix_refdata_selected %>%
    select(-Type, -SubType) %>%
    pivot_longer(cols = -MutationType, names_to = "Sample", values_to = "Mutations") %>%
    mutate(Cancer_Type = cancer_type_user)


  ## Tumor Overall Mutational Burden
  if ("all" %in% args$fn || "tmb" %in% args$fn) {
    fnTime <- proc.time()
    tryCatch(
      {
        print("Tumor Mutational Burden")
        tumorMutationalBurden(genomesize, tmbPath, exposure_refdata_selected)
        output[["tmbPath"]] <- tmbPath
      },
      error = function(e) {
        errors[["tmbError"]] <<- e$message
        print(e)
      }
    )
    print(paste0("Tumor Mutational Burden Runtime: ", (proc.time() - fnTime)[["elapsed"]]))
  }

  # Tumor Mutational Burden separated by signatures
  if ("all" %in% args$fn || "tmbSig" %in% args$fn) {
    fnTime <- proc.time()
    tryCatch(
      {
        print("Tumor Mutational Burden Separated by Signatures")
        mutationalSignatureBurdenSeparated(genomesize, cancer_type_user, signaturePath, exposure_refdata_selected)
        output[["signaturePath"]] <- signaturePath
      },
      error = function(e) {
        errors[["signatureError"]] <<- e$message
        print(e)
      }
    )
    print(paste0("Tumor Mutational Burden Separated by Signatures Runtime: ", (proc.time() - fnTime)[["elapsed"]]))
  }
  # Mutational signature burden across cancer types
  if (length(burden)) {
    fnTime <- proc.time()
    tryCatch(
      {
        print("Mutational Signature Burden Across Cancer Types")
        mutationalSignatureBurdenAcrossCancer(burden$signatureName, genomesize, burdenPath, exposure_refdata_selected)
        output[["burdenPath"]] <- burdenPath
      },
      error = function(e) {
        errors[["burdenError"]] <<- e$message
        print(e)
      }
    )
    print(paste0("Mutational Signature Burden Across Cancer Types Runtime: ", (proc.time() - fnTime)[["elapsed"]]))
  }
  # Mutational Signature Association
  if (length(association)) {
    fnTime <- proc.time()
    tryCatch(
      {
        print("Mutational Signature Association")
        mutationalSignatureAssociation(common$useCancerType, cancer_type_user, association$both, association$signatureName1, association$signatureName2, associationPath, exposure_refdata_selected)
        output[["associationPath"]] <- associationPath
      },
      error = function(e) {
        errors[["associationError"]] <<- e$message
        print(e)
      }
    )
    print(paste0("Mutational Signature Association Runtime: ", (proc.time() - fnTime)[["elapsed"]]))
  }
  # Evaluating the Performance of Mutational Signature Decomposition --------
  if ("all" %in% args$fn || "decomposition" %in% args$fn) {
    fnTime <- proc.time()
    tryCatch(
      {
        print("Evaluating the Performance of Mutational Signature Decomposition")
        mutationalSignatureDecomposition(decompositionPath, decompositionData, exposure_refdata_selected, signature_refsets_selected, seqmatrix_refdata_selected)
        output[["decompositionPath"]] <- decompositionPath
        output[["decompositionData"]] <- decompositionData
      },
      error = function(e) {
        errors[["decompositionError"]] <<- e$message
        print(e)
      }
    )
    print(paste0("Evaluating the Performance of Mutational Signature Decomposition Runtime: ", (proc.time() - fnTime)[["elapsed"]]))
  }
  # Landscape of Mutational Signature Activity
  if ("all" %in% args$fn | "landscape" %in% args$fn) {
    fnTime <- proc.time()
    tryCatch(
      {
        print("Landscape of Mutational Signature Activity")
        varDataPath <- ""
        if (stringi::stri_length(landscape$variableFile) > 0) {
          varDataPath <- file.path(paste0(config$wd, "/", config$id), landscape$variableFile$signatureFile)
        }
        mutationalSignatureLandscape(cancer_type_user, varDataPath, landscapePath, exposure_refdata_selected, signature_refsets_selected, seqmatrix_refdata_selected)
        output[["landscapePath"]] <- landscapePath
      },
      error = function(e) {
        errors[["landscapeError"]] <<- e$message
        print(e)
      }
    )
    print(paste0("Landscape of Mutational Signature Activity Runtime: ", (proc.time() - fnTime)[["elapsed"]]))
  }
  # Prevalence plot
  if ("all" %in% args$fn | "prevalence" %in% args$fn) {
    fnTime <- proc.time()
    tryCatch(
      {
        print("Prevalence of Mutational Signature")
        mutationalSignaturePrevalence(prevalence$mutation, cancer_type_user, prevalencePath, exposure_refdata_selected)
        output[["prevalencePath"]] <- prevalencePath
      },
      error = function(e) {
        errors[["prevalenceError"]] <<- e$message
        print(e)
      }
    )
    print(paste0("Prevalence of Mutational Signature Runtime: ", (proc.time() - fnTime)[["elapsed"]]))
  }

  # Individual plot
  if (length(individual)) {
    fnTime <- proc.time()
    tryCatch(
      {
        print("Mutational Signature in Individual Sample")
        mutationalSignatureIndividual(individual$sample, common$cancerType, individualPath, exposure_refdata_selected, signature_refsets_selected, seqmatrix_refdata_selected)
        output[["individualPath"]] <- individualPath
      },
      error = function(e) {
        errors[["individualError"]] <<- e$message
        print(e)
      }
    )
    print(paste0("Prevalence of Mutational Signature Runtime: ", (proc.time() - fnTime)[["elapsed"]]))
  }


  print(paste0("Total Runtime: ", (proc.time() - totalTime)[["elapsed"]]))
  return(c(output, errors))
}

msLandscape <- function(args, config) {
  source("services/R/Sigvisualfunc.R")

  exposureData <- args$exposureData %>%
    select(Sample = sample, Signature_name = signatureName, Exposure = exposure) %>%
    pivot_wider(id_cols = Sample, names_from = Signature_name, values_from = Exposure)

  signatureData <- args$signatureData %>%
    select(MutationType = mutationType, Signature_name = signatureName, Contribution = contribution) %>%
    pivot_wider(id_cols = MutationType, names_from = Signature_name, values_from = Contribution) %>%
    arrange(MutationType) # have to sort the mutationtype

  seqmatrixData <- args$seqmatrixData %>%
    filter(if ("profile" %in% names(.)) profile == args$signatureData$profile[1] else TRUE) %>%
    filter(if ("matrix" %in% names(.)) matrix == args$signatureData$matrix[1] else TRUE) %>%
    select(MutationType = mutationType, Sample = sample, Mutations = mutations) %>%
    pivot_wider(id_cols = MutationType, names_from = Sample, values_from = Mutations) %>%
    arrange(MutationType) ## have to sort the mutationtype

  decompsite_input <- calculate_similarities(orignal_genomes = seqmatrixData, signature = signatureData, signature_activaties = exposureData)

  transformExposure <- args$exposureData %>%
    select(sample, signatureName, exposure) %>%
    pivot_wider(id_cols = sample, names_from = signatureName, values_from = exposure) %>%
    select_if(~ !is.numeric(.) || sum(.) > 0) %>%
    filter(rowAny(across(where(is.numeric), ~ .x > 0)))

  if (!is.data.frame(decompsite_input)) {
    cosineData <- transformExposure %>%
      select(sample) %>%
      mutate(similarity = NA_real_)
  } else {
    cosineData <- decompsite_input %>% select(sample = Sample_Names, similarity = Cosine_similarity)
  }

  # calculate hierarchy
  transformExposure <- transformExposure %>%
    select(where(~ is.character(.x) || sum(.x) != 0)) %>%
    janitor::adorn_percentages("row")
  mdata <- as.matrix(transformExposure[, -1])
  rownames(mdata) <- transformExposure$sample
  clustern <- ifelse(dim(mdata)[1] < 10, 2L, 5)
  cluster <- factoextra::hcut(mdata, k = clustern, hc_func = "hclust", hc_metric = "euclidean", hc_method = "ward.D2", stand = TRUE)

  # create ggplot plotly dendrogram
  dendrogramPlot <- ggdendro::ggdendrogram(cluster)
  dendrogram_json <- fromJSON(plotly::plotly_json(dendrogramPlot))

  # sort according to hierarchy order
  exposureData <- args$exposureData %>%
    arrange(factor(sample, levels = cluster$labels[cluster$order])) %>%
    filter(exposure > 0)
  cosineData <- cosineData %>%
    arrange(factor(sample, levels = cluster$labels[cluster$order]))

  return(list(cosineData = cosineData, exposureData = exposureData, dendrogram = dendrogram_json))
}

msDecomposition <- function(args, config) {
  source("services/R/Sigvisualfunc.R")

  exposureData <- args$exposureData %>%
    mutate(Sample = paste0(cancer, "@", sample)) %>%
    select(Sample, Signature_name = signatureName, Exposure = exposure) %>%
    pivot_wider(id_cols = Sample, names_from = Signature_name, values_from = Exposure)

  signatureData <- args$signatureData %>%
    select(MutationType = mutationType, Signature_name = signatureName, Contribution = contribution) %>%
    pivot_wider(id_cols = MutationType, names_from = Signature_name, values_from = Contribution) %>%
    arrange(MutationType) # have to sort the mutationtype

  seqmatrixData <- args$seqmatrixData %>%
    filter(if ("profile" %in% names(.)) profile == args$signatureData$profile[1] else TRUE) %>%
    filter(if ("matrix" %in% names(.)) matrix == args$signatureData$matrix[1] else TRUE) %>%
    mutate(Sample = paste0(cancer, "@", sample)) %>%
    select(MutationType = mutationType, Sample, Mutations = mutations) %>%
    pivot_wider(id_cols = MutationType, names_from = Sample, values_from = Mutations) %>%
    arrange(MutationType) ## have to sort the mutationtype

  decompsite_input <- calculate_similarities(orignal_genomes = seqmatrixData, signature = signatureData, signature_activaties = exposureData)
  #  before <- calculate_similarities(orignal_genomes = seqmatrixData, signature = signatureData, signature_activaties = exposureData)

  # decompsite_input <- decompsite_input %>% separate(col = Sample_Names, into = c("cancer", "sample"), sep = "@")

  # return(list(data = decompsite_input))
  decompsite <- decompsite_input %>% separate(col = Sample_Names, into = c("Cancer_Type", "Sample"), sep = "@")

  fealist <- c("Cancer_Type", "Sample", "Cosine_similarity", "100-L1_Norm_%", "100-L2_Norm_%", "KL_Divergence", "Correlation")
  decompsite2 <- decompsite %>%
    select(one_of(fealist)) %>%
    pivot_longer(cols = -c(Cancer_Type, Sample))
  mtmp <- decompsite %>%
    group_by(Cancer_Type) %>%
    summarise(m = median(Cosine_similarity, na.rm = TRUE))
  mtmp2 <- decompsite2 %>%
    group_by(Cancer_Type, name) %>%
    summarise(m = median(value, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(m = if_else(name == "KL_Divergence", 1 - m, m)) %>%
    group_by(name) %>%
    arrange(desc(m)) %>%
    mutate(Seq = seq_along(Cancer_Type)) %>%
    ungroup()
  ntmp <- decompsite %>%
    count(Cancer_Type) %>%
    mutate(Caner_type2 = paste0(Cancer_Type, " (", n, ")")) %>%
    left_join(mtmp) %>%
    arrange(m)

  data <- decompsite2 %>%
    left_join(mtmp2) %>%
    left_join(ntmp) %>%
    select(cancer = Cancer_Type, sample = Sample, name, value, m, n, Seq)

  return(list(data = data, download = decompsite))
}

cosineSimilarity <- function(args, config) {
  source("services/R/Sigvisualfunc.R")

  signatureData1 <- args$signatureData1 %>%
    select(Signature_name = signatureName, MutationType = mutationType, Contribution = contribution) %>%
    pivot_wider(names_from = Signature_name, values_from = Contribution)

  signatureData2 <- args$signatureData2 %>%
    select(Signature_name = signatureName, MutationType = mutationType, Contribution = contribution) %>%
    pivot_wider(names_from = Signature_name, values_from = Contribution)

  cosSim <- cos_sim_df(signatureData1, signatureData2)

  colnames(cosSim)[1] <- "sample"
  cosSimMatrix <- as.matrix(cosSim[, -1])
  rownames(cosSimMatrix) <- cosSim[[1]]

  hc.sample <- hclust(dist(t(cosSimMatrix)), method = "complete")
  signatureOrder <- rownames(t(cosSimMatrix))[hc.sample$order]
  hc.sample <- hclust(dist(cosSimMatrix), method = "complete")
  sampleOrder <- rownames(cosSimMatrix)[hc.sample$order]

  data <- cosSim %>%
    pivot_longer(-1, names_to = "signature", values_to = "similarity") %>%
    select(sample, signature, similarity)
  data$signature <- factor(data$signature, levels = signatureOrder)
  data$sample <- factor(data$sample, levels = sampleOrder)

  return(list(original = cosSim, data = data, sampleOrder = sampleOrder, signatureOrder = signatureOrder))
}
