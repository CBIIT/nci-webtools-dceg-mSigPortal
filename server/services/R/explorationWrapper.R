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




# Util Functions for retrieving data
# get dataframe with column and filter arguments
getReferenceSignatureData <- function(args, s3Data, localData, bucket) {
  s3load(paste0(s3Data, 'Signature/signature_refsets.RData'), bucket)

  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()

    columns = unlist(args$columns, use.names = FALSE)
    filters = args$filters
    # get selected columns
    data = signature_refsets %>% select(columns) %>% unique()
    # apply filters
    if (length(filters) > 0) {
      for (i in 1:length(names(filters))) {
        data = data %>% filter(get(names(filters)[[i]]) == filters[[i]])
      }
    }

    output = list('data' = data)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

# retrieve signature names filtered by cancer type
getSignatureNames <- function(args, s3Data, localData, bucket) {
  s3load(paste0(s3Data, 'Exposure/exposure_refdata.RData'), bucket)

  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    exposure_refdata_selected <- exposure_refdata %>% filter(Study == args$study, Dataset == args$strategy, Signature_set_name == args$refSignatureSet)

    # available siganture name, Dropdown list for all the signature name
    signature_name_avail <- exposure_refdata_selected %>% filter(Cancer_Type == args$cancerType, Exposure > 0) %>% pull(Signature_name) %>% unique()

    output = list(data = signature_name_avail)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

# retrieve sample names filtered by cancer type
getSampleNames <- function(args, s3Data, localData, bucket) {
  s3load(paste0(s3Data, 'Exposure/exposure_refdata.RData'), bucket)

  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    exposure_refdata_selected <- exposure_refdata %>% filter(Study == args$study, Dataset == args$strategy, Signature_set_name == args$refSignatureSet)

    # available siganture name, Dropdown list for all the signature name
    sampleNames <- exposure_refdata_selected %>% filter(Cancer_Type == args$cancerType) %>% pull(Sample) %>% unique()

    output = list(data = sampleNames)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), auto_unbox = TRUE))
  })
}

exposureDownload <- function(study, strategy, refSignatureSet, cancerType, projectID, pythonOutput, rootDir, savePath, s3Data, localData, bucket) {
  s3load(paste0(s3Data, 'Exposure/exposure_refdata.RData'), bucket)

  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    exposure_refdata_selected <- exposure_refdata %>% filter(Study == study, Dataset == strategy, Signature_set_name == refSignatureSet, Cancer_Type == cancerType)

    dfile_name <- paste(study, strategy, cancerType, str_remove_all(str_remove_all(refSignatureSet, '\\('), '\\)'), 'exposure_data.txt.gz', sep = '_')
    dfile_name <- str_replace_all(dfile_name, ' +', '_')
    filepath = paste0(savePath, '/', dfile_name)
    exposure_refdata_input <- exposure_refdata_selected
    exposure_refdata_input %>%
      select(Sample, Signature_name, Exposure) %>%
      pivot_wider(names_from = Signature_name, values_from = Exposure) %>%
      write_delim(file = filepath, delim = '\t', col_names = T)


    output = list(path = filepath, filename = dfile_name)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), auto_unbox = TRUE))
  })
}

# Signature Explore -------------------------------------------------------
# section 1: Current reference signatures in mSigPortal -------------------
referenceSignatures <- function(projectID, pythonOutput, rootDir, savePath, s3Data, localData, bucket) {
  source('services/R/Sigvisualfunc.R')

  s3load(paste0(s3Data, 'Signature/signature_refsets.RData'), bucket)

  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    plotPath = paste0(savePath, '/reference_signatures.svg')

    nsig_data <- signature_refsets %>%
      group_by(Profile, Signature_set_name) %>%
      summarise(N = n_distinct(Signature_name)) %>%
      ungroup() %>%
      mutate(Profile = factor(Profile, levels = c("SBS96", "SBS192", "SBS1536", "ID83", "DBS78", "RS32")))

    nsig_data <- nsig_data %>% left_join(
      nsig_data %>% group_by(Profile) %>% summarise(Total = sum(N))
      ) %>%
      mutate(freq = N / Total) %>%
      mutate(N2 = if_else(freq > 0.02, as.character(N), ""))

    # put the follow pie-chart on website
    signature_piechart(nsig_data, sigsetcolor, output_plot = plotPath)

    output = list('plotPath' = plotPath)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

# section 2: Mutational signature profile  --------------------------------------------------------------
mutationalProfiles <- function(signatureSource, profileName, refSignatureSet, experimentalStrategy, signatureName, projectID, pythonOutput, rootDir, savePath, s3Data, localData, bucket) {
  source('services/R/Sigvisualfunc.R')

  s3load(paste0(s3Data, 'Signature/signature_refsets.RData'), bucket)

  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()

    path_profile <- paste0(s3Data, 'Signature/Reference_Signature_Profiles_SVG/')
    signature_profile_files <- signature_refsets %>% select(Source, Profile, Signature_set_name, Dataset, Signature_name) %>% unique() %>% mutate(Path = str_replace_all(Signature_set_name, " ", "_"), Path = str_remove_all(Path, "[()]"), Path = paste0(path_profile, Path, "/", Signature_name, ".svg"))
    svgfile_selected <- signature_profile_files %>%
      filter(Source == signatureSource, Profile == profileName, Signature_set_name == refSignatureSet, Dataset == experimentalStrategy, Signature_name == signatureName) %>% pull(Path)

    # fix filename
    splitPath = strsplit(svgfile_selected, '/')[[1]]
    profileType = str_extract_all(profileName, "[aA-zZ]+")
    matrixSize = str_extract_all(profileName, "[0-9]+")
    filenamePrefix = paste0(profileType, '_', matrixSize, '_plots_mSigPortal_')
    splitPath[length(splitPath)] = paste0(filenamePrefix, splitPath[length(splitPath)])
    newPath = paste0(splitPath, collapse = "/")

    output = list('plotPath' = newPath)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

# section3: Cosine similarities among mutational signatures -------------------------
cosineSimilarity <- function(profileName, refSignatureSet1, refSignatureSet2, projectID, pythonOutput, rootDir, savePath, s3Data, localData, bucket) {
  # The parameters will be “Matrix Size”, “Reference Signature Set1” and “Reference Signature Set2”. 
  source('services/R/Sigvisualfunc.R')

  s3load(paste0(s3Data, 'Signature/signature_refsets.RData'), bucket)

  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    plotPath = paste0(savePath, 'signature_cos_sim_res.svg')
    txtPath = paste0(savePath, 'signature_cos_sim_res.txt')

    signature_refsets %>% filter(Profile == profileName) %>% pull(Signature_set_name) %>% unique()
    sigrefset1_data <- signature_refsets %>%
      filter(Profile == profileName, Signature_set_name == refSignatureSet1) %>%
      select(Signature_name, MutationType, Contribution) %>%
      pivot_wider(names_from = Signature_name, values_from = Contribution)

    sigrefset2_data <- signature_refsets %>%
      filter(Profile == profileName, Signature_set_name == refSignatureSet2) %>%
      select(Signature_name, MutationType, Contribution) %>%
      pivot_wider(names_from = Signature_name, values_from = Contribution)

    cos_sim_res = cos_sim_df(sigrefset1_data, sigrefset2_data)

    # put this heatmap on the web
    plot_cosine_heatmap_df(cos_sim_res, cluster_rows = TRUE, plot_values = FALSE, output_plot = plotPath)
    # add a link bellow the heatmap
    cos_sim_res %>% write_delim(txtPath, delim = '\t', col_names = T)

    output = list('plotPath' = plotPath, 'txtPath' = txtPath)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

# section4: Mutational signatures comparisons
## A comparison of two reference signatures
# There will be five parameters: “Profile Type”,  “Reference Signature Set1”, “Signature Name1”, “Reference Signature Set2”, “Signature Name2”;
mutationalSignatureComparison <- function(profileName, refSignatureSet1, signatureName1, refSignatureSet2, signatureName2, projectID, pythonOutput, rootDir, savePath, s3Data, localData, bucket) {
  # The parameters will be “Matrix Size”, “Reference Signature Set1” and “Reference Signature Set2”. 
  source('services/R/Sigvisualfunc.R')

  s3load(paste0(s3Data, 'Signature/signature_refsets.RData'), bucket)

  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    plotPath = paste0(savePath, 'mutationalSignatureComparison.svg')

    profile_name_input <- "SBS96" # profile type
    signatureset_name1 <- "COSMIC v3 Signatures (SBS)"
    signatureset_name2 <- "COSMIC v3 Signatures (SBS)"
    signature_name1 <- "SBS1"
    signature_name2 <- "SBS5"

    signature_refsets %>% filter(Profile == profileName) %>% pull(Signature_set_name) %>% unique()
    profile1 <- signature_refsets %>%
      filter(Profile == profileName, Signature_set_name == refSignatureSet1) %>%
      select(Signature_name, MutationType, Contribution) %>%
      pivot_wider(names_from = Signature_name, values_from = Contribution) %>%
      select(MutationType, one_of(signatureName1))

    profile2 <- signature_refsets %>%
      filter(Profile == profileName, Signature_set_name == refSignatureSet2) %>%
      select(Signature_name, MutationType, Contribution) %>%
      pivot_wider(names_from = Signature_name, values_from = Contribution) %>%
      select(MutationType, one_of(signatureName2))

    plot_compare_profiles_diff(profile1, profile2, condensed = FALSE, output_plot = plotPath)

    output = list('plotPath' = plotPath)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

# Exposure Explore -------------------------------------------------------
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

  error = signature_association(data = data_input, signature_name_input1 = signatureName1, signature_name_input2 = signatureName2, signature_both = both, output_plot = plotPath)
  if (!is.null(error)) stop(error)
}

mutationalSignatureDecomposition <- function(plotPath, s3Data, localData, exposure_refdata_selected, signature_refsets_selected, seqmatrix_refdata_selected) {
  exposure_refdata_input <- exposure_refdata_selected %>% mutate(Sample = paste0(Cancer_Type, "@", Sample)) %>%
    select(Sample, Signature_name, Exposure) %>%
    pivot_wider(id_cols = Sample, names_from = Signature_name, values_from = Exposure, values_fill = 0)

  signature_refsets_input <- signature_refsets_selected %>%
    select(MutationType, Signature_name, Contribution) %>%
    pivot_wider(id_cols = MutationType, names_from = Signature_name, values_from = Contribution) %>%
    arrange(MutationType) # have to sort the mutationtype

  seqmatrix_refdata_input <- seqmatrix_refdata_selected %>% mutate(Sample = paste0(Cancer_Type, "@", Sample)) %>%
    select(MutationType, Sample, Mutations) %>%
    pivot_wider(id_cols = MutationType, names_from = Sample, values_from = Mutations) %>%
    arrange(MutationType) ## have to sort the mutation type

  decompsite_input <- drop_na(calculate_similarities(orignal_genomes = seqmatrix_refdata_input, signature = signature_refsets_input, signature_activaties = exposure_refdata_input))

  if (!is.data.frame(decompsite_input)) {
    stop('Evaluating step failed due to missing the data')
  } else {
    decompsite_input <- decompsite_input %>% separate(col = Sample_Names, into = c('Cancer_Type', 'Sample'), sep = '@')
    decompsite_distribution(decompsite = decompsite_input, output_plot = plotPath) # put the distribution plot online.
    decompsite_input %>% write_delim(s3Data, delim = '\t', col_names = T) ## put the link to download this table
  }
}

mutationalSignatureLandscape <- function(cancerType, vars3Data, localData, plotPath, exposure_refdata, signature_refsets, seqmatrix_refdata) {
  exposure_refdata_input <- exposure_refdata %>% filter(Cancer_Type == cancerType) %>%
      select(Sample, Signature_name, Exposure) %>%
      pivot_wider(id_cols = Sample, names_from = Signature_name, values_from = Exposure)

  signature_refsets_input <- signature_refsets %>%
      select(MutationType, Signature_name, Contribution) %>%
      pivot_wider(id_cols = MutationType, names_from = Signature_name, values_from = Contribution) %>%
      arrange(MutationType) # have to sort the mutationtype

  seqmatrix_refdata_input <- seqmatrix_refdata %>% filter(Cancer_Type == cancerType) %>%
      select(MutationType, Sample, Mutations) %>%
      pivot_wider(id_cols = MutationType, names_from = Sample, values_from = Mutations) %>%
      arrange(MutationType) ## have to sort the mutationtype

  decompsite_input <- calculate_similarities(orignal_genomes = seqmatrix_refdata_input, signature = signature_refsets_input, signature_activaties = exposure_refdata_input)


  data_input <- exposure_refdata %>%
      filter(Cancer_Type == cancerType) %>%
      select(Sample, Signature_name, Exposure) %>%
      pivot_wider(id_cols = Sample, names_from = Signature_name, values_from = Exposure) %>%
      rename(Samples = Sample)

  data_input <- data_input %>% select_if(~!is.numeric(.) || sum(.) > 0)
  data_input <- data_input %>% filter(rowAny(across(where(is.numeric), ~ .x > 0)))

  sigdata <- data_input

  if (!is.data.frame(decompsite_input)) {
    cosinedata <- sigdata %>% select(Samples) %>% mutate(Similarity = NA_real_)
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
    vardata_input <- read_delim(vars3Data, localData, delim = '\t', col_names = T)

    vardata1_input <- vardata_input %>% select(1:2)
    colnames(vardata1_input) <- c('Samples', 'Study')
    vardata1_cat_input <- if_else(is.character(vardata1_input$Study), TRUE, FALSE)

    if (dim(vardata_input)[2] > 2) {
      vardata2_input <- vardata_input %>% select(1, 3)
      colnames(vardata2_input) <- c('Samples', 'Purity')
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
  sigdata <- sigdata %>% select_if(~!is.numeric(.) || sum(.) > 0)

  nsams <- sigdata %>% pivot_longer(cols = -Samples) %>% filter(value > 100) %>% dim()

  if (nsams[1] > 0) {
    prevalence_plot(sigdata = sigdata, nmutation = mutation, output_plot = plotPath)
  } else {
    stop(paste0('No signature in any sample with number of mutation larger than ', mutation))
  }
}

mutationalSignatureIndividual <- function(sample, cancerType, plotPath, exposure_refdata, signature_refsets, seqmatrix_refdata) {
  require(scales)
  exposure_refdata_input <- exposure_refdata %>% filter(Sample == sample) %>% select(Signature_name, Exposure)
  signature_refsets_input <- signature_refsets %>% select(Signature_name, MutationType, Contribution)
  seqmatrix_refdata_input <- seqmatrix_refdata %>% filter(Sample == sample) %>% select(MutationType, Mutations)

  plot_individual_samples(exposure_refdata_input = exposure_refdata_input, signature_refsets_input = signature_refsets_input, seqmatrix_refdata_input = seqmatrix_refdata_input, condensed = FALSE, output_plot = plotPath)
}

exposurePublic <- function(fn, common, burden = '{}', association = '{}', landscape = '{}', prevalence = '{}', individual = '{}', projectID, pythonOutput, rootDir, savePath, s3Data, localData, bucket) {
  tryCatch({
    source('services/R/Sigvisualfunc.R')
    con <- textConnection('stdout', 'wr', local = TRUE)
    sink(con, type = "message")
    sink(con, type = "output")

    totalTime = proc.time()

    s3load(paste0(s3Data, 'Signature/signature_refsets.RData'), bucket)
    s3load(paste0(s3Data, 'Seqmatrix/seqmatrix_refdata_subset_files.RData'), bucket)
    s3load(paste0(s3Data, 'Exposure/exposure_refdata.RData'), bucket)

    output = list()
    errors = list()
    tmbPath = paste0(savePath, 'tumorMutationalBurden.svg')
    signaturePath = paste0(savePath, 'tmbSeparated.svg')
    burdenPath = paste0(savePath, 'burdenAcrossCancer.svg')
    associationPath = paste0(savePath, 'mutationalSignatureAssociation.svg')
    decompositionPath = paste0(savePath, 'mutationalSignatureDecomposition.svg')
    decompositionData = paste0(savePath, 'mutationalSignatureDecomposition.txt')
    landscapePath = paste0(savePath, 'landscapeMutationalSignature.svg')
    prevalencePath = paste0(savePath, 'prevalenceMutationalSignature.svg')
    individualPath = paste0(savePath, 'msIndividual.svg')

    # parse arguments
    common = fromJSON(common)
    burden = fromJSON(burden)
    association = fromJSON(association)
    landscape = fromJSON(landscape)
    prevalence = fromJSON(prevalence)
    individual = fromJSON(individual)
    exposure_refdata_selected <- exposure_refdata %>% filter(Study == common$study, Dataset == common$strategy, Signature_set_name == common$refSignatureSet)

    genome <- case_when(
    common$study == "PCAWG" ~ "GRCh37",
    common$study == "TCGA" ~ "GRCh37",
    TRUE ~ "GRCh37"
  )
    genomesize = genome2size(genome)

    signature_refsets_selected <- signature_refsets %>%
      filter(Signature_set_name == common$refSignatureSet)

    # seqmatrix_refdata_selected <- seqmatrix_refdata %>% filter(Study == common$study, Dataset == common$strategy, Profile == signature_refsets_selected$Profile[1])
    if (common$useCancerType) {
      seqmatrixFile <- seqmatrix_refdata_subset_files %>% filter(Study == common$study, Dataset == common$strategy, Cancer_Type == common$cancerType) %>% pull(file)
    } else {
      seqmatrixFile <- paste0(common$study, '_', common$strategy, '_seqmatrix_refdata.RData')
    }

    seqmatrix_refdata_selected = NULL

    file <- get_object(paste0(s3Data, 'Seqmatrix/', seqmatrixFile), bucket)
    seqmatrix_refdata_selected <- get(load(rawConnection(file)))

    seqmatrix_refdata_selected = seqmatrix_refdata_selected %>% filter(Profile == signature_refsets_selected$Profile[1])
    # Tumor Overall Mutational Burden
    if ('all' %in% fn || 'tmb' %in% fn) {
      fnTime = proc.time()
      tryCatch({
        print('Tumor Mutational Burden')
        tumorMutationalBurden(genomesize, tmbPath, exposure_refdata_selected)
        output[['tmbPath']] = tmbPath
      }, error = function(e) {
        errors[['tmbError']] <<- e$message
        print(e)
      })
      print(paste0("Tumor Mutational Burden Runtime: ", (proc.time() - fnTime)[['elapsed']]))
    }

    # Tumor Mutational Burden separated by signatures
    if ('all' %in% fn || 'tmbSig' %in% fn) {
      fnTime = proc.time()
      tryCatch({
        print('Tumor Mutational Burden Separated by Signatures')
        mutationalSignatureBurdenSeparated(genomesize, common$cancerType, signaturePath, exposure_refdata_selected)
        output[['signaturePath']] = signaturePath
      }, error = function(e) {
        errors[['signatureError']] <<- e$message
        print(e)
      })
      print(paste0("Tumor Mutational Burden Separated by Signatures Runtime: ", (proc.time() - fnTime)[['elapsed']]))
    }

    # Mutational signature burden across cancer types
    if (length(burden)) {
      fnTime = proc.time()
      tryCatch({
        print('Mutational Signature Burden Across Cancer Types')
        mutationalSignatureBurdenAcrossCancer(burden$signatureName, genomesize, burdenPath, exposure_refdata_selected)
        output[['burdenPath']] = burdenPath
      }, error = function(e) {
        errors[['burdenError']] <<- e$message
        print(e)
      })
      print(paste0("Mutational Signature Burden Across Cancer Types Runtime: ", (proc.time() - fnTime)[['elapsed']]))
    }

    # Mutational Signature Association
    if (length(association)) {
      fnTime = proc.time()
      tryCatch({
        print('Mutational Signature Association')
        mutationalSignatureAssociation(common$useCancerType, common$cancerType, association$both, association$signatureName1, association$signatureName2, associationPath, exposure_refdata_selected)
        output[['associationPath']] = associationPath
      }, error = function(e) {
        errors[['associationError']] <<- e$message
        print(e)
      })
      print(paste0("Mutational Signature Association Runtime: ", (proc.time() - fnTime)[['elapsed']]))
    }

    # Evaluating the Performance of Mutational Signature Decomposition --------
    if ('all' %in% fn || 'decomposition' %in% fn) {
      fnTime = proc.time()
      tryCatch({
        print('Evaluating the Performance of Mutational Signature Decomposition')
        mutationalSignatureDecomposition(decompositionPath, decompositionData, exposure_refdata_selected, signature_refsets_selected, seqmatrix_refdata_selected)
        output[['decompositionPath']] = decompositionPath
        output[['decompositionData']] = decompositionData
      }, error = function(e) {
        errors[['decompositionError']] <<- e$message
        print(e)
      })
      print(paste0("Evaluating the Performance of Mutational Signature Decomposition Runtime: ", (proc.time() - fnTime)[['elapsed']]))
    }

    # Landscape of Mutational Signature Activity
    if ('all' %in% fn | 'landscape' %in% fn) {
      fnTime = proc.time()
      tryCatch({
        print('Landscape of Mutational Signature Activity')
        varDataPath = ''
        if (stringi::stri_length(landscape$variableFile) > 0) {
          varDataPath = file.path(rootDir, landscape$variableFile)
        }
        mutationalSignatureLandscape(common$cancerType, vars3Data, localData, landscapePath, exposure_refdata_selected, signature_refsets_selected, seqmatrix_refdata_selected)
        output[['landscapePath']] = landscapePath
      }, error = function(e) {
        errors[['landscapeError']] <<- e$message
        print(e)
      })
      print(paste0("Landscape of Mutational Signature Activity Runtime: ", (proc.time() - fnTime)[['elapsed']]))
    }

    # Prevalence plot
    if ('all' %in% fn | 'prevalence' %in% fn) {
      fnTime = proc.time()
      tryCatch({
        print('Prevalence of Mutational Signature')
        mutationalSignaturePrevalence(prevalence$mutation, common$cancerType, prevalencePath, exposure_refdata_selected)
        output[['prevalencePath']] = prevalencePath
      }, error = function(e) {
        errors[['prevalenceError']] <<- e$message
        print(e)
      })
      print(paste0("Prevalence of Mutational Signature Runtime: ", (proc.time() - fnTime)[['elapsed']]))
    }

    # Individual plot
    if (length(individual)) {
      fnTime = proc.time()
      tryCatch({
        print('Mutational Signature in Individual Sample')
        mutationalSignatureIndividual(individual$sample, common$cancerType, individualPath, exposure_refdata_selected, signature_refsets_selected, seqmatrix_refdata_selected)
        output[['individualPath']] = individualPath
      }, error = function(e) {
        errors[['individualError']] <<- e$message
        print(e)
      })
      print(paste0("Prevalence of Mutational Signature Runtime: ", (proc.time() - fnTime)[['elapsed']]))
    }

  }, error = function(e) {
    print(e)
  }, finally = {
    print(paste0("TotalRuntime: ", (proc.time() - totalTime)[['elapsed']]))
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output, 'errors' = errors), pretty = TRUE, auto_unbox = TRUE))
  })
}

exposureUser <- function(fn, files, common, burden = '{}', association = '{}', landscape = '{}', prevalence = '{}', individual = '{}', projectID, pythonOutput, rootDir, savePath, s3Data, localData, bucket) {
  tryCatch({
    source('services/R/Sigvisualfunc.R')
    con <- textConnection('stdout', 'wr', local = TRUE)
    sink(con, type = "message")
    sink(con, type = "output")

    output = list()
    errors = list()
    totalTime = proc.time()
    tmbPath = paste0(savePath, 'tumorMutationalBurden.svg')
    signaturePath = paste0(savePath, 'tmbSeparated.svg')
    burdenPath = paste0(savePath, 'burdenAcrossCancer.svg')
    associationPath = paste0(savePath, 'mutationalSignatureAssociation.svg')
    decompositionPath = paste0(savePath, 'mutationalSignatureDecomposition.svg')
    decompositionData = paste0(savePath, 'mutationalSignatureDecomposition.txt')
    landscapePath = paste0(savePath, 'landscapeMutationalSignature.svg')
    prevalencePath = paste0(savePath, 'prevalenceMutationalSignature.svg')
    individualPath = paste0(savePath, 'msIndividual.svg')

    # parse arguments
    common = fromJSON(common)
    burden = fromJSON(burden)
    association = fromJSON(association)
    landscape = fromJSON(landscape)
    prevalence = fromJSON(prevalence)
    individual = fromJSON(individual)
    files = fromJSON(files)

    exposure_refdata_selected <- read_delim(file.path(rootDir, files$exposureFile), delim = '\t', col_names = T)
    seqmatrix_refdata_selected <- read_delim(file.path(rootDir, files$matrixFile), delim = '\t', col_names = T)

    if (stringi::stri_length(files$signatureFile) > 0) {
      # if using user uploaded signature file
      signature_refsets_selected <- read_delim(file.path(rootDir, files$signatureFile), delim = '\t', col_names = T)
    } else {
      # else use public signature data
      s3load(paste0(s3Data, 'Signature/signature_refsets.RData'), bucket)

      signature_refsets_selected <- signature_refsets %>%
        filter(Signature_set_name == common$refSignatureSet) %>%
        select(MutationType, Signature_name, Contribution) %>%
        pivot_wider(names_from = Signature_name, values_from = Contribution)
    }

    genomesize <- genome2size(common$genome)
    cancer_type_user <- "Input"

    ## format the input as suggest by the public dataset ##
    colnames(exposure_refdata_selected)[1] <- "Sample"
    colnames(seqmatrix_refdata_selected)[1] <- "MutationType"
    colnames(signature_refsets_selected)[1] <- "MutationType"
    seqmatrix_refdata_selected <- seqmatrix_refdata_selected %>% select(1, any_of(exposure_refdata_selected$Sample)) %>% profile_format_df()
    signature_refsets_selected <- signature_refsets_selected %>% select(MutationType, any_of(colnames(exposure_refdata_selected))) %>% profile_format_df()

    exposure_refdata_selected <- exposure_refdata_selected %>% pivot_longer(cols = -Sample, names_to = "Signature_name", values_to = "Exposure") %>% mutate(Cancer_Type = cancer_type_user)
    signature_refsets_selected <- signature_refsets_selected %>% select(-Type, - SubType) %>% pivot_longer(cols = -MutationType, names_to = "Signature_name", values_to = "Contribution")
    seqmatrix_refdata_selected <- seqmatrix_refdata_selected %>% select(-Type, - SubType) %>% pivot_longer(cols = -MutationType, names_to = "Sample", values_to = "Mutations") %>% mutate(Cancer_Type = cancer_type_user)


    ## Tumor Overall Mutational Burden
    if ('all' %in% fn || 'tmb' %in% fn) {
      fnTime = proc.time()
      tryCatch({
        print('Tumor Mutational Burden')
        tumorMutationalBurden(genomesize, tmbPath, exposure_refdata_selected)
        output[['tmbPath']] = tmbPath
      }, error = function(e) {
        errors[['tmbError']] <<- e$message
        print(e)
      })
      print(paste0("Tumor Mutational Burden Runtime: ", (proc.time() - fnTime)[['elapsed']]))
    }

    # Tumor Mutational Burden separated by signatures
    if ('all' %in% fn || 'tmbSig' %in% fn) {
      fnTime = proc.time()
      tryCatch({
        print('Tumor Mutational Burden Separated by Signatures')
        mutationalSignatureBurdenSeparated(genomesize, cancer_type_user, signaturePath, exposure_refdata_selected)
        output[['signaturePath']] = signaturePath
      }, error = function(e) {
        errors[['signatureError']] <<- e$message
        print(e)
      })
      print(paste0("Tumor Mutational Burden Separated by Signatures Runtime: ", (proc.time() - fnTime)[['elapsed']]))
    }
    # Mutational signature burden across cancer types
    if (length(burden)) {
      fnTime = proc.time()
      tryCatch({
        print('Mutational Signature Burden Across Cancer Types')
        mutationalSignatureBurdenAcrossCancer(burden$signatureName, genomesize, burdenPath, exposure_refdata_selected)
        output[['burdenPath']] = burdenPath
      }, error = function(e) {
        errors[['burdenError']] <<- e$message
        print(e)
      })
      print(paste0("Mutational Signature Burden Across Cancer Types Runtime: ", (proc.time() - fnTime)[['elapsed']]))
    }
    # Mutational Signature Association
    if (length(association)) {
      fnTime = proc.time()
      tryCatch({
        print('Mutational Signature Association')
        mutationalSignatureAssociation(common$useCancerType, cancer_type_user, association$both, association$signatureName1, association$signatureName2, associationPath, exposure_refdata_selected)
        output[['associationPath']] = associationPath
      }, error = function(e) {
        errors[['associationError']] <<- e$message
        print(e)
      })
      print(paste0("Mutational Signature Association Runtime: ", (proc.time() - fnTime)[['elapsed']]))
    }
    # Evaluating the Performance of Mutational Signature Decomposition --------
    if ('all' %in% fn || 'decomposition' %in% fn) {
      fnTime = proc.time()
      tryCatch({
        print('Evaluating the Performance of Mutational Signature Decomposition')
        mutationalSignatureDecomposition(decompositionPath, decompositionData, exposure_refdata_selected, signature_refsets_selected, seqmatrix_refdata_selected)
        output[['decompositionPath']] = decompositionPath
        output[['decompositionData']] = decompositionData
      }, error = function(e) {
        errors[['decompositionError']] <<- e$message
        print(e)
      })
      print(paste0("Evaluating the Performance of Mutational Signature Decomposition Runtime: ", (proc.time() - fnTime)[['elapsed']]))
    }
    # Landscape of Mutational Signature Activity
    if ('all' %in% fn | 'landscape' %in% fn) {
      fnTime = proc.time()
      tryCatch({
        print('Landscape of Mutational Signature Activity')
        varDataPath = ''
        if (stringi::stri_length(landscape$variableFile) > 0) {
          varDataPath = file.path(rootDir, landscape$variableFile$signatureFile)
        }
        mutationalSignatureLandscape(cancer_type_user, vars3Data, localData, landscapePath, exposure_refdata_selected, signature_refsets_selected, seqmatrix_refdata_selected)
        output[['landscapePath']] = landscapePath
      }, error = function(e) {
        errors[['landscapeError']] <<- e$message
        print(e)
      })
      print(paste0("Landscape of Mutational Signature Activity Runtime: ", (proc.time() - fnTime)[['elapsed']]))
    }
    # Prevalence plot
    if ('all' %in% fn | 'prevalence' %in% fn) {
      fnTime = proc.time()
      tryCatch({
        print('Prevalence of Mutational Signature')
        mutationalSignaturePrevalence(prevalence$mutation, cancer_type_user, prevalencePath, exposure_refdata_selected)
        output[['prevalencePath']] = prevalencePath
      }, error = function(e) {
        errors[['prevalenceError']] <<- e$message
        print(e)
      })
      print(paste0("Prevalence of Mutational Signature Runtime: ", (proc.time() - fnTime)[['elapsed']]))
    }

    # Individual plot
    if (length(individual)) {
      fnTime = proc.time()
      tryCatch({
        print('Mutational Signature in Individual Sample')
        mutationalSignatureIndividual(individual$sample, common$cancerType, individualPath, exposure_refdata_selected, signature_refsets_selected, seqmatrix_refdata_selected)
        output[['individualPath']] = individualPath
      }, error = function(e) {
        errors[['individualError']] <<- e$message
        print(e)
      })
      print(paste0("Prevalence of Mutational Signature Runtime: ", (proc.time() - fnTime)[['elapsed']]))
    }

  }, error = function(e) {
    print(e)
  }, finally = {
    print(paste0("Total Runtime: ", (proc.time() - totalTime)[['elapsed']]))
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output, 'errors' = errors), pretty = TRUE, auto_unbox = TRUE))
  })
}
