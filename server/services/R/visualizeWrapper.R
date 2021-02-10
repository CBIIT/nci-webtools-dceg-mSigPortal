library(tidyverse)
library(cowplot)
library(ggsci)
library(ggrepel)
library(factoextra)
library(scales)
library(ggpubr)
library(jsonlite)

# get all Reference Signature Set options using profile_name (profile type and matrix size)
getReferenceSignatureSets <- function(profileType, dataPath) {
  load(paste0(dataPath, 'Signature/signature_refsets.RData'))

  profile_name <- if_else(profileType == "SBS", "SBS96", if_else(profileType == "DBS", "DBS78", if_else(profileType == "ID", "ID83", NA_character_)))
  signatureSets <- signature_refsets %>% filter(Profile == profile_name) %>% pull(Signature_set_name) %>% unique()

  return(signatureSets)
}

# get list of signatures in the selected signature set
getSignaturesR <- function(profileType, signatureSetName, dataPath) {
  load(paste0(dataPath, 'Signature/signature_refsets.RData'))
  profile_name <- if_else(profileType == "SBS", "SBS96", if_else(profileType == "DBS", "DBS78", if_else(profileType == "ID", "ID83", NA_character_)))
  signature_refsets_input <- signature_refsets %>% filter(Profile == profile_name, Signature_set_name == signatureSetName)
  refsig <- signature_refsets_input %>%
      select(Signature_name, MutationType, Contribution) %>%
      pivot_wider(names_from = Signature_name, values_from = Contribution)
  signatures <- colnames(refsig[1, -1])

  return(signatures)
}

# get list of options for study, cancer type, and experimental strategy
getPublicDataOptions <- function(dataPath) {
  load(paste0(dataPath, 'Seqmatrix/seqmatrix_refdata_info.RData'))
  return(toJSON(seqmatrix_refdata_info2, pretty = TRUE, auto_unbox = TRUE))
}

# get public svgFiles and load into session
getPublicData <- function(study, cancerType, experimentalStrategy, dataPath) {
  load(paste0(dataPath, 'Seqmatrix/seqmatrix_refdata_info.RData'))

  svgfiles <- seqmatrix_refdata_info %>% mutate(Path = paste0(dataPath, 'Seqmatrix/', Path))
  if (cancerType == 'PanCancer') {
    svgfiles_public <- svgfiles %>% filter(Study == study, Dataset == experimentalStrategy)
  } else {
    svgfiles_public <- svgfiles %>% filter(Study == study, Cancer_Type == cancerType, Dataset == experimentalStrategy)
  }

  return(toJSON(svgfiles_public, pretty = TRUE, auto_unbox = TRUE))
}

# Tumor Mutation Burden ---------------------------------------------------
profilerSummary <- function(matrixList, projectID, pythonOutput, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    plotPath = paste0(savePath, 'profilerSummary.svg')

    matrixJSON = fromJSON(matrixList)
    matrixList = data.frame(matrixJSON$data)
    colnames(matrixList) = matrixJSON$columns

    data_input <- tibble()

    for (i in 1:dim(matrixList)[1]) {
      matrixfile_selected <- matrixList$Path[i]
      data_input_tmp <- read_delim(matrixfile_selected, delim = '\t')
      if (dim(data_input_tmp)[1] > 0) {
        data_input_tmp <- data_input_tmp %>% pivot_longer(cols = -MutationType) %>%
        group_by(name) %>%
        summarise(Mutations = sum(value, na.rm = TRUE)) %>%
        mutate(Profile = paste0(matrixList$Profile_Type[i], matrixList$Matrix_Size[i])) %>%
        select(Sample = name, Profile, Mutations)
        data_input <- bind_rows(data_input, data_input_tmp)
      }
    }
    profile_heatmap_plot(data = data_input, output_plot = plotPath)

    output = list('plotPath' = plotPath)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

profilerSummaryPublic <- function(study, cancerType, experimentalStrategy, projectID, pythonOutput, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'Seqmatrix/seqmatrix_refdata_subset_files.RData'))
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    plotPath = paste0(savePath, 'profilerSummaryPublic.svg')

    publicDataFile <- seqmatrix_refdata_subset_files %>% filter(Study == study, Cancer_Type == cancerType, Dataset == experimentalStrategy) %>% pull(file)
    seqmatrix_refdata_public <- get(load(paste0(dataPath, 'Seqmatrix/', publicDataFile)))

    data_input <- seqmatrix_refdata_public %>%
      group_by(Sample, Profile) %>%
      summarise(Mutations = sum(Mutations, na.rm = TRUE)) %>%
      ungroup()

    profile_heatmap_plot(data = data_input, output_plot = plotPath)

    output = list('plotPath' = plotPath)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

### Cosine Similarity tab ###
# section 1: Cosine similarity within samples #
cosineSimilarityWithin <- function(matrixFile, projectID, pythonOutput, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    plotPath = paste0(savePath, 'cos_sim_within.svg')
    txtPath = paste0(savePath, 'cos_sim_within.txt')

    data_input <- read_delim(matrixFile, delim = '\t')
    data_input <- data_input %>% select_if(~!is.numeric(.) || sum(.) > 0)
    # Heatmap of cosine similarity within samples  and put on the web---------------------------
    cos_sim_res1 = cos_sim_df(data_input, data_input)
    plot_cosine_heatmap_df(cos_sim_res1, cluster_rows = TRUE, plot_values = FALSE, output_plot = plotPath)
    cos_sim_res1 %>% write_delim(txtPath, delim = '\t', col_names = T)

    output = list('plotPath' = plotPath, 'txtPath' = txtPath)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

cosineSimilarityWithinPublic <- function(profileType, matrixSize, study, cancerType, experimentalStrategy, projectID, pythonOutput, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'Seqmatrix/seqmatrix_refdata_subset_files.RData'))
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    plotPath = paste0(savePath, 'cos_sim_within.svg')
    txtPath = paste0(savePath, 'cos_sim_within.txt')

    publicDataFile <- seqmatrix_refdata_subset_files %>% filter(Study == study, Cancer_Type == cancerType, Dataset == experimentalStrategy) %>% pull(file)
    seqmatrix_refdata_public <- get(load(paste0(dataPath, 'Seqmatrix/', publicDataFile)))

    data_input <- seqmatrix_refdata_public %>% filter(Profile == paste0(profileType, matrixSize)) %>%
      select(MutationType, Sample, Mutations) %>%
      pivot_wider(id_cols = MutationType, names_from = Sample, values_from = Mutations)

    data_input <- data_input %>% select_if(~!is.numeric(.) || sum(.) > 0)
    # Heatmap of cosine similarity within samples  and put on the web---------------------------
    cos_sim_res1 = cos_sim_df(data_input, data_input)
    plot_cosine_heatmap_df(cos_sim_res1, cluster_rows = TRUE, plot_values = FALSE, output_plot = plotPath)
    cos_sim_res1 %>% write_delim(txtPath, delim = '\t', col_names = T)

    output = list('plotPath' = plotPath, 'txtPath' = txtPath)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

# section 2: Cosine similarity  to reference signatures 
# Two parameters need: Profile Type, Reference Signature Set
# Profile Type only support SBS, DBS, ID
cosineSimilarityRefSig <- function(profileType, signatureSetName, matrixFile, projectID, pythonOutput, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'Signature/signature_refsets.RData'))
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    plotPath = paste0(savePath, 'cos_sim_refsig.svg')
    txtPath = paste0(savePath, 'cos_sim_refsig.txt')

    profile_name <- if_else(profileType == "SBS", "SBS96", if_else(profileType == "DBS", "DBS78", if_else(profileType == "ID", "ID83", NA_character_)))
    signature_refsets_data <- signature_refsets %>% filter(Profile == profile_name, Signature_set_name == signatureSetName)
    refsig <- signature_refsets_data %>%
      select(Signature_name, MutationType, Contribution) %>%
      pivot_wider(names_from = Signature_name, values_from = Contribution)

    data_input <- read_delim(matrixFile, delim = '\t')
    data_input <- data_input %>% select_if(~!is.numeric(.) || sum(.) > 0)

    # Heatmap of cosine similarity to reference set signature and put on the web---------------------------
    cos_sim_res2 = cos_sim_df(data_input, refsig)
    plot_cosine_heatmap_df(cos_sim_res2, cluster_rows = TRUE, plot_values = FALSE, output_plot = plotPath)
    cos_sim_res2 %>% write_delim(txtPath, delim = '\t', col_names = T)

    output = list('plotPath' = plotPath, 'txtPath' = txtPath)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

cosineSimilarityRefSigPublic <- function(profileType, signatureSetName, study, cancerType, experimentalStrategy, projectID, pythonOutput, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'Seqmatrix/seqmatrix_refdata_subset_files.RData'))
  load(paste0(dataPath, 'Signature/signature_refsets.RData'))
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    plotPath = paste0(savePath, 'cos_sim_refsig.svg')
    txtPath = paste0(savePath, 'cos_sim_refsig.txt')

    profile_name <- if_else(profileType == "SBS", "SBS96", if_else(profileType == "DBS", "DBS78", if_else(profileType == "ID", "ID83", NA_character_)))
    signature_refsets_data <- signature_refsets %>% filter(Profile == profile_name, Signature_set_name == signatureSetName)
    refsig <- signature_refsets_data %>%
      select(Signature_name, MutationType, Contribution) %>%
      pivot_wider(names_from = Signature_name, values_from = Contribution)

    publicDataFile <- seqmatrix_refdata_subset_files %>% filter(Study == study, Cancer_Type == cancerType, Dataset == experimentalStrategy) %>% pull(file)
    seqmatrix_refdata_public <- get(load(paste0(dataPath, 'Seqmatrix/', publicDataFile)))

    data_input <- seqmatrix_refdata_public %>% filter(Profile == profile_name) %>%
      select(MutationType, Sample, Mutations) %>%
      pivot_wider(id_cols = MutationType, names_from = Sample, values_from = Mutations)
    data_input <- data_input %>% select_if(~!is.numeric(.) || sum(.) > 0)

    # Heatmap of cosine similarity to reference set signature and put on the web---------------------------
    cos_sim_res2 = cos_sim_df(data_input, refsig)
    plot_cosine_heatmap_df(cos_sim_res2, cluster_rows = TRUE, plot_values = FALSE, output_plot = plotPath)
    cos_sim_res2 %>% write_delim(txtPath, delim = '\t', col_names = T)

    output = list('plotPath' = plotPath, 'txtPath' = txtPath)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

# section 3: Cosine similarity to Public data -----------------------------
# find the common profile between data and seqmatrix
cosineSimilarityPublic <- function(matrixFile, study, cancerType, profileName, projectID, pythonOutput, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'Seqmatrix/seqmatrix_refdata_subset_files.RData'))
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    plotPath = paste0(savePath, 'cos_sim_public.svg')
    txtPath = paste0(savePath, 'cos_sim_public.txt')


    ## input data
    data_input <- read_delim(matrixFile, delim = '\t')
    data_input <- data_input %>% select_if(~!is.numeric(.) || sum(.) > 0)

    ## seqmatrix data from public data
    publicDataFile <- seqmatrix_refdata_subset_files %>% filter(Study == study, Cancer_Type == cancerType) %>% pull(file)
    publicData <- get(load(paste0(dataPath, 'Seqmatrix/', publicDataFile)))

    sigmatrix_data <- publicData %>% filter(Profile == profileName) %>%
    select(Sample, MutationType, Mutations) %>%
    pivot_wider(names_from = Sample, values_from = Mutations)

    # Heatmap of cosine similarity to public seqmatrix data and put on the web---------------------------
    cos_sim_res3 = cos_sim_df(data_input, sigmatrix_data)
    plot_cosine_heatmap_df(cos_sim_res3, cluster_rows = TRUE, plot_values = FALSE, output_plot = plotPath)
    cos_sim_res3 %>% write_delim(txtPath, delim = '\t', col_names = T)

    output = list('plotPath' = plotPath, 'txtPath' = txtPath)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}


### Profile Comparison tab ###
# section 1: Cosine similarity within samples 
# three parameters need: Profile Type, Sample Name1 and Sample Name2
profileComparisonWithin <- function(profileType, sampleName1, sampleName2, matrixFile, projectID, pythonOutput, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    plotPath = paste0(savePath, 'pro_com_within.svg')
    error = ''

    data_input <- read_delim(matrixFile, delim = '\t')
    data_input <- data_input %>% select_if(~!is.numeric(.) || sum(.) > 0)

    profile1 <- data_input %>% select(MutationType, one_of(sampleName1))
    profile2 <- data_input %>% select(MutationType, one_of(sampleName2))

    if (dim(profile1)[2] != 2) {
      error = paste0('Sample: ', sampleName1, ' has no Data. Try a different sample.')
      stop(error)
    }
    if (dim(profile2)[2] != 2) {
      error = paste0('Sample: ', sampleName2, ' has no Data. Try a different sample.')
      stop(error)
    }

    plot_compare_profiles_diff(profile1, profile2, condensed = FALSE, output_plot = plotPath)

    output = list('plotPath' = plotPath)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output, 'error' = error), pretty = TRUE, auto_unbox = TRUE))
  })
}

profileComparisonWithinPublic <- function(profileType, sampleName1, sampleName2, study, cancerType, experimentalStrategy, projectID, pythonOutput, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'Seqmatrix/seqmatrix_refdata_subset_files.RData'))
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    plotPath = paste0(savePath, 'pro_com_within.svg')

    matrix_size <- if_else(profileType == "SBS", "96", if_else(profileType == "DBS", "78", if_else(profileType == "ID", "83", NA_character_)))
    profile_name <- paste0(profileType, matrix_size)

    publicDataFile <- seqmatrix_refdata_subset_files %>% filter(Study == study, Cancer_Type == cancerType, Dataset == experimentalStrategy) %>% pull(file)
    seqmatrix_refdata_public <- get(load(paste0(dataPath, 'Seqmatrix/', publicDataFile)))

    data_input <- seqmatrix_refdata_public %>% filter(Profile == profile_name) %>%
      select(MutationType, Sample, Mutations) %>%
      filter(Sample %in% c(sampleName1, sampleName2)) %>%
      pivot_wider(id_cols = MutationType, names_from = Sample, values_from = Mutations)

    profile1 <- data_input %>% select(MutationType, one_of(sampleName1))
    profile2 <- data_input %>% select(MutationType, one_of(sampleName2))
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

# section 2: Comparison to reference signatures # 
# four parameters need: “Profile Type”, “Sample Name”, “Reference Signature Set” and “Compare Single Signature or Combined Signatures” # 
profileComparisonRefSig <- function(profileType, sampleName, signatureSetName, compare, matrixFile, projectID, pythonOutput, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'Signature/signature_refsets.RData'))
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    plotPath = paste0(savePath, 'pro_com_refsig.svg')
    error = ''

    profile_name <- if_else(profileType == "SBS", "SBS96", if_else(profileType == "DBS", "DBS78", if_else(profileType == "ID", "ID83", NA_character_)))

    signature_refsets_input <- signature_refsets %>% filter(Profile == profile_name, Signature_set_name == signatureSetName)
    refsig <- signature_refsets_input %>%
      select(Signature_name, MutationType, Contribution) %>%
      pivot_wider(names_from = Signature_name, values_from = Contribution)

    if (str_detect(compare, ";")) {
      profile2 <- signature_sum_operation(sigdatabase = signature_refsets_input, sigsetname = signatureSetName, formulax = compare, outsigname = "Reconstructed")
      profile_names = c("Original", "Reconstructed")
    } else {
      profile2 <- refsig %>% select(MutationType, one_of(compare))
    }

    data_input <- read_delim(matrixFile, delim = '\t')
    data_input <- data_input %>% select_if(~!is.numeric(.) || sum(.) > 0)

    profile1 <- data_input %>% select(MutationType, one_of(sampleName))

    if (dim(profile1)[2] != 2) {
      error = paste0('Sample: ', sampleName, ' has no Data. Try a different sample.')
      stop(error)
    }

    profile_names = c(colnames(profile1)[2], colnames(profile2)[2])
    plot_compare_profiles_diff(profile1, profile2, condensed = FALSE, profile_names = profile_names, output_plot = plotPath)

    output = list('plotPath' = plotPath)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output, 'error' = error), pretty = TRUE, auto_unbox = TRUE))
  })
}

profileComparisonRefSigPublic <- function(profileType, sampleName, signatureSetName, compare, study, cancerType, experimentalStrategy, projectID, pythonOutput, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'Signature/signature_refsets.RData'))
  load(paste0(dataPath, 'Seqmatrix/seqmatrix_refdata_subset_files.RData'))
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    plotPath = paste0(savePath, 'pro_com_refsig.svg')

    profile_name <- if_else(profileType == "SBS", "SBS96", if_else(profileType == "DBS", "DBS78", if_else(profileType == "ID", "ID83", NA_character_)))

    signature_refsets_input <- signature_refsets %>% filter(Profile == profile_name, Signature_set_name == signatureSetName)
    refsig <- signature_refsets_input %>%
      select(Signature_name, MutationType, Contribution) %>%
      pivot_wider(names_from = Signature_name, values_from = Contribution)

    if (str_detect(compare, ";")) {
      profile2 <- signature_sum_operation(sigdatabase = signature_refsets_input, sigsetname = signatureSetName, formulax = compare, outsigname = "Reconstructed")
      profile_names = c("Original", "Reconstructed")
    } else {
      profile2 <- refsig %>% select(MutationType, one_of(compare))
    }
    publicDataFile <- seqmatrix_refdata_subset_files %>% filter(Study == study, Cancer_Type == cancerType, Dataset == experimentalStrategy) %>% pull(file)
    seqmatrix_refdata_public <- get(load(paste0(dataPath, 'Seqmatrix/', publicDataFile)))

    profile1 <- seqmatrix_refdata_public %>% filter(Profile == profile_name) %>%
      select(MutationType, Sample, Mutations) %>%
      filter(Sample %in% c(sampleName)) %>%
      pivot_wider(id_cols = MutationType, names_from = Sample, values_from = Mutations)
    profile_names = c(colnames(profile1)[2], colnames(profile2)[2])
    plot_compare_profiles_diff(profile1, profile2, condensed = FALSE, profile_names = profile_names, output_plot = plotPath)

    output = list('plotPath' = plotPath)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

# section 3: Profile Comparison to Public data ----------------------------
profileComparisonPublic <- function(profileName, matrixFile, userSample, study, cancerType, publicSample, projectID, pythonOutput, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'Seqmatrix/seqmatrix_refdata_subset_files.RData'))
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    plotPath = paste0(savePath, 'pro_com_public.svg')

    ## input data
    data_input <- read_delim(matrixFile, delim = '\t')
    data_input <- data_input %>% select_if(~!is.numeric(.) || sum(.) > 0)
    profile1 <- data_input %>% select(MutationType, one_of(userSample))

    ## seqmatrix data from public data
    publicDataFile <- seqmatrix_refdata_subset_files %>% filter(Study == study, Cancer_Type == cancerType) %>% pull(file)
    publicData <- get(load(paste0(dataPath, 'Seqmatrix/', publicDataFile)))

    profile2 <- publicData %>% filter(Profile == profileName) %>%
      select(Sample, MutationType, Mutations) %>%
      filter(Sample == publicSample) %>%
      pivot_wider(names_from = Sample, values_from = Mutations)

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

# Motif analysis ----------------------------------------------------------
###parameters:
mutationalPattern <- function(matrixFile, proportion, pattern, projectID, pythonOutput, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    plotPath = paste0(savePath, 'mpea.svg')
    txtPath = paste0(savePath, 'mpea.txt')

    data_input <- read_delim(matrixFile, delim = '\t')
    data_input <- data_input %>% select_if(~!is.numeric(.) || sum(.) > 0) %>%
      pivot_longer(cols = -MutationType) %>%
      mutate(Study = "Input") %>%
      select(Study, Sample = name, MutationType, Mutations = value) %>%
      mutate(Type = str_sub(MutationType, 3, 5), SubType1 = str_sub(MutationType, 1, 1), SubType2 = str_sub(str_sub(MutationType, 7, 7))) %>%
      select(Study, Sample, MutationType, Type, SubType1, SubType2, Mutations)

    content_all_tmp <- content_extraction(data_input)
    data_tmp <- content_all_tmp %>%
    filter(N1 > proportion) %>%
    count(Pattern, sort = T) %>%
    mutate(Type = "Frequency of Mutational Pattern") %>% select(Type, Pattern, n)


    context_plot(data = data_input, pattern = pattern, output_plot = plotPath)
    context_plot(data = data_input, pattern = pattern, data_return = TRUE) %>%
    arrange(desc(N1)) %>% write_delim(txtPath, delim = '\t', col_names = T)
    if (dim(data_tmp) > 0) {
      barPath = paste0(savePath, 'barchart.svg')
      barchart_plot2(data = data_tmp, plot_width = 16, plot_height = 5, output_plot = barPath)
      output = list('barPath' = barPath, 'plotPath' = plotPath, 'txtPath' = txtPath)
    } else {
      print(paste0('No mutational pattern with proportion of mutations large than', proportion))
      output = list('plotPath' = plotPath, 'txtPath' = txtPath)
    }

  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

# Motif analysis ----------------------------------------------------------
###parameters:
mutationalPatternPublic <- function(study, cancerType, experimentalStrategy, proportion, pattern, projectID, pythonOutput, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'Others/content_data_all.RData'))
  load(paste0(dataPath, 'Seqmatrix/seqmatrix_refdata_subset_files.RData'))
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    barPath = paste0(savePath, 'barchart.svg')
    plotPath = paste0(savePath, 'mpea.svg')
    txtPath = paste0(savePath, 'mpea.txt')

    data_tmp <- content_data_all %>%
      filter(N1 > proportion, str_detect(Study, paste0("^", study, "@"))) %>%
      count(Pattern, sort = T) %>%
      mutate(Type = "Frequency of Mutational Pattern") %>% select(Type, Pattern, n)

    barchart_plot2(data = data_tmp, plot_width = 16, plot_height = 5, output_plot = barPath)

    if (cancerType != "PanCancer") {
      publicDataFile <- seqmatrix_refdata_subset_files %>% filter(Study == study, Cancer_Type == cancerType, Dataset == experimentalStrategy) %>% pull(file)
      seqmatrix_refdata <- get(load(paste0(dataPath, 'Seqmatrix/', publicDataFile))) %>% filter(Study == study)
    } else {
      seqmatrix_refdata <- get(load(paste0(dataPath, 'Seqmatrix/seqmatrix_refdata.RData'))) %>% filter(Study == study)
    }
    data_input <- seqmatrix_refdata %>%
      filter(Profile == "SBS96") %>%
      mutate(Study = paste0(Study, "@", Cancer_Type)) %>%
      select(Study, Sample, MutationType, Mutations) %>%
      mutate(Type = str_sub(MutationType, 3, 5), SubType1 = str_sub(MutationType, 1, 1), SubType2 = str_sub(str_sub(MutationType, 7, 7))) %>%
      select(Study, Sample, MutationType, Type, SubType1, SubType2, Mutations)

    context_plot(data = data_input, pattern = pattern, output_plot = plotPath)
    context_plot(data = data_input, pattern = pattern, data_return = TRUE) %>%
    arrange(desc(N1)) %>% write_delim(txtPath, delim = '\t', col_names = T)

    output = list('barPath' = barPath, 'plotPath' = plotPath, 'txtPath' = txtPath)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

### “Principal Component Analysis” ###
# Two parameters need: Profile Type, Reference Signature Set
# Profile Type only support SBS, DBS, ID
pca <- function(profileType, signatureSetName, matrixFile, projectID, pythonOutput, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'Signature/signature_refsets.RData'))
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    pca1 = paste0(savePath, 'pca1.svg')
    pca2 = paste0(savePath, 'pca2.svg')
    pca3 = paste0(savePath, 'pca3.svg')
    heatmap = paste0(savePath, 'heatmap.svg')
    pca2Data = paste0(savePath, 'pca2_data.txt')
    pca3Data = paste0(savePath, 'pca3_data.txt')
    heatmapData = paste0(savePath, 'heatmap_data.txt')



    profile_name <- if_else(profileType == "SBS", "SBS96", if_else(profileType == "DBS", "DBS78", if_else(profileType == "ID", "ID83", NA_character_)))

    data_input <- read_delim(matrixFile, delim = '\t')
    data_input <- data_input %>% select_if(~!is.numeric(.) || sum(.) > 0)

    # PCA plot ----------------------------------------------------------------
    mdata_input <- as.matrix(data_input[, -1])
    rownames(mdata_input) <- data_input$MutationType

    res.pca <- prcomp(t(mdata_input), scale = FALSE, center = FALSE)

    xleng <- dim(res.pca$x)[2] * 0.25 + 1
    xleng <- if_else(xleng > 4, 4, xleng)
    yleng <- 4
    pcap1 <- fviz_eig(res.pca, ncp = 10, main = "", addlabels = T)
    ggsave(filename = pca1, plot = pcap1, width = xleng, height = yleng)

    if ((dim(data_input)[2] - 1) > 35) {
      pcap2 <- fviz_pca_ind(res.pca, axes = c(1, 2), geom = "point", col.ind = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
    } else {
      pcap2 <- fviz_pca_ind(res.pca, axes = c(1, 2), col.ind = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
    }
    ggsave(filename = pca2, plot = pcap2, width = 10, height = 7)

    ## a link to download pca data1
    res.pca$x %>% as.data.frame() %>% rownames_to_column(var = 'Sample') %>% write_delim(pca2Data, delim = '\t', col_names = T)

    pcap3 <- fviz_pca_var(res.pca, axes = c(1, 2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE # Avoid text overlapping
    )
    ggsave(filename = pca3, plot = pcap3, width = 10, height = 7)


    # a link to download pca data2
    res.pca$rotation %>% as.data.frame() %>% rownames_to_column(var = 'Sample') %>% write_delim(pca3Data, delim = '\t', col_names = T)

    # parameters:  Reference Signature Set
    #all the option of Reference Signature Set
    ## only work if the profile included in singatures dataset
    if (profile_name %in% unique(signature_refsets$Profile)) {
      signature_refsets_input <- signature_refsets %>% filter(Profile == profile_name, Signature_set_name == signatureSetName)
      sigref_data <- signature_refsets_input %>%
        select(Signature_name, MutationType, Contribution) %>%
        pivot_wider(names_from = Signature_name, values_from = Contribution)

      ## heatmap between PCs and signatures
      sigpca_data <- res.pca$rotation %>% as.data.frame() %>% rownames_to_column(var = "MutationType")
      cos_sim_res4 = cos_sim_df(sigpca_data, sigref_data)
      # put this heatmap on the web
      plot_cosine_heatmap_df(cos_sim_res4, cluster_rows = TRUE, plot_values = FALSE, output_plot = heatmap)
      # a link to download the cosine similarity bellow the plot
      # you could rename the file name if you need
      cos_sim_res4 %>% write_delim(heatmapData, delim = '\t', col_names = T)
    }

    output = list(
      'pca1' = pca1,
      'pca2' = pca2,
      'pca3' = pca3,
      'heatmap' = heatmap,
      'pca2Data' = pca2Data,
      'pca3Data' = pca3Data,
      'heatmapData' = heatmapData
     )
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

pcaPublic <- function(profileType, signatureSetName, study, cancerType, experimentalStrategy, projectID, pythonOutput, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'Signature/signature_refsets.RData'))
  load(paste0(dataPath, 'Seqmatrix/seqmatrix_refdata_subset_files.RData'))
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    pca1 = paste0(savePath, 'pca1.svg')
    pca2 = paste0(savePath, 'pca2.svg')
    pca3 = paste0(savePath, 'pca3.svg')
    heatmap = paste0(savePath, 'heatmap.svg')
    pca2Data = paste0(savePath, 'pca2_data.txt')
    pca3Data = paste0(savePath, 'pca3_data.txt')
    heatmapData = paste0(savePath, 'heatmap_data.txt')



    profile_name <- if_else(profileType == "SBS", "SBS96", if_else(profileType == "DBS", "DBS78", if_else(profileType == "ID", "ID83", NA_character_)))
    matrix_size <- if_else(profileType == "SBS", "96", if_else(profileType == "DBS", "78", if_else(profileType == "ID", "83", NA_character_)))

    publicDataFile <- seqmatrix_refdata_subset_files %>% filter(Study == study, Cancer_Type == cancerType, Dataset == experimentalStrategy) %>% pull(file)
    seqmatrix_refdata_public <- get(load(paste0(dataPath, 'Seqmatrix/', publicDataFile)))

    data_input <- seqmatrix_refdata_public %>% filter(Profile == profile_name) %>%
      select(MutationType, Sample, Mutations) %>%
      pivot_wider(id_cols = MutationType, names_from = Sample, values_from = Mutations)
    data_input <- data_input %>% select_if(~!is.numeric(.) || sum(.) > 0)

    # PCA plot ----------------------------------------------------------------
    mdata_input <- as.matrix(data_input[, -1])
    rownames(mdata_input) <- data_input$MutationType

    res.pca <- prcomp(t(mdata_input), scale = FALSE, center = FALSE)

    xleng <- dim(res.pca$x)[2] * 0.25 + 1
    xleng <- if_else(xleng > 4, 4, xleng)
    yleng <- 4
    pcap1 <- fviz_eig(res.pca, ncp = 10, main = "", addlabels = T)
    ggsave(filename = pca1, plot = pcap1, width = xleng, height = yleng)

    if ((dim(data_input)[2] - 1) > 35) {
      pcap2 <- fviz_pca_ind(res.pca, axes = c(1, 2), geom = "point", col.ind = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
    } else {
      pcap2 <- fviz_pca_ind(res.pca, axes = c(1, 2), col.ind = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
    }
    ggsave(filename = pca2, plot = pcap2, width = 10, height = 7)

    ## a link to download pca data1
    res.pca$x %>% as.data.frame() %>% rownames_to_column(var = 'Sample') %>% write_delim(pca2Data, delim = '\t', col_names = T)

    pcap3 <- fviz_pca_var(res.pca, axes = c(1, 2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE # Avoid text overlapping
    )
    ggsave(filename = pca3, plot = pcap3, width = 10, height = 7)


    # a link to download pca data2
    res.pca$rotation %>% as.data.frame() %>% rownames_to_column(var = 'Sample') %>% write_delim(pca3Data, delim = '\t', col_names = T)

    # parameters:  Reference Signature Set
    #all the option of Reference Signature Set
    ## only work if the profile included in singatures dataset
    if (profile_name %in% unique(signature_refsets$Profile)) {
      signature_refsets_input <- signature_refsets %>% filter(Profile == profile_name, Signature_set_name == signatureSetName)
      sigref_data <- signature_refsets_input %>%
        select(Signature_name, MutationType, Contribution) %>%
        pivot_wider(names_from = Signature_name, values_from = Contribution)

      ## heatmap between PCs and signatures
      sigpca_data <- res.pca$rotation %>% as.data.frame() %>% rownames_to_column(var = "MutationType")
      cos_sim_res4 = cos_sim_df(sigpca_data, sigref_data)
      # put this heatmap on the web
      plot_cosine_heatmap_df(cos_sim_res4, cluster_rows = TRUE, plot_values = FALSE, output_plot = heatmap)
      # a link to download the cosine similarity bellow the plot
      # you could rename the file name if you need
      cos_sim_res4 %>% write_delim(heatmapData, delim = '\t', col_names = T)
    }

    output = list(
      'pca1' = pca1,
      'pca2' = pca2,
      'pca3' = pca3,
      'heatmap' = heatmap,
      'pca2Data' = pca2Data,
      'pca3Data' = pca3Data,
      'heatmapData' = heatmapData
     )
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

#  section 2 PCA together with public data --------------------------------
pcaWithPublic <- function(matrixFile, study, cancerType, profileName, projectID, pythonOutput, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'Signature/signature_refsets.RData'))
  load(paste0(dataPath, 'Seqmatrix/seqmatrix_refdata_subset_files.RData'))
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    pca1 = paste0(savePath, 'pca1_with_public_.svg')
    pca2 = paste0(savePath, 'pca2_with_public_.svg')
    pca3 = paste0(savePath, 'pca3_with_public_.svg')
    pca2Data = paste0(savePath, 'pca2_data_with_public_.txt')
    pca3Data = paste0(savePath, 'pca3_data_with_public_.txt')


    data_input1 <- read_delim(matrixFile, delim = '\t')
    data_input1 <- data_input1 %>% select_if(~!is.numeric(.) || sum(.) > 0)

    ## seqmatrix data from public data
    publicDataFile <- seqmatrix_refdata_subset_files %>% filter(Study == study, Cancer_Type == cancerType) %>% pull(file)
    publicData <- get(load(paste0(dataPath, 'Seqmatrix/', publicDataFile)))

    data_input2 <- publicData %>% filter(Profile == profileName) %>%
      select(Sample, MutationType, Mutations) %>%
      pivot_wider(names_from = Sample, values_from = Mutations)

    data_input <- data_input1 %>% left_join(data_input2)

    indcolors <- c(rep("Input", length(colnames(data_input1)[-1])), rep(cancerType, length(colnames(data_input2)[-1])))
    names(indcolors) <- c(colnames(data_input1)[-1], colnames(data_input2)[-1])

    # PCA plot ----------------------------------------------------------------
    mdata_input <- as.matrix(data_input[, -1])
    rownames(mdata_input) <- data_input$MutationType

    res.pca <- prcomp(t(mdata_input), scale = FALSE, center = FALSE)
    xleng <- dim(res.pca$x)[2] * 0.25 + 1
    xleng <- if_else(xleng > 4, 4, xleng)
    yleng <- 4
    pcap1 <- fviz_eig(res.pca, ncp = 10, main = "", addlabels = T)
    ggsave(filename = pca1, plot = pcap1, width = xleng, height = yleng)


    if ((dim(data_input)[2] - 1) > 35) {
      pcap2 <- fviz_pca_ind(res.pca, axes = c(1, 2), geom = "point", col.ind = indcolors, repel = TRUE) + scale_color_brewer(palette = "Set1")
    } else {
      pcap2 <- fviz_pca_ind(res.pca, axes = c(1, 2), col.ind = indcolors, repel = TRUE) + scale_color_brewer(palette = "Set1")
    }
    ggsave(filename = pca2, plot = pcap2, width = 10, height = 7)

    ## a link to download pca data1
    res.pca$x %>% as.data.frame() %>% rownames_to_column(var = 'Sample') %>% write_delim(pca2Data, delim = '\t', col_names = T)

    pcap3 <- fviz_pca_var(res.pca, axes = c(1, 2),
                        col.var = "contrib", # Color by contributions to the PC
                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                        repel = FALSE # Avoid text overlapping
  )
    ggsave(filename = pca3, plot = pcap3, width = 10, height = 7)

    # a link to download pca data2
    res.pca$rotation %>% as.data.frame() %>% rownames_to_column(var = 'Sample') %>% write_delim(pca3Data, delim = '\t', col_names = T)


    output = list(
      'pca1' = pca1,
      'pca2' = pca2,
      'pca3' = pca3,
      'pca2Data' = pca2Data,
      'pca3Data' = pca3Data
     )
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

# Kataegis Identificaiton -----------------------------------------------------------
#parameters for this function: Sample Name, Highlight Kataegis Mutations, Minimum Number of Mutations, Maximum Distance, Chromosome
# Sample Name: same as other module
# Highlight Kataegis Mutations: logical, default FALSE
# Minimum Number of Mutations: default 5
# Maximum Distance: 1000
# Chromosome: default null; can choose from chr1:chr22, chrX and chrY
kataegis <- function(sample, highlight, min, max, chromosome, projectID, pythonOutput, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    errors = ''
    kataegis = paste0(savePath, 'kataegis.svg')
    kataegisData = ''
    kataegisJSON = ''

    mutationPath = paste0(savePath, '../')
    mutationFile = list.files(mutationPath, pattern = "_mSigPortal_SNV")
    mutation_file = paste0(mutationPath, mutationFile)

    if (chromosome == 'None') chromosome = NULL

    if (file.exists(mutation_file)) {
      mutation_data <- read_delim(file = mutation_file, delim = '\t', col_names = FALSE)
      colnames(mutation_data) <- c('project', 'sample', 'type', 'genome_build', 'mutation_type', 'chr', 'pos', 'end', 'ref', 'alt', 'source')

      genome_build <- mutation_data$genome_build[1]
      mutdata <- mutation_data %>% filter(sample == sample) %>% dplyr::select(chr, pos, ref, alt)

      kataegis_result <- kataegis_rainfall_plot(mutdata, sample_name = sample, genome_build = genome_build, reference_data_folder = paste0(dataPath, 'Others'), chromsome = chromosome, kataegis_highligh = highlight, min.mut = min, max.dis = max, filename = kataegis)

      if (is.data.frame(kataegis_result)) {
        kataegisData = paste0(savePath, 'kataegisData.txt')
        kataegisJSON = toJSON(kataegis_result, auto_unbox = TRUE)
        kataegis_result %>% write_delim(kataegisData, delim = '\t', col_names = T)
      }
    }
    else {
      stop("Mutation File Unavailable")
    }

    output = list(
      'plotPath' = kataegis,
      'txtPath' = kataegisData
     )
  }, error = function(e) {
    errors <<- e$message
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output, 'errors' = errors, 'data' = kataegisJSON), pretty = TRUE, auto_unbox = TRUE))
  })
}