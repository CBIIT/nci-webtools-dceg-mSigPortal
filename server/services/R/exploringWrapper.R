library(tidyverse)
library(cowplot)
library(hrbrthemes)
library(ggsci)
library(ggrepel)
library(ggforce)
library(ggtext)
library(jsonlite)
library(stringr)




# Util Functions
# get dataframe with column and filter arguments
getReferenceSignatureData <- function(args, dataPath) {
  load(paste0(dataPath, 'Signature/signature_refsets.RData'))
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

# Signature Explore -------------------------------------------------------
# section 1: Current reference signatures in mSigPortal -------------------
referenceSignatures <- function(projectID, pythonOutput, rootDir, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'Signature/signature_refsets.RData'))
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
mutationalProfiles <- function(signatureSource, profileName, refSignatureSet, experimentalStrategy, signatureName, projectID, pythonOutput, rootDir, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'Signature/signature_refsets.RData'))
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()

    path_profile <- paste0(dataPath, 'Signature/Reference_Signature_Profiles_SVG/')
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
cosineSimilarity <- function(profileName, refSignatureSet1, refSignatureSet2, projectID, pythonOutput, rootDir, savePath, dataPath) {
  # The parameters will be “Matrix Size”, “Reference Signature Set1” and “Reference Signature Set2”. 
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'Signature/signature_refsets.RData'))
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
mutationalSignatureComparison <- function(profileName, refSignatureSet1, signatureName1, refSignatureSet2, signatureName2, projectID, pythonOutput, rootDir, savePath, dataPath) {
  # The parameters will be “Matrix Size”, “Reference Signature Set1” and “Reference Signature Set2”. 
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'Signature/signature_refsets.RData'))
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
  if (useCancer == TRUE)
    signature_association(data = data_input, cancer_type_input = cancerType, signature_name_input1 = signatureName1, signature_name_input2 = signatureName2, signature_both = both, output_plot = plotPath)
  else
    signature_association(data = data_input, cancer_type_input = NULL, signature_name_input1 = signatureName1, signature_name_input2 = signatureName2, signature_both = both, output_plot = plotPath)
}

mutationalSignatureDecomposition <- function(plotPath, dataPath, exposure_refdata, signature_refsets, seqmatrix_refdata) {
  exposure_refdata_input <- exposure_refdata %>% mutate(Sample = paste0(Cancer_Type, "@", Sample)) %>%
      select(Sample, Signature_name, Exposure) %>%
      pivot_wider(id_cols = Sample, names_from = Signature_name, values_from = Exposure)

  signature_refsets_input <- signature_refsets %>%
      select(MutationType, Signature_name, Contribution) %>%
      pivot_wider(id_cols = MutationType, names_from = Signature_name, values_from = Contribution) %>%
      arrange(MutationType) # have to sort the mutationtype

  seqmatrix_refdata_input <- seqmatrix_refdata %>% mutate(Sample = paste0(Cancer_Type, "@", Sample)) %>%
      select(MutationType, Sample, Mutations) %>%
      pivot_wider(id_cols = MutationType, names_from = Sample, values_from = Mutations) %>%
      arrange(MutationType) ## have to sort the mutation type


  decompsite_input <- calculate_similarities(orignal_genomes = seqmatrix_refdata_input, signature = signature_refsets_input, signature_activaties = exposure_refdata_input)
  decompsite_input <- decompsite_input %>% separate(col = Sample_Names, into = c('Cancer_Type', 'Sample'), sep = '@')
  decompsite_input %>% write_delim(dataPath, delim = '\t', col_names = T) ## put the link to download this table

  decompsite_distribution(decompsite = decompsite_input, output_plot = plotPath) # put the distribution plot online.
}

mutationalSignatureLandscape <- function(cancerType, varDataPath, plotPath, exposure_refdata, signature_refsets, seqmatrix_refdata) {
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

  cosinedata <- decompsite_input %>% select(Samples = Sample_Names, Similarity = Cosine_similarity)

  data_input <- exposure_refdata %>%
      filter(Cancer_Type == cancerType) %>%
      select(Sample, Signature_name, Exposure) %>%
      pivot_wider(id_cols = Sample, names_from = Signature_name, values_from = Exposure) %>%
      rename(Samples = Sample)

  data_input <- data_input %>% select_if(~!is.numeric(.) || sum(.) > 0)

  sigdata <- data_input
  ## two parameters to add the two bars: vardata1, vardata1_cat, vardata2, vardata2_cat
  # studydata <- data_input %>% select(Samples) %>% mutate(Study=if_else((seq_along(Samples) %% 2 ==0), "A","B"))
  # puritydata <-  data_input %>% select(Samples) %>% mutate(Purity=0)
  # puritydata$Purity <- runif(n=length(puritydata$Purity), min=1e-12, max=.9999999999)
  # highlight <-  c('SP124389','SP124273')
  # Exposure_Clustering(sigdata = sigdata,studydata = studydata,studyname = "VAR1",puritydata = puritydata,purityname = 'VAR2',cosinedata = cosinedata,clustern=5,output_plot = 'tmp.svg' )
  # Exposure_Clustering(sigdata = sigdata,cosinedata = cosinedata,clustern=5,output_plot = 'tmp.svg' )

  # parameter: Cancer Type, Vardata_input_file
  if (stringi::stri_length(varDataPath) > 0) {
    vardata_input <- read_delim(varDataPath, delim = '\t', col_names = T)

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

  prevalence_plot(sigdata = sigdata, nmutation = mutation, output_plot = plotPath)
}

exposurePublic <- function(fn, common, across = '{}', association = '{}', landscape = '{}', prevalence = '{}', projectID, pythonOutput, rootDir, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'Signature/signature_refsets.RData'))
  load(paste0(dataPath, 'Seqmatrix/seqmatrix_refdata.RData'))
  load(paste0(dataPath, 'Exposure/exposure_refdata.RData'))
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    tumorPath = paste0(savePath, 'tumorMutationalBurden.svg')
    burdenSeparatedPath = paste0(savePath, 'burdenSeparatedPath.svg')
    burdenAcrossPath = paste0(savePath, 'burdenAcrossCancer.svg')
    associationPath = paste0(savePath, 'mutationalSignatureAssociation.svg')
    decompositionPath = paste0(savePath, 'mutationalSignatureDecomposition.svg')
    decompositionData = paste0(savePath, 'mutationalSignatureDecomposition.txt')
    landscapePath = paste0(savePath, 'landscapeMutationalSignature.svg')
    prevalencePath = paste0(savePath, 'prevalenceMutationalSignature.svg')

    # parse arguments
    common = fromJSON(common)
    across = fromJSON(across)
    association = fromJSON(association)
    landscape = fromJSON(landscape)
    prevalence = fromJSON(prevalence)

    exposure_refdata_selected <- exposure_refdata %>% filter(Study == common$study, Dataset == common$strategy, Signature_set_name == common$refSignatureSet)

    genome <- case_when(
    common$study == "PCAWG" ~ "GRCh37",
    common$study == "TCGA" ~ "GRCh37",
    TRUE ~ NA_character_
  )
    genomesize = genome2size(genome)

    signature_refsets_selected <- signature_refsets %>%
      filter(Signature_set_name == common$refSignatureSet)

    seqmatrix_refdata_selected <- seqmatrix_refdata %>% filter(Study == common$study, Dataset == common$strategy, Profile == signature_refsets_selected$Profile[1])


    ## Tumor Overall Mutational Burden
    if ('all' %in% fn) {
      print('Tumor Mutational Burden')
      tumorMutationalBurden(genomesize, tumorPath, exposure_refdata_selected)
      output[['tumorPath']] = tumorPath
    }

    # Tumor Mutational Burden separated by signatures
    if ('all' %in% fn || 'separated' %in% fn) {
      print('Tumor Mutational Burden Separated by Signatures')
      mutationalSignatureBurdenSeparated(genomesize, common$cancerType, burdenSeparatedPath, exposure_refdata_selected)
      output[['burdenSeparatedPath']] = burdenSeparatedPath
    }

    # Mutational signature burden across cancer types
    if ('all' %in% fn || 'across' %in% fn) {
      print('Mutational Signature Burden Across Cancer Types')
      mutationalSignatureBurdenAcrossCancer(across$signatureName, genomesize, burdenAcrossPath, exposure_refdata_selected)
      output[['burdenAcrossPath']] = burdenAcrossPath
    }
    # Mutational Signature Association
    if ('all' %in% fn | 'association' %in% fn) {
      print('Evaluating the Performance of Mutational Signature Decomposition')
      mutationalSignatureAssociation(association$useCancer, common$cancerType, association$both, association$signatureName1, association$signatureName2, associationPath, exposure_refdata_selected)
      output[['associationPath']] = associationPath
    }
    # Evaluating the Performance of Mutational Signature Decomposition --------
    if ('all' %in% fn) {
      print('Mutational Signature Association')
      mutationalSignatureDecomposition(decompositionPath, decompositionData, exposure_refdata_selected, signature_refsets_selected, seqmatrix_refdata_selected)
      output[['decompositionPath']] = decompositionPath
    }
    # Landscape of Mutational Signature Activity
    if ('all' %in% fn | 'landscape' %in% fn) {
      print('Landscape of Mutational Signature Activity')
      varDataPath = ''
      if (stringi::stri_length(landscape$variableFile) > 0) {
        varDataPath = file.path(rootDir, landscape$variableFile)
      }
      mutationalSignatureLandscape(common$cancerType, varDataPath, landscapePath, exposure_refdata_selected, signature_refsets_selected, seqmatrix_refdata_selected)
      output[['landscapePath']] = landscapePath
    }
    # Prevalence plot
    if ('all' %in% fn | 'prevalence' %in% fn) {
      print('Prevalence of Mutational Signature')
      mutationalSignaturePrevalence(prevalence$mutation, common$cancerType, prevalencePath, exposure_refdata_selected)
      output[['prevalencePath']] = prevalencePath
    }

  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

exposureUser <- function(fn, files, common, across = '{}', association = '{}', landscape = '{}', prevalence = '{}', projectID, pythonOutput, rootDir, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')

  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    tumorPath = paste0(savePath, 'tumorMutationalBurden.svg')
    burdenSeparatedPath = paste0(savePath, 'burdenSeparatedPath.svg')
    burdenAcrossPath = paste0(savePath, 'burdenAcrossCancer.svg')
    associationPath = paste0(savePath, 'mutationalSignatureAssociation.svg')
    decompositionPath = paste0(savePath, 'mutationalSignatureDecomposition.svg')
    decompositionData = paste0(savePath, 'mutationalSignatureDecomposition.txt')
    landscapePath = paste0(savePath, 'landscapeMutationalSignature.svg')
    prevalencePath = paste0(savePath, 'prevalenceMutationalSignature.svg')

    # parse arguments
    common = fromJSON(common)
    across = fromJSON(across)
    association = fromJSON(association)
    landscape = fromJSON(landscape)
    prevalence = fromJSON(prevalence)
    files = fromJSON(files)

    exposure_refdata_selected <- read_delim(file.path(rootDir, files$exposureFile), delim = '\t', col_names = T)
    seqmatrix_refdata_selected <- read_delim(file.path(rootDir, files$matrixFile), delim = '\t', col_names = T)

    if (stringi::stri_length(files$signatureFile) > 0) {
      # if using user uploaded signature file
      signature_refsets_selected <- read_delim(file.path(rootDir, files$signatureFile), delim = '\t', col_names = T)
    } else {
      # else use public signature data
      load(paste0(dataPath, 'Signature/signature_refsets.RData'))
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
    if ('all' %in% fn) {
      print('Tumor Mutational Burden')
      tumorMutationalBurden(genomesize, tumorPath, exposure_refdata_selected)
      output[['tumorPath']] = tumorPath
    }

    # Tumor Mutational Burden separated by signatures
    if ('all' %in% fn || 'separated' %in% fn) {
      print('Tumor Mutational Burden Separated by Signatures')
      mutationalSignatureBurdenSeparated(genomesize, cancer_type_user, burdenSeparatedPath, exposure_refdata_selected)
      output[['burdenSeparatedPath']] = burdenSeparatedPath
    }
    # Mutational signature burden across cancer types
    if ('all' %in% fn || 'across' %in% fn) {
      print('Mutational Signature Burden Across Cancer Types')
      mutationalSignatureBurdenAcrossCancer(across$signatureName, genomesize, burdenAcrossPath, exposure_refdata_selected)
      output[['burdenAcrossPath']] = burdenAcrossPath
    }
    # Mutational Signature Association
    if ('all' %in% fn | 'association' %in% fn) {
      print('Evaluating the Performance of Mutational Signature Decomposition')
      mutationalSignatureAssociation(association$useCancer, cancer_type_user, association$both, association$signatureName1, association$signatureName2, associationPath, exposure_refdata_selected)
      output[['associationPath']] = associationPath
    }
    # Evaluating the Performance of Mutational Signature Decomposition --------
    if ('all' %in% fn) {
      print('Mutational Signature Association')
      mutationalSignatureDecomposition(decompositionPath, decompositionData, exposure_refdata_selected, signature_refsets_selected, seqmatrix_refdata_selected)
      output[['decompositionPath']] = decompositionPath
    }
    # Landscape of Mutational Signature Activity
    if ('all' %in% fn | 'landscape' %in% fn) {
      print('Landscape of Mutational Signature Activity')
      varDataPath = ''
      if (stringi::stri_length(landscape$variableFile) > 0) {
        varDataPath = file.path(rootDir, landscape$variableFile$signatureFile)
      }
      mutationalSignatureLandscape(cancer_type_user, varDataPath, landscapePath, exposure_refdata_selected, signature_refsets_selected, seqmatrix_refdata_selected)
      output[['landscapePath']] = landscapePath
    }
    # Prevalence plot
    if ('all' %in% fn | 'prevalence' %in% fn) {
      print('Prevalence of Mutational Signature')
      mutationalSignaturePrevalence(prevalence$mutation, cancer_type_user, prevalencePath, exposure_refdata_selected)
      output[['prevalencePath']] = prevalencePath
    }

  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}
