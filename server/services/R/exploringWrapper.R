library(tidyverse)
library(cowplot)
library(hrbrthemes)
library(ggsci)
library(ggrepel)
library(ggforce)
# library(ggtext) requires R 3.5
library(jsonlite)
library(stringr)




# Util Functions
# get dataframe with column and filter arguments
getReferenceSignatureData <- function(args, dataPath) {
  load(paste0(dataPath, 'signature_refsets.RData'))
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
referenceSignatures <- function(projectID, pythonOutput, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'signature_refsets.RData'))
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
mutationalProfiles <- function(signatureSource, profileName, refSignatureSet, experimentalStrategy, signatureName, projectID, pythonOutput, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'signature_refsets.RData'))
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()

    path_profile <- paste0(dataPath, 'Reference_Signature_Profiles_SVG/')
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
cosineSimilarity <- function(profileName, refSignatureSet1, refSignatureSet2, projectID, pythonOutput, savePath, dataPath) {
  # The parameters will be “Matrix Size”, “Reference Signature Set1” and “Reference Signature Set2”. 
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'signature_refsets.RData'))
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