library(tidyverse)
library(cowplot)
library(hrbrthemes)
library(ggsci)
library(ggrepel)
library(ggforce)
# library(ggtext) requires R 3.5
library(jsonlite)




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
mutationalSigPro <- function(projectID, pythonOutput, savePath, dataPath) {
  source('services/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'signature_refsets.RData'))
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()

    path_profile <- paste0(dataPath, 'Reference_Signature_Profiles_SVG/')
    signature_profile_files <- signature_refsets %>% select(Source, Profile, Signature_set_name, Dataset, Signature_name) %>% unique() %>% mutate(Path = str_replace_all(Signature_set_name, " ", "_"), Path = str_remove_all(Path, "[()]"), Path = paste0(path_profile, Path, "/", Signature_name, ".svg"))
    # print(signature_profile_files$Path)
    ## multiple selected profiles
    signature_source_input <- "Reference_signatures"
    profile_name_input <- "SBS96"
    signatureset_name_input <- "COSMIC v3 Signatures (SBS)"
    dataset_input <- "WGS"
    signature_name_input <- "SBS1"
    svgfile_selected <- signature_profile_files %>% filter(Source == signature_source_input, Profile == profile_name_input, Signature_set_name == signatureset_name_input, Dataset == dataset_input, Signature_name == signature_name_input) %>% pull(Path)
    # print(svgfile_selected)

    output = list('plotPath' = svgfile_selected)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

