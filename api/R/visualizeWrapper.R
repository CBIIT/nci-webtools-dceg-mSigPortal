library(tidyverse)
library(hrbrthemes)
library(cowplot)
library(ggsci)
library(ggrepel)
library(factoextra)
library(jsonlite)


# source msigportal function ----------------------------------------------
source('api/R/Sigvisualfunc.R')
# load reference signatures files -----------------------------------------
load('api/R/signature_refsets.RData')

# get all Reference Signature Set options using profileName (matrix)
getSignatureReferenceSets <- function(profileType) {
  profileName <- if_else(profileType == "SBS", "SBS96", if_else(profileType == "DBS", "DBS78", if_else(profileType == "ID", "ID83", NA_character_)))
  return(signature_refsets %>% filter(Profile == profileName) %>% pull(Signature_set_name) %>% unique())
}

### Cosine Similarity tab ###
# section 1: Cosine similarity within samples #
# two parameters need: Profile Type and Matrix Size #
cosineSimilarityWithin <- function(profileType, matrixSize, projectID, pythonOutput, savePath) {
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    paths = list()
    plotPath = paste0(savePath, '/cos_sim_within.svg')
    txtPath = paste0(savePath, '/cos_sim_within.txt')
    data_input <- read_delim(paste0(pythonOutput, '/', profileType, '/', projectID, '.', matrixSize, '.all'), delim = '\t')
    if ("False" %in% colnames(data_input)) {
      data_input <- data_input %>% select(-False)
    }
    # Heatmap of cosine similarity within samples  and put on the web---------------------------
    cos_sim_res1 = cos_sim_df(data_input, data_input)
    plot_cosine_heatmap_df(cos_sim_res1, cluster_rows = TRUE, plot_values = FALSE)
    ggsave(plotPath)
    cos_sim_res1 %>% write_delim(txtPath, delim = '\t', col_names = T)

    paths = list('plotPath' = plotPath[1], 'txtPath' = txtPath[1])
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'paths' = paths), pretty = TRUE, auto_unbox = TRUE))
  })
}

# section 2: Cosine similarity  to reference signatures 
# Two parameters need: Profile Type, Reference Signature Set
# Profile Type only support SBS, DBS, ID
cosineSimilarityRefSig <- function(profileType, signatureSetName, projectID, pythonOutput, savePath) {
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    paths = list()
    plotPath = paste0(savePath, 'cos_sim_refsig.svg')
    txtPath = paste0(savePath, 'cos_sim_refsig.txt')
    matrix <- if_else(profileType == "SBS", "SBS96", if_else(profileType == "DBS", "DBS78", if_else(profileType == "ID", "ID83", NA_character_)))
    signature_refsets_input <- signature_refsets %>% filter(Profile == matrix, Signature_set_name == signatureSetName)
    refsig <- signature_refsets_input %>%
      select(Signature_name, MutationType, Contribution) %>%
      pivot_wider(names_from = Signature_name, values_from = Contribution)

    data_input <- read_delim(paste0(pythonOutput, '/', profileType, '/', projectID, '.', matrix, '.all'), delim = '\t')
    if ("False" %in% colnames(data_input)) {
      data_input <- data_input %>% select(-False)
    }
    # Heatmap of cosine similarity to reference set signature and put on the web---------------------------
    cos_sim_res2 = cos_sim_df(data_input, refsig)
    plot_cosine_heatmap_df(cos_sim_res2, cluster_rows = TRUE, plot_values = FALSE)
    ggsave(plotPath)
    cos_sim_res2 %>% write_delim(txtPath, delim = '\t', col_names = T)

    paths = list('plotPath' = plotPath, 'txtPath' = txtPath)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'paths' = paths), pretty = TRUE, auto_unbox = TRUE))
  })
}

### Profile Comparison tab ###
# section 1: Cosine similarity within samples 
# three parameters need: Profile Type, Sample Name1 and Sample Name2
profileComparisonWithin <- function(profileType, sampleName1, sampleName2, projectID, pythonOutput, savePath) {
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    paths = list()
    plotPath = paste0(savePath, '/cos_sim_within.svg')
    matrixsize <- if_else(profileType == "SBS", "SBS96", if_else(profileType == "DBS", "DBS78", if_else(profileType == "ID", "ID83", NA_character_)))
    data_input <- read_delim(paste0(pythonOutput, '/', profileType, '/', projectID, '.', matrixsize, '.all'), delim = '\t')
    # data_input <- data_input %>% select(-False,-seen)
    profile1 <- data_input %>% select(MutationType, contains(sampleName1))
    profile2 <- data_input %>% select(MutationType, contains(sampleName2))
    plot_compare_profiles_diff(profile1, profile2, condensed = FALSE)
    ggsave(plotPath)

    paths = list('plotPath' = plotPath)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'paths' = paths), pretty = TRUE, auto_unbox = TRUE))
  })
}

# section 2: Comparison to reference signatures # 
# four parameters need: “Profile Type”, “Sample Name”, “Reference Signature Set” and “Compare Single Signature or Combined Signatures” # 
profileComparisonRefSig <- function(profileType, sampleName, signatureSetName, compare, projectID, pythonOutput, savePath) {
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    paths = list()
    plotPath = paste0(savePath, 'cos_sim_refsig.svg')
    profilename <- if_else(profileType == "SBS", "SBS96", if_else(profileType == "DBS", "DBS78", if_else(profileType == "ID", "ID83", NA_character_)))
    matrixsize <- profilename
    data_input <- read_delim(paste0(pythonOutput, '/', profileType, '/', projectID, '.', matrixsize, '.all'), delim = '\t')
    profile1 <- data_input %>% select(MutationType, one_of(profileType))
    signature_refsets_input <- signature_refsets %>% filter(Profile == profilename, Signature_set_name == signatureSetName)
    refsig <- signature_refsets_input %>%
      select(Signature_name, MutationType, Contribution) %>%
      pivot_wider(names_from = Signature_name, values_from = Contribution)

    user_input = compare # user_input <- "SBS5" #user_input <- "0.8*SBS5;0.1*SBS1;0.1*SBS3"

    if (str_detect(user_input, ";")) {
      profile2 <- signature_sum_operation(sigdatabase = signature_refsets_input, sigsetname = signatureSetName, formulax = user_input, outsigname = "Reconstructed")
      profile_names = c("Original", "Reconstructed")
    } else {
      profile2 <- refsig %>% select(MutationType, one_of(user_input))
      profile_names = c(colnames(profile1)[2], colnames(profile2)[2])
    }

    plot_compare_profiles_diff(profile1, profile2, condensed = FALSE, profile_names = profile_names)
    ggsave(plotPath)

    paths = list('plotPath' = plotPath)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'paths' = paths), pretty = TRUE, auto_unbox = TRUE))
  })

}