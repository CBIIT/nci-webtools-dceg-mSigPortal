library(tidyverse)
library(hrbrthemes)
import_roboto_condensed()
library(cowplot)
library(ggsci)
library(ggrepel)
library(factoextra)


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
cosineSimilarity <- function(profileType1, matrixSize, profileType2, signatureSetName, projectID, pythonOutput, savePath) {
  tryCatch({
    stdout <- vector('character')
    con <- textConnection('stdout', 'wr', local = TRUE)
    sink(con, type = "message")
    sink(con, type = "output")

    # section 1: Cosine similarity within samples #
    # two parameters need: Profile Type and Matrix Size #
    data_input <- read_delim(paste0(pythonOutput, '/', profileType1, '/', projectID, '.', matrixSize, '.all'), delim = '\t')

    # data_input <- data_input %>% select(-False, - seen)
    # Heatmap of cosine similarity within samples  and put on the web---------------------------
    cos_sim_res1 = cos_sim_df(data_input, data_input)
    plot_cosine_heatmap_df(cos_sim_res1, cluster_rows = TRUE, plot_values = FALSE)
    ggsave(paste0(savePath, 'cos_sim_within.svg'))
    cos_sim_res1 %>% write_delim(paste0(savePath, 'cos_sim_within.txt'), delim = '\t', col_names = T)


    # section 2: Cosine similarity  to reference signatures #
    # Two parameters need: Profile Type, Reference Signature Set
    # Profile Type only support SBS, DBS, ID
    matrixSize2 <- if_else(profileType2 == "SBS", "SBS96", if_else(profileType2 == "DBS", "DBS78", if_else(profileType2 == "ID", "ID83", NA_character_)))
    signature_refsets_input <- signature_refsets %>% filter(Profile == matrixSize2, Signature_set_name == signatureSetName)
    sigref <- signature_refsets_input %>%
      select(Signature_name, MutationType, Contribution) %>%
      pivot_wider(names_from = Signature_name, values_from = Contribution)

    data_input <- read_delim(paste0(pythonOutput, '/', profileType2, '/', projectID, '.', matrixSize2, '.all'), delim = '\t')
    # data_input <- data_input %>% select(-False, - seen)

    # Heatmap of cosine similarity to reference set signature and put on the web---------------------------
    cos_sim_res2 = cos_sim_df(data_input, sigref)
    plot_cosine_heatmap_df(cos_sim_res2, cluster_rows = TRUE, plot_values = FALSE)
    ggsave(paste0(savePath, 'cos_sim_sigref.svg'))
    cos_sim_res2 %>% write_delim(paste0(savePath, 'cos_sim_sigref.txt', delim = '\t', col_names = T))

  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(stdout)
  })
}
