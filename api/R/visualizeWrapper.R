library(tidyverse)
library(cowplot)
library(ggsci)
library(ggrepel)
library(hrbrthemes)

import_roboto_condensed()
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
generatePlots <- function(profileType, matrixSize, signatureSetName, sampleName, signatureName, formula, projectID, pythonOutput, savePath) {
  tryCatch({
    stdout <- vector('character')
    con <- textConnection('stdout', 'wr', local = TRUE)
    sink(con, type = "message")
    sink(con, type = "output")


    # section 1: Cosine similarity within samples #
    # two parameters need: Profile Type and Matrix Size #
    if (profileType == "SBS") {
      data_input <- read_delim(paste0(pythonOutput, '/SBS/', projectID, '.', matrixSize, '.all'), delim = '\t')
    }
    else if (profileType == "DBS") {
      data_input <- read_delim(paste0(pythonOutput, '/DBS/', projectID, '.', matrixSize, '.all'), delim = '\t')
    }
    else if (profileType == "ID") {
      data_input <- read_delim(paste0(pythonOutput, '/ID/', projectID, '.', matrixSize, '.all'), delim = '\t')
    }

    ## remove the false data from the result. Don't know why the web generated two additonal columns: False and seen, i think you trying to find the bug and remove.
    data_input <- data_input %>% select(-False, - seen)
    # Heatmap of cosine similarity within samples  and put on the web---------------------------
    cos_sim_res1 = cos_sim_df(data_input, data_input)
    plot_cosine_heatmap_df(cos_sim_res1, cluster_rows = TRUE, plot_values = FALSE)
    ggsave(paste0(savePath, 'cosineSimilarityBetweenSamples.svg'))
    # a link to download the cosine similarity bellow the plot
    # you could rename the file name if you need
    cos_sim_res1 %>% write_delim('cos_sim_res1.txt', delim = '\t', col_names = T)


    # section 2: Cosine similarity  to reference signatures #
    # Two parameters need: Profile Type, Reference Signature Set
    # Profile Type only support SBS, DBS, ID


    # on the web, dropdown list for the input sample and (specific signature or a formula), put the barplot on the web
    # sample_input <- sampleName

    if (!is.null(signatureName)) {
      # option1: specific signature
      profile1 <- data_input %>% select(MutationType, one_of(sampleName))
      profile2 <- sigref %>% select(MutationType, one_of(signatureName))
      # return(c(profile1, profile2))
      plot_compare_profiles_diff(profile1, profile2, condensed = FALSE)
      ggsave(paste0(savePath, 'barplot1.svg'))
    } else if (!is.null(formula)) {
      # option2: formula input
      formula_input <- "0.8*SBS5;0.1*SBS1;0.1*SBS3"
      profile1 <- data_input %>% select(MutationType, one_of(sampleName))
      profile2 <- signature_sum_operation(sigdatabase = signature_refsets_input, sigsetname = signatureSetName, formulax = formula_input, outsigname = "Reconstructed")
      plot_compare_profiles_diff(profile1, profile2, condensed = FALSE, profile_names = c("Original", "Reconstructed"))
      ggsave(paste0(savePath, 'barplot2.svg'))
    }

    # PCA plot ----------------------------------------------------------------

    data_input
    mdata_input <- as.matrix(data_input[, -1])
    rownames(mdata_input) <- data_input$MutationType

    res.pca <- prcomp(t(mdata_input), scale = FALSE, center = FALSE)
    fviz_eig(res.pca, ncp = 10)
    ggsave(paste0(savePath, 'pca1.svg'))

    if ((dim(data_input)[2] - 1) > 20) {
      fviz_pca_ind(res.pca, axes = c(1, 2), geom = "point", col.ind = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
      ggsave(paste0(savePath, 'pca2.svg'))
    } else {
      fviz_pca_ind(res.pca, axes = c(1, 2), col.ind = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
      ggsave(paste0(savePath, 'pca2.svg'))
    }


    fviz_pca_var(res.pca, axes = c(1, 2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE # Avoid text overlapping
  )
    ggsave(paste0(savePath, 'pca3.svg'))

    ## output PCAs
    res.pca$x
    res.pca$rotation

    ## heatmap between PCs and signatures
    sigpca <- res.pca$rotation %>% as.data.frame() %>% rownames_to_column(var = "MutationType")
    cos_sim_res = cos_sim_df(sigpca, sigref)

    # put this heatmap on the web
    plot_cosine_heatmap_df(cos_sim_res, cluster_rows = TRUE, plot_values = FALSE)
    ggsave(paste0(savePath, 'heatmap2.svg'))
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(stdout)
  })
  # t-SNE plots, coming soon ###
}

