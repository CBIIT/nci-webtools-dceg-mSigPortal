library(tidyverse)
library(cowplot)
library(hrbrthemes)
library(ggsci)
library(ggrepel)
library(factoextra)


# source msigportal function ----------------------------------------------
source('api/R/Sigvisualfunc.R')

# load reference signatures files -----------------------------------------
load('api/R/signature_refsets.RData')

generatePlots <- function(profileName, signatureSetName, sampleName, signatureName, formula, projectID, savePath) {
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con)
  tryCatch({

    # on the web, we need the dropdown list for signature Profile Type and reference set name --------
    matrix <- if_else(profileName == "SBS", "SBS96", if_else(profileName == "DBS", "DBS78", if_else(profileName == "ID", "ID83", NA_character_)))

    signature_refsets_input <- signature_refsets %>% filter(Profile == matrix, Signature_set_name == signatureSetName)


    # load SBS-96/DBS-78/ID-83 catalog file generated from previous result according to the profileName -----------
    if (profileName == "SBS") {
      data_input <- read_delim(paste0(savePath, '/SBS/', projectID, '.', matrix, '.all'), delim = '\t')
    }
    else if (profileName == "DBS") {
      data_input <- read_delim(paste0(savePath, '/DBS/', projectID, '.', matrix, '.all'), delim = '\t')
    }
    else if (profileName == "ID") {
      data_input <- read_delim(paste0(savePath, '/ID/', projectID, '.', matrix, '.all'), delim = '\t')
    }

    ## remove the false data from the result. Don't know why the web generated two additonal columns: False and seen, i think you trying to find the bug and remove.
    # data_input <- data_input %>% select(-seen)

    # Heatmap of cosine similarity to reference set ---------------------------

    sigref <- signature_refsets_input %>%
  select(Signature_name, MutationType, Contribution) %>%
  pivot_wider(names_from = Signature_name, values_from = Contribution)

    cos_sim_res = cos_sim_df(data_input, sigref)
    # return(cos_sim_res)
    # put this heatmap on the web
    plot_cosine_heatmap_df(cos_sim_res, cluster_rows = TRUE, plot_values = FALSE)
    ggsave('heatmap1.svg')


    # on the web, dropdown list for the input sample and (specific signature or a formula), put the barplot on the web
    # sample_input <- sampleName

    if (!is.null(signatureName)) {
      # option1: specific signature
      profile1 <- data_input %>% select(MutationType, one_of(sampleName))
      profile2 <- sigref %>% select(MutationType, one_of(signatureName))
      # return(c(profile1, profile2))
      plot_compare_profiles_diff(profile1, profile2, condensed = FALSE)
      ggsave('barplot1.svg')
    } else if (!is.null(formula)) {
      # option2: formula input
      formula_input <- "0.8*SBS5;0.1*SBS1;0.1*SBS3"
      profile1 <- data_input %>% select(MutationType, one_of(sampleName))
      profile2 <- signature_sum_operation(sigdatabase = signature_refsets_input, sigsetname = signatureSetName, formulax = formula_input, outsigname = "Reconstructed")
      plot_compare_profiles_diff(profile1, profile2, condensed = FALSE, profile_names = c("Original", "Reconstructed"))
      ggsave('barplot2.svg')
    }

    # PCA plot ----------------------------------------------------------------

    data_input
    mdata_input <- as.matrix(data_input[, -1])
    rownames(mdata_input) <- data_input$MutationType

    res.pca <- prcomp(t(mdata_input), scale = FALSE, center = FALSE)
    fviz_eig(res.pca, ncp = 10)
    ggsave('pca1.svg')

    if ((dim(data_input)[2] - 1) > 20) {
      fviz_pca_ind(res.pca, axes = c(1, 2), geom = "point", col.ind = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
      ggsave('pca2.svg')
    } else {
      fviz_pca_ind(res.pca, axes = c(1, 2), col.ind = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
      ggsave('pca2.svg')
    }


    fviz_pca_var(res.pca, axes = c(1, 2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE # Avoid text overlapping
  )
    ggsave('pca3.svg')

    ## output PCAs
    res.pca$x
    res.pca$rotation

    ## heatmap between PCs and signatures
    sigpca <- res.pca$rotation %>% as.data.frame() %>% rownames_to_column(var = "MutationType")
    cos_sim_res = cos_sim_df(sigpca, sigref)

    # put this heatmap on the web
    plot_cosine_heatmap_df(cos_sim_res, cluster_rows = TRUE, plot_values = FALSE)
    ggsave('heatmap2.svg')

  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    return(stdout)
  })
  # t-SNE plots, coming soon ###
}






