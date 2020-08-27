library(tidyverse)
library(hrbrthemes)
library(cowplot)
library(ggsci)
library(ggrepel)
library(factoextra)
library(jsonlite)

# get all Reference Signature Set options using profileName (matrix)
getReferenceSignatureSets <- function(profileType, dataPath) {
  load(paste0(dataPath, 'signature_refsets.RData'))

  profileName <- if_else(profileType == "SBS", "SBS96", if_else(profileType == "DBS", "DBS78", if_else(profileType == "ID", "ID83", NA_character_)))
  signatureSets <- signature_refsets %>% filter(Profile == profileName) %>% pull(Signature_set_name) %>% unique()

  return(signatureSets)
}

# get list of signatures in the selected signature set
getSignatures <- function(profileType, signatureSetName, dataPath) {
  load(paste0(dataPath, 'signature_refsets.RData'))
  profileName <- if_else(profileType == "SBS", "SBS96", if_else(profileType == "DBS", "DBS78", if_else(profileType == "ID", "ID83", NA_character_)))
  signature_refsets_input <- signature_refsets %>% filter(Profile == profileName, Signature_set_name == signatureSetName)
  refsig <- signature_refsets_input %>%
      select(Signature_name, MutationType, Contribution) %>%
      pivot_wider(names_from = Signature_name, values_from = Contribution)
  signatures <- colnames(refsig[1, -1])

  return(signatures)
}

# get svg list for public data
mutationProfilesPublic <- function(dataPath) {
  load(paste0(dataPath, 'signature_refsets.RData'))
  load(paste0(dataPath, 'seqmatrix_refdata_info.RData'))



}

### Cosine Similarity tab ###
# section 1: Cosine similarity within samples #
# two parameters need: Profile Type and Matrix Size #
cosineSimilarityWithin <- function(profileType, matrixSize, projectID, pythonOutput, savePath, dataPath) {
  source('api/R/Sigvisualfunc.R')
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    plotPath = paste0(savePath, '/cos_sim_within.svg')
    txtPath = paste0(savePath, '/cos_sim_within.txt')

    data_input <- read_delim(file.path(pythonOutput, profileType, paste0(projectID, '.', profileType, matrixSize, '.all')), delim = '\t')
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
cosineSimilarityRefSig <- function(profileType, signatureSetName, projectID, pythonOutput, savePath, dataPath) {
  source('api/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'signature_refsets.RData'))
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    plotPath = paste0(savePath, '/cos_sim_refsig.svg')
    txtPath = paste0(savePath, '/cos_sim_refsig.txt')

    matrix <- if_else(profileType == "SBS", "SBS96", if_else(profileType == "DBS", "DBS78", if_else(profileType == "ID", "ID83", NA_character_)))
    signature_refsets_input <- signature_refsets %>% filter(Profile == matrix, Signature_set_name == signatureSetName)
    refsig <- signature_refsets_input %>%
      select(Signature_name, MutationType, Contribution) %>%
      pivot_wider(names_from = Signature_name, values_from = Contribution)

    data_input <- read_delim(paste0(pythonOutput, '/', profileType, '/', projectID, '.', matrix, '.all'), delim = '\t')

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

### Profile Comparison tab ###
# section 1: Cosine similarity within samples 
# three parameters need: Profile Type, Sample Name1 and Sample Name2
profileComparisonWithin <- function(profileType, sampleName1, sampleName2, projectID, pythonOutput, savePath, dataPath) {
  source('api/R/Sigvisualfunc.R')
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    plotPath = paste0(savePath, '/pro_com_within.svg')

    matrixsize <- if_else(profileType == "SBS", "SBS96", if_else(profileType == "DBS", "DBS78", if_else(profileType == "ID", "ID83", NA_character_)))
    data_input <- read_delim(paste0(pythonOutput, '/', profileType, '/', projectID, '.', matrixsize, '.all'), delim = '\t')
    profile1 <- data_input %>% select(MutationType, contains(sampleName1))
    profile2 <- data_input %>% select(MutationType, contains(sampleName2))
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
profileComparisonRefSig <- function(profileType, sampleName, signatureSetName, compare, projectID, pythonOutput, savePath, dataPath) {
  source('api/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'signature_refsets.RData'))
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    plotPath = paste0(savePath, '/pro_com_refsig.svg')

    profilename <- if_else(profileType == "SBS", "SBS96", if_else(profileType == "DBS", "DBS78", if_else(profileType == "ID", "ID83", NA_character_)))
    matrixsize <- profilename
    data_input <- read_delim(paste0(pythonOutput, '/', profileType, '/', projectID, '.', matrixsize, '.all'), delim = '\t')
    profile1 <- data_input %>% select(MutationType, one_of(sampleName))
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

### “Principal Component Analysis” ###
# Two parameters need: Profile Type, Reference Signature Set
# Profile Type only support SBS, DBS, ID
pca <- function(profileType, signatureSetName, projectID, pythonOutput, savePath, dataPath) {
  source('api/R/Sigvisualfunc.R')
  load(paste0(dataPath, 'signature_refsets.RData'))
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()
    pca1 = paste0(savePath, '/pca1.svg')
    pca2 = paste0(savePath, '/pca2.svg')
    pca3 = paste0(savePath, '/pca3.svg')
    heatmap = paste0(savePath, '/heatmap.svg')
    pca2Data = paste0(savePath, '/pca2_data.txt')
    pca3Data = paste0(savePath, '/pca3_data.txt')
    heatmapData = paste0(savePath, '/heatmap_data.txt')



    profilename <- if_else(profileType == "SBS", "SBS96", if_else(profileType == "DBS", "DBS78", if_else(profileType == "ID", "ID83", NA_character_)))
    matrixsize <- profilename

    signature_refsets_input <- signature_refsets %>% filter(Profile == profilename, Signature_set_name == signatureSetName)
    refsig <- signature_refsets_input %>%
      select(Signature_name, MutationType, Contribution) %>%
      pivot_wider(names_from = Signature_name, values_from = Contribution)

    data_input <- read_delim(paste0(pythonOutput, '/', profileType, '/', projectID, '.', matrixsize, '.all'), delim = '\t')

    # PCA plot ----------------------------------------------------------------
    mdata_input <- as.matrix(data_input[, -1])
    rownames(mdata_input) <- data_input$MutationType

    res.pca <- prcomp(t(mdata_input), scale = FALSE, center = FALSE)

    xleng <- dim(res.pca$x)[2] * 0.5 + 1
    yleng <- 4
    pcap1 <- fviz_eig(res.pca, ncp = 10, main = "")
    ggsave(filename = pca1, plot = pcap1, width = xleng, height = yleng)

    if ((dim(data_input)[2] - 1) > 20) {
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

    ## heatmap between PCs and signatures
    sigpca <- res.pca$rotation %>% as.data.frame() %>% rownames_to_column(var = "MutationType")
    cos_sim_res3 = cos_sim_df(sigpca, refsig)
    # put this heatmap on the web
    plot_cosine_heatmap_df(cos_sim_res3, cluster_rows = TRUE, plot_values = FALSE, output_plot = heatmap)

    # a link to download the cosine similarity bellow the plot 
    # you could rename the file name if you need
    cos_sim_res3 %>% write_delim(heatmapData, delim = '\t', col_names = T)

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