library(tidyverse)
# library(cowplot)
# library(ggsci)
# library(ggrepel)
# library(ggforce)
# library(ggtext)
# library(ggpubr)
library(jsonlite)
# library(stringr)
library(aws.s3)


loadData <- function(args, s3Data, localData, bucket) {
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()

    # load exposure data files
    exposure_data_file <- paste0(s3Data, 'Exposure/', args$study, "_", args$strategy, '_exposure_refdata.RData')
    association_data_file <- paste0(s3Data, 'Association/', args$study, '_vardata.RData')

    tryCatch({
      s3load(exposure_data_file, bucket)
      s3load(association_data_file, bucket)
    }, error = function(e) {
      msg = "ERROR: Exposure or assoicaiton variable data are not avaiable for selected study. please check input or select other study"
      print(msg)
      print(e)
      stop(msg)
    })

    exposure_refdata_selected <- exposure_refdata %>% filter(Signature_set_name == args$rsSet, Cancer_Type == args$cancer)

    exposure_refdata_selected <- exposure_refdata_selected %>%
        select(Sample, Signature_name, Signature_exposure = Exposure) %>%
        group_by(Sample) %>%
        mutate(Signature_exposure_ratio = Signature_exposure / sum(Signature_exposure), Signature_exposure_cat = if_else(Signature_exposure > 100, "Observed", "Not_observed")) %>%
        ungroup() %>%
        mutate(Signature_exposure_cat = factor(Signature_exposure_cat, level = c("Not_observed", "Observed")))

    vardata_refdata_selected <- vardata_refdata %>% filter(Cancer_Type == args$cancer)

    # overlapped samples
    osamples <- intersect(unique(vardata_refdata_selected$Sample), unique(exposure_refdata_selected$Sample))

    vardata_refdata_selected <- vardata_refdata_selected %>% filter(Sample %in% osamples)
    exposure_refdata_selected <- exposure_refdata_selected %>% filter(Sample %in% osamples)

    # expsorue variant list
    Exposure_varlist <- colnames(exposure_refdata_selected)[-c(1:2)]

    output = list('expVarList' = Exposure_varlist)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}


loadCollapse <- function(args, s3Data, localData, bucket) {
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    output = list()

    # load exposure data files
    exposure_data_file <- paste0(s3Data, 'Exposure/', args$study, "_", args$strategy, '_exposure_refdata.RData')
    association_data_file <- paste0(s3Data, 'Association/', args$study, '_vardata.RData')

    tryCatch({
      s3load(exposure_data_file, bucket)
      s3load(association_data_file, bucket)
    }, error = function(e) {
      msg = "ERROR: Exposure or assoicaiton variable data are not avaiable for selected study. please check input or select other study"
      print(msg)
      print(e)
      stop(msg)
    })

    exposure_refdata_selected <- exposure_refdata %>% filter(Signature_set_name == args$rsSet, Cancer_Type == args$cancer)

    exposure_refdata_selected <- exposure_refdata_selected %>%
        select(Sample, Signature_name, Signature_exposure = Exposure) %>%
        group_by(Sample) %>%
        mutate(Signature_exposure_ratio = Signature_exposure / sum(Signature_exposure), Signature_exposure_cat = if_else(Signature_exposure > 100, "Observed", "Not_observed")) %>%
        ungroup() %>%
        mutate(Signature_exposure_cat = factor(Signature_exposure_cat, level = c("Not_observed", "Observed")))

    vardata_refdata_selected <- vardata_refdata %>% filter(Cancer_Type == args$cancer)

    # overlapped samples
    osamples <- intersect(unique(vardata_refdata_selected$Sample), unique(exposure_refdata_selected$Sample))

    vardata_refdata_selected <- vardata_refdata_selected %>% filter(Sample %in% osamples)
    exposure_refdata_selected <- exposure_refdata_selected %>% filter(Sample %in% osamples)

    exposure_refdata_selected <- exposure_refdata_selected %>% select(Sample, args$expVar)

    vardata_refdata_selected <- vardata_refdata_selected %>%
      filter(data_source == args$dataSource, data_type == args$dataType, variable_name == args$assocVar)

    if (unique(vardata_refdata_selected$variable_value_type) == "numeric") { vardata_refdata_selected$variable_value <- as.numeric(vardata_refdata_selected$variable_value) }

    vardata_refdata_selected <- vardata_refdata_selected %>%
    pivot_wider(id_cols = Sample, names_from = variable_name, values_from = variable_value)

    data_input <- left_join(vardata_refdata_selected, exposure_refdata_selected)

    ## dropdown list for collapse_var1 and collapse_var2 
    collapse_var1_list <- levels(data_input[[2]])
    collapse_var2_list <- levels(data_input[[3]])


    output = list('collapseVar1' = collapse_var1_list, 'collapseVar2' = collapse_var2_list)
  }, error = function(e) {
    print(e)
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

calculate <- function(args, projectID, rootDir, savePath, s3Data, localData, bucket) {
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  tryCatch({
    source('services/R/Sigvisualfunc.R')
    output = list()
    plotPath = paste0(savePath, "/association_result.svg")
    dataPath = paste0(savePath, '/asssociation_data.txt')

    # load exposure data files
    exposure_data_file <- paste0(s3Data, 'Exposure/', args$study, "_", args$strategy, '_exposure_refdata.RData')
    association_data_file <- paste0(s3Data, 'Association/', args$study, '_vardata.RData')

    tryCatch({
      s3load(exposure_data_file, bucket)
      s3load(association_data_file, bucket)
    }, error = function(e) {
      msg = "ERROR: Exposure or assoicaiton variable data are not avaiable for selected study. please check input or select other study"
      print(msg)
      print(e)
      stop(msg)
    })

    exposure_refdata_selected <- exposure_refdata %>% filter(Signature_set_name == args$rsSet, Cancer_Type == args$cancer)

    exposure_refdata_selected <- exposure_refdata_selected %>%
        select(Sample, Signature_name, Signature_exposure = Exposure) %>%
        group_by(Sample) %>%
        mutate(Signature_exposure_ratio = Signature_exposure / sum(Signature_exposure), Signature_exposure_cat = if_else(Signature_exposure > 100, "Observed", "Not_observed")) %>%
        ungroup() %>%
        mutate(Signature_exposure_cat = factor(Signature_exposure_cat, level = c("Not_observed", "Observed")))

    vardata_refdata_selected <- vardata_refdata %>% filter(Cancer_Type == args$cancer)

    # overlapped samples
    osamples <- intersect(unique(vardata_refdata_selected$Sample), unique(exposure_refdata_selected$Sample))

    vardata_refdata_selected <- vardata_refdata_selected %>% filter(Sample %in% osamples)
    exposure_refdata_selected <- exposure_refdata_selected %>% filter(Sample %in% osamples)

    # expsorue variant list
    Exposure_varlist <- colnames(exposure_refdata_selected)[-c(1:2)]

    exposure_refdata_selected <- exposure_refdata_selected %>% select(Sample, args$expVar)
    vardata_refdata_selected <- vardata_refdata_selected %>%
      filter(data_source == args$dataSource, data_type == args$dataType, variable_name == args$assocVar)

    if (unique(vardata_refdata_selected$variable_value_type) == "numeric") { vardata_refdata_selected$variable_value <- as.numeric(vardata_refdata_selected$variable_value) }

    vardata_refdata_selected <- vardata_refdata_selected %>%
      pivot_wider(id_cols = Sample, names_from = variable_name, values_from = variable_value)

    data_input <- left_join(vardata_refdata_selected, exposure_refdata_selected)

    mSigPortal_associaiton(data = data_input, Var1 = args$assocVar, Var2 = args$expVar, type = args$testType,
      xlab = args$xlab, ylab = args$ylab, filter_zero1 = args$filter1, filter_zero2 = args$filter2,
      log1 = args$log2_1, log2 = args$log2_2, collapse_var1 = args$collapse1,
      collapse_var2 = args$collapse2, output_plot = plotPath)

    ## asssociation_data.txt will output as download text file. 
    data_input %>% write_delim(file = dataPath, delim = '\t', col_names = T, na = '')

    output = list('plotPath' = plotPath, 'dataPath' = dataPath)
  }, error = function(e) {
    print(e)
    output[['error']] <<- e
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list('stdout' = stdout, 'output' = output), pretty = TRUE, auto_unbox = TRUE))
  })
}


