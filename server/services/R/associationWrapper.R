library(tidyverse)
library(jsonlite)
library(aws.s3)

# capture console output for all functions called in wrapper
wrapper <- function(fn, args, dataArgs) {
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  output = list()

  tryCatch({
    output = get(fn)(args, dataArgs)
  }, error = function(e) {
    output <<- append(output, list(uncaughtError = paste0(deparse(e$call), ': ', e$message)))
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list(stdout = stdout, output = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

getAssocVarData <- function(args, dataArgs) {
  setwd(dataArgs$wd)
  fullDataPath = paste0(dataArgs$savePath, 'vardata_refdata_selected.txt')

  association_data_file <- paste0(dataArgs$s3Data, 'Association/', args$study, '_vardata.RData')

  tryCatch({
    s3load(association_data_file, dataArgs$bucket)
  }, error = function(e) {
    return(list(error = "ERROR: association variable data are not avaiable for selected study. please check input or select another study"))
  })

  ## extract the variable information
  vardata_refdata_selected <- vardata_refdata %>% filter(Cancer_Type == args$cancer)
  # clist will be used for the Assocaition Variable Data and Select Variables.
  clist <- vardata_refdata_selected %>% select(data_source, data_type, variable_name, variable_value_type) %>% unique() %>% arrange(variable_name)

  # save vardata_refdata_selected for download
  vardata_refdata_selected %>% write_delim(file = fullDataPath, delim = '\t', col_names = T, na = '')

  return(list(assocVarData = clist, fullDataPath = fullDataPath))
}

getExpVarData <- function(args, dataArgs) {
  # load exposure data files
  exposure_data_file <- paste0(dataArgs$s3Data, 'Exposure/', args$study, "_", args$strategy, '_exposure_refdata.RData')
  association_data_file <- paste0(dataArgs$s3Data, 'Association/', args$study, '_vardata.RData')

  tryCatch({
    s3load(exposure_data_file, dataArgs$bucket)
    s3load(association_data_file, dataArgs$bucket)
  }, error = function(e) {
    return(list(error = "ERROR: Exposure or association variable data are not avaiable for selected study. please check input or select another study"))
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

  exposure_refdata_selected <- exposure_refdata_selected %>% filter(Sample %in% osamples)

  # expsorue variable list
  Exposure_varlist <- colnames(exposure_refdata_selected)[-c(1:2)]

  return(list(expVarList = Exposure_varlist))
}


loadCollapse <- function(args, dataArgs) {
  # load exposure data files
  exposure_data_file <- paste0(dataArgs$s3Data, 'Exposure/', args$study, "_", args$strategy, '_exposure_refdata.RData')
  association_data_file <- paste0(dataArgs$s3Data, 'Association/', args$study, '_vardata.RData')
  s3load(exposure_data_file, dataArgs$bucket)
  s3load(association_data_file, dataArgs$bucket)

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

  exposure_refdata_selected <- exposure_refdata_selected %>% select(Sample, Signature_name, args$expVariable)

  vardata_refdata_selected <- vardata_refdata_selected %>%
      filter(data_source == args$variableSource, data_type == args$variableType, variable_name == args$assocVariable)

  if (unique(vardata_refdata_selected$variable_value_type) == "numeric") { vardata_refdata_selected$variable_value <- as.numeric(vardata_refdata_selected$variable_value) }

  vardata_refdata_selected <- vardata_refdata_selected %>%
      pivot_wider(id_cols = Sample, names_from = variable_name, values_from = variable_value)

  ### combined dataset
  data_input <- left_join(vardata_refdata_selected, exposure_refdata_selected) %>% select(-Sample)

  ## dropdown list for collapse_var1 and collapse_var2
  collapse_var1_list <- levels(data_input[[2]])
  collapse_var2_list <- levels(data_input[[3]])

  return(list(collapseVar1 = collapse_var1_list, collapseVar2 = collapse_var2_list))
}

univariable <- function(args, dataArgs) {
  source('services/R/Sigvisualfunc.R')
  setwd(dataArgs$wd)
  plotPath = paste0(dataArgs$savePath, 'association_result.svg')
  dataPath = paste0(dataArgs$savePath, 'asssociation_data.txt')
  assocTablePath = paste0(dataArgs$savePath, 'asssociation_test.txt')

  # load exposure data files
  exposure_data_file <- paste0(dataArgs$s3Data, 'Exposure/', args$study, "_", args$strategy, '_exposure_refdata.RData')
  association_data_file <- paste0(dataArgs$s3Data, 'Association/', args$study, '_vardata.RData')
  s3load(exposure_data_file, dataArgs$bucket)
  s3load(association_data_file, dataArgs$bucket)

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

  exposure_refdata_selected <- exposure_refdata_selected %>% select(Sample, Signature_name, args$variable2$name)

  vardata_refdata_selected <- vardata_refdata_selected %>%
      filter(data_source == args$variableSource, data_type == args$variableType, variable_name == args$variable1$name)

  if (unique(vardata_refdata_selected$variable_value_type) == "numeric") { vardata_refdata_selected$variable_value <- as.numeric(vardata_refdata_selected$variable_value) }

  vardata_refdata_selected <- vardata_refdata_selected %>%
      pivot_wider(id_cols = Sample, names_from = variable_name, values_from = variable_value)

  tryCatch({
    ## check data integration
    vardata_refdata_selected <- exposure_refdata_selected %>% select(Sample) %>% unique() %>% left_join(vardata_refdata_selected)
    ## including NA
    if (length(unique(vardata_refdata_selected[[2]])) == 1) {
      stop(paste0("mSigPortal Association failed: the selected variable name ", args$variable1$name, " have only unique value: ", unique(vardata_refdata_selected[[2]]), '.'))
    }
    tmpdata <- vardata_refdata_selected
    colnames(tmpdata)[2] <- 'Variable'
    tmpvalue <- tmpdata %>% count(Variable) %>% filter(n < 2) %>% dim() %>% .[[1]]

    if (tmpvalue != 0) {
      stop(paste0("mSigPortal Association failed: the selected variable name ", args$variable1$name, " have not enough obsevations for both levels."))
    }
  }, error = function(e) {
    return(list(error = e$message))
  })

  ### combined dataset
  data_input <- left_join(vardata_refdata_selected, exposure_refdata_selected) %>% select(-Sample)

  ## association test by group of signature name
  assocTable <- mSigPortal_associaiton_group(data = data_input, Group_Var = "Signature_name",
    Var1 = args$variable1$name, Var2 = args$variable2$name, type = args$testType,
    filter1 = args$variable1$filter, filter2 = args$variable2$filter,
    log1 = args$variable1$log2, log2 = args$variable2$log2, collapse_var1 = args$variable1$collapse,
    collapse_var2 = args$variable2$collapse)

  assocTable %>% write_delim(file = assocTablePath, delim = '\t', col_names = T, na = '')
  ## put result as a short table above the figure

  signature_name_list <- unique(assocTable[[1]]) ## dropdown list for the signature name
  signature_name_input <- if_else(args$signature != '', args$signature, signature_name_list[1]) ## by default, select the first signature name

  data_input <- data_input %>% filter(Signature_name == signature_name_input) %>% select(-Signature_name)


  mSigPortal_associaiton(data = data_input, Var1 = args$variable1$name, Var2 = args$variable2$name, type = args$testType,
    xlab = args$xlab, ylab = args$ylab, filter1 = args$variable1$filter, filter2 = args$variable2$filter,
    log1 = args$variable1$log2, log2 = args$variable2$log2, collapse_var1 = args$variable1$collapse,
    collapse_var2 = args$variable2$collapse, output_plot = plotPath)

  ## asssociation_data.txt will output as download text file.
  data_input %>% write_delim(file = dataPath, delim = '\t', col_names = T, na = '')

  return(list(plotPath = plotPath, dataPath = dataPath, assocTablePath = assocTablePath, dataTable = assocTable, signatureOptions = signature_name_list))
}


