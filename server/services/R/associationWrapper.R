library(tidyverse)
library(jsonlite)
library(aws.s3)
library(scales)

# capture console output for all functions called in wrapper
wrapper <- function(fn, args, config) {
  stdout <- vector('character')
  con <- textConnection('stdout', 'wr', local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  output = list()

  tryCatch({
    output = get(fn)(args, config)
  }, error = function(e) {
    print(e)
    output <<- append(output, list(uncaughtError = e$message))
  }, finally = {
    sink(con)
    sink(con)
    return(toJSON(list(stdout = stdout, output = output), pretty = TRUE, auto_unbox = TRUE))
  })
}

getAssocVarData <- function(args, config) {
  setwd(config$wd)
  fullDataPath = paste0(config$savePath, 'vardata_refdata_selected.txt')

  association_data_file <- paste0(config$s3Data, 'Association/', args$study, '_', args$strategy, '_', args$cancer, '_vardata.RData')

  tryCatch({
    s3load(association_data_file, config$bucket)
  }, error = function(e) {
    stop("ERROR: association variable data are not avaiable for selected study. please check input or select another study")
  })

  ## extract the variable information
  vardata_refdata_selected <- vardata_refdata %>% filter(Cancer_Type == args$cancer)
  # clist will be used for the Assocaition Variable Data and Select Variables.
  clist <- vardata_refdata_selected %>% select(data_source, data_type, variable_name, variable_value_type) %>% unique() %>% arrange(variable_name)

  # save vardata_refdata_selected for download
  vardata_refdata_selected %>% write_delim(file = fullDataPath, delim = '\t', col_names = T, na = '')
  return(list(assocVarData = clist, fullDataPath = fullDataPath))
}

getExpVarData <- function(args, config) {
  # load exposure data files
  exposure_data_file <- paste0(config$s3Data, 'Exposure/', args$study, "_", args$strategy, '_exposure_refdata.RData')
  association_data_file <- paste0(config$s3Data, 'Association/', args$study, '_', args$strategy, '_', args$cancer, '_vardata.RData')

  tryCatch({
    s3load(exposure_data_file, config$bucket)
    s3load(association_data_file, config$bucket)
  }, error = function(e) {
    stop("ERROR: Exposure or association variable data are not avaiable for selected study. please check input or select another study")
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


loadCollapse <- function(args, config) {
  source('services/R/Sigvisualfunc.R')
  # load exposure data files
  exposure_data_file <- paste0(config$s3Data, 'Exposure/', args$study, "_", args$strategy, '_exposure_refdata.RData')
  association_data_file <- paste0(config$s3Data, 'Association/', args$study, '_', args$strategy, '_', args$cancer, '_vardata.RData')
  s3load(exposure_data_file, config$bucket)
  s3load(association_data_file, config$bucket)

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

  exposure_refdata_selected <- exposure_refdata_selected %>% select(Sample, Signature_name, args$expName)

  vardata_refdata_selected <- vardata_refdata_selected %>%
    filter(data_source == args$source, data_type == args$type, variable_name == args$assocName)

  if (unique(vardata_refdata_selected$variable_value_type) == "numeric") {
    vardata_refdata_selected$variable_value <- as.numeric(vardata_refdata_selected$variable_value)
  }

  vardata_refdata_selected <- vardata_refdata_selected %>%
    pivot_wider(id_cols = Sample, names_from = variable_name, values_from = variable_value)

  ## check data integration
  vardata_refdata_selected <- exposure_refdata_selected %>% select(Sample) %>% unique() %>% left_join(vardata_refdata_selected)
  ## including NA
  if (length(unique(vardata_refdata_selected[[2]])) == 1) {
    error = paste0("mSigPortal Association failed: the selected variable name ", args$assocName, " only has unique value: ", unique(vardata_refdata_selected[[2]]), '.')
    return(list(error = error))
  }
  tmpdata <- vardata_refdata_selected
  colnames(tmpdata)[2] <- 'Variable'
  # tmpvalue <- tmpdata %>% count(Variable) %>% filter(n < 2) %>% dim() %>% .[[1]]

  # if (tmpvalue != 0) {
  #   error = paste0("mSigPortal Association failed: the selected variable name ", args$assocName, " does not have enough obsevations for both levels.")
  #   return(list(error = error))
  # }

  ## revise for category only...
    
    if(is.character(tmpdata$Variable)){
      tmpvalue <- tmpdata %>% count(Variable) %>% filter(n>2) %>% dim() %>% .[[1]]
      if(tmpvalue == 0){
        stop(paste0("mSigPortal Association failed: the selected variable name ",Association_varinput_name," have not enough obsevations for both levels."))
      }
    }

  ### combined dataset
  data_input <- left_join(vardata_refdata_selected, exposure_refdata_selected) %>% select(-Sample)
  data_input <- validate_vardf(data_input)

  ## dropdown list for collapse_var1 and collapse_var2
  collapse_var1_list <- levels(data_input[[args$assocName]])
  # collapse_var2_list <- unique(data_input[[args$expName]])

  return(list(collapseVar1 = collapse_var1_list)) #, collapseVar2 = collapse_var2_list))
}

univariable <- function(args, config) {
  source('services/R/Sigvisualfunc.R')
  setwd(config$wd)
  plotPath = paste0(config$savePath, 'association_result.svg')
  dataPath = paste0(config$savePath, 'asssociation_data.txt')
  assocTablePath = paste0(config$savePath, 'asssociation_test.txt')

  # load exposure data files
  exposure_data_file <- paste0(config$s3Data, 'Exposure/', args$study, "_", args$strategy, '_exposure_refdata.RData')
  association_data_file <- paste0(config$s3Data, 'Association/', args$study, '_', args$strategy, '_', args$cancer, '_vardata.RData')
  s3load(exposure_data_file, config$bucket)
  s3load(association_data_file, config$bucket)

 

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

  exposure_refdata_selected <- exposure_refdata_selected %>% select(Sample, Signature_name, args$exposureVar$name)

  vardata_refdata_selected <- vardata_refdata_selected %>%
    filter(data_source == args$associationVar$source, data_type == args$associationVar$type, variable_name == args$associationVar$name)

  if (unique(vardata_refdata_selected$variable_value_type) == "numeric") { vardata_refdata_selected$variable_value <- as.numeric(vardata_refdata_selected$variable_value) }

  vardata_refdata_selected <- vardata_refdata_selected %>%
    pivot_wider(id_cols = Sample, names_from = variable_name, values_from = variable_value)

  ## check data integration
  vardata_refdata_selected <- exposure_refdata_selected %>% select(Sample) %>% unique() %>% left_join(vardata_refdata_selected)
  ## including NA
  if (length(unique(vardata_refdata_selected[[2]])) == 1) {
    error = paste0("mSigPortal Association failed: the selected variable name ", args$assocName, " have only unique value: ", unique(vardata_refdata_selected[[2]]), '.')
    return(list(error = error))
  }
  tmpdata <- vardata_refdata_selected
  colnames(tmpdata)[2] <- 'Variable'
  # tmpvalue <- tmpdata %>% count(Variable) %>% filter(n < 2) %>% dim() %>% .[[1]]

  # if (tmpvalue != 0) {
  #   error = paste0("mSigPortal Association failed: the selected variable name ", args$assocName, " have not enough obsevations for both levels.")
  #   return(list(error = error))
  # }

  ## revise for category only...
    
    if(is.character(tmpdata$Variable)){
      tmpvalue <- tmpdata %>% count(Variable) %>% filter(n>2) %>% dim() %>% .[[1]]
      if(tmpvalue == 0){
        stop(paste0("mSigPortal Association failed: the selected variable name ",Association_varinput_name," have not enough obsevations for both levels."))
      }
    }

  ### combined dataset
  data_input <- left_join(vardata_refdata_selected, exposure_refdata_selected) %>% select(-Sample)

  ## association test by group of signature name
  # assocTable <- mSigPortal_associaiton_group(data = data_input, Group_Var = "Signature_name",
  #   Var1 = args$associationVar$name, Var2 = args$exposureVar$name, type = args$associationVar$type,
  #   filter1 = args$associationVar$filter, filter2 = args$exposureVar$filter,
  #   collapse_var1 = args$associationVar$collapse, collapse_var2 = NULL,
  #   log1 = args$associationVar$log2, log2 = args$exposureVar$log2
  #   )
assocTable <- mSigPortal_associaiton_group(data = data_input, Group_Var = "Signature_name",
    Var1 = args$associationVar$name, Var2 = args$exposureVar$name, type = args$testType,
    filter1 = args$associationVar$filter, filter2 = args$exposureVar$filter,
    collapse_var1 = args$associationVar$collapse, collapse_var2 = NULL,
    log1 = args$associationVar$log2, log2 = args$exposureVar$log2
    )


  assocTable %>% write_delim(file = assocTablePath, delim = '\t', col_names = T, na = '')
  ## put result as a short table above the figure

  signature_name_list <- unique(assocTable[[1]]) ## dropdown list for the signature name
  signature_name_input <- if_else(args$signature != '', args$signature, signature_name_list[1]) ## by default, select the first signature name

  data_input <- data_input %>% filter(Signature_name == signature_name_input) %>% select(-Signature_name)

  # mSigPortal_associaiton(data = data_input, Var1 = args$associationVar$name, Var2 = args$exposureVar$name, type = args$associationVar$type,
  #   xlab = args$xlab, ylab = args$ylab,
  #   filter1 = args$associationVar$filter, filter2 = args$exposureVar$filter,
  #   collapse_var1 = args$associationVar$collapse, collapse_var2 = NULL,
  #   log1 = args$associationVar$log2, log2 = args$exposureVar$log2,
  #   output_plot = plotPath)

   mSigPortal_associaiton(data = data_input, Var1 = args$associationVar$name, Var2 = args$exposureVar$name, type = args$testType,
    xlab = args$xlab, ylab = args$ylab,
    filter1 = args$associationVar$filter, filter2 = args$exposureVar$filter,
    collapse_var1 = args$associationVar$collapse, collapse_var2 = NULL,
    log1 = args$associationVar$log2, log2 = args$exposureVar$log2,
    output_plot = plotPath)

   ## asssociation_data.txt will output as download text file.
  data_input %>% write_delim(file = dataPath, delim = '\t', col_names = T, na = '')


  return(list(plotPath = plotPath, dataPath = dataPath, assocTablePath = assocTablePath, dataTable = assocTable, signatureOptions = signature_name_list))
}

loadCollapseMulti <- function(args, config) {
  require(broom)
  require(purrr)
  source('services/R/Sigvisualfunc.R')
  setwd(config$wd)
  # parse into array of objects
  associationVars = split(args$associationVars, 1:nrow(args$associationVars))

  # load exposure data files
  exposure_data_file <- paste0(config$s3Data, 'Exposure/', args$study, "_", args$strategy, '_exposure_refdata.RData')
  association_data_file <- paste0(config$s3Data, 'Association/', args$study, '_', args$strategy, '_', args$cancer, '_vardata.RData')
  s3load(exposure_data_file, config$bucket)
  s3load(association_data_file, config$bucket)

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

  exposure_refdata_selected <- exposure_refdata_selected %>% select(Sample, Signature_name, args$exposureVar$name)

  ### add more parameters according to user's input
  vardata_refdata_selected <- multivariable_inputs(vardata_refdata_selected, associationVars)
  data_input <- left_join(exposure_refdata_selected, vardata_refdata_selected) %>% select(-Sample)
  data_input <- validate_vardf(data_input)

  collapseOptions = map(associationVars, function(assocVar) levels(data_input[[assocVar$name]]))

  return(collapseOptions)
}

multivariable <- function(args, config) {
  require(broom)
  source('services/R/Sigvisualfunc.R')
  setwd(config$wd)
  plotPath = paste0(config$savePath, 'association_result.svg')
  dataPath = paste0(config$savePath, 'asssociation_data.txt')
  assocTablePath = paste0(config$savePath, 'asssociation_test.txt')
  # parse into array of objects
  associationVars = split(args$associationVars, 1:nrow(args$associationVars))

  # load exposure data files
  exposure_data_file <- paste0(config$s3Data, 'Exposure/', args$study, "_", args$strategy, '_exposure_refdata.RData')
  association_data_file <- paste0(config$s3Data, 'Association/', args$study, '_', args$strategy, '_', args$cancer, '_vardata.RData')
  s3load(exposure_data_file, config$bucket)
  s3load(association_data_file, config$bucket)

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

  exposure_refdata_selected <- exposure_refdata_selected %>% select(Sample, Signature_name, args$exposureVar$name)

  ### add more parameters according to user's input
  vardata_refdata_selected <- multivariable_inputs(vardata_refdata_selected, associationVars)
  data_input <- left_join(exposure_refdata_selected, vardata_refdata_selected) %>% select(-Sample)

  ## change variable name if detected special chacters ##
  colnames(data_input)[-c(1:2)] <- str_replace_all(str_replace(str_replace_all(colnames(data_input)[-c(1:2)], "[^[:alnum:]_ ]*", ""), "^[^[:alpha:]]*", ""), "  *", "_")

  rformula = paste0(args$exposureVar$name, " ~ ", paste0(colnames(data_input)[-c(1:2)], collapse = ' + '))
  ## regressionby group of signature name
  assocTable <- mSigPortal_associaiton_group(data = data_input, Group_Var = "Signature_name", type = args$testType, regression = TRUE, formula = rformula)
  assocTable %>% write_delim(file = assocTablePath, delim = '\t', col_names = T, na = '')
  # put result as a short table above the figure

  signature_name_list <- unique(assocTable[[1]]) ## drop-down list for the signature name
  signature_name_input <- if_else(args$signature != '', args$signature, signature_name_list[1]) ## by default, select the first signature name
  data_input <- data_input %>% filter(Signature_name == signature_name_input) %>% select(-Signature_name)

  mSigPortal_associaiton(data = data_input, type = args$testType, regression = TRUE, formula = rformula, output_plot = plotPath)

  data_input %>% write_delim(file = dataPath, delim = '\t', col_names = T, na = '')

  return(list(plotPath = plotPath, dataPath = dataPath, assocTablePath = assocTablePath, dataTable = assocTable, signatureOptions = signature_name_list))
}


