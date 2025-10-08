library(tidyverse)
library(GenomicRanges)
library(Biostrings)
library(SATS)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(jsonlite)

# Load the refitting function
source("/app/refitting-service/refitting.R")

# capture console output for all functions called in wrapper
wrapper <- function(fn, args, config = list()) {
  stdout <- vector("character")
  con <- textConnection("stdout", "wr", local = TRUE)
  sink(con, type = "message")
  sink(con, type = "output")

  output <- list()

  tryCatch(
    {
      output <- get(paste0("msigportal.", fn))(args, config)
    },
    error = function(e) {
      print(e)
      output <<- append(output, list(uncaughtError = e$message))
    },
    finally = {
      sink(type = "message")
      sink(type = "output")
      close(con)
    }
  )

  return(append(output, list(stdout = stdout)))
}

# Main refitting function that processes form inputs and calls the R function
msigportal.refitSBS <- function(args, config = list()) {
  
  # Extract form inputs
  maf_file <- args$mafFile
  genomic_file <- args$genomicFile
  clinical_file <- args$clinicalFile
  output_dir <- args$outputDir %||% "/tmp/refitting_output"
  
  # Form parameters
  signature_type <- args$signatureType %||% "SBS"
  reference_genome <- args$referenceGenome %||% "hg19"
  job_name <- args$jobName %||% paste0("refitting_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  email <- args$email
  
  # Validate signature type
  if (signature_type != "SBS") {
    stop("Currently only SBS signature type is supported")
  }
  
  # Validate required files
  if (is.null(maf_file) || is.null(genomic_file) || is.null(clinical_file)) {
    stop("Missing required files: MAF file, Genomic file, and Clinical file are all required")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Generate output filename
  out_file <- paste0(job_name, "_H_Burden_est.csv")
  
  # Call the main refitting function
  tryCatch({
    results <- run_sbs_refitting(
      maf_file = maf_file,
      genomic_file = genomic_file,
      clinical_file = clinical_file,
      output_dir = output_dir,
      genome = reference_genome,
      save_csv = TRUE,
      out_file = out_file,
      match_on_oncotree = FALSE
    )
    
    # Format output for the targeted sequencing table
    formatted_output <- list(
      success = TRUE,
      jobName = job_name,
      email = email,
      signatureType = signature_type,
      referenceGenome = reference_genome,
      outputFile = file.path(output_dir, out_file),
      
      # Results for targeted sequencing table
      h_burden_results = list(
        data = results$H_Burden,
        total_samples = ncol(results$V),
        total_signatures = length(unique(results$H_Burden$Signature)),
        samples_processed = ncol(results$V)
      ),
      
      # Additional results data
      panel_context = if (!is.null(results$Panel_context)) as.data.frame(results$Panel_context) else NULL,
      v_matrix = if (!is.null(results$V)) as.data.frame(results$V) else NULL,
      l_matrix = if (!is.null(results$L)) as.data.frame(results$L) else NULL,
      
      # Summary information
      summary = list(
        analysis_type = "SBS_Refitting",
        genome_version = reference_genome,
        total_samples = ncol(results$V),
        unique_signatures = length(unique(results$H_Burden$Signature)),
        total_signature_activities = nrow(results$H_Burden),
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        output_csv_path = file.path(output_dir, out_file)
      ),
      
      # Metadata for the frontend
      metadata = list(
        processing_completed = TRUE,
        files_processed = list(
          maf_file = basename(maf_file),
          genomic_file = basename(genomic_file),
          clinical_file = basename(clinical_file)
        ),
        common_files_used = list(
          signature_reference = "Alex_Sigs_TMB_check_V3.4_SBS2_13 together.csv",
          cancer_dictionary = "0 Cancer_Dictionary_BZ.csv"
        )
      )
    )
    
    return(formatted_output)
    
  }, error = function(e) {
    # Return error information
    return(list(
      success = FALSE,
      error = paste("SBS Refitting Error:", e$message),
      jobName = job_name,
      email = email,
      signatureType = signature_type,
      referenceGenome = reference_genome,
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ))
  })
}

# Function to validate inputs before processing
msigportal.validateRefittingInputs <- function(args, config = list()) {
  
  validation_results <- list(
    valid = TRUE,
    errors = list(),
    warnings = list()
  )
  
  # Check required files
  required_files <- c("mafFile", "genomicFile", "clinicalFile")
  
  for (file_param in required_files) {
    if (is.null(args[[file_param]])) {
      validation_results$errors <- append(validation_results$errors, 
                                         paste("Missing required file:", file_param))
      validation_results$valid <- FALSE
    } else if (!file.exists(args[[file_param]])) {
      validation_results$errors <- append(validation_results$errors, 
                                         paste("File not found:", args[[file_param]]))
      validation_results$valid <- FALSE
    }
  }
  
  # Validate signature type
  if (!is.null(args$signatureType) && args$signatureType != "SBS") {
    validation_results$errors <- append(validation_results$errors, 
                                       "Currently only SBS signature type is supported")
    validation_results$valid <- FALSE
  }
  
  # Validate genome parameter
  if (!is.null(args$referenceGenome) && !args$referenceGenome %in% c("hg19", "hg38")) {
    validation_results$errors <- append(validation_results$errors, 
                                       "Invalid reference genome. Must be 'hg19' or 'hg38'")
    validation_results$valid <- FALSE
  }
  
  # Validate email format
  if (!is.null(args$email) && !grepl("^[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\\.[A-Za-z]{2,}$", args$email)) {
    validation_results$errors <- append(validation_results$errors, 
                                       "Invalid email format")
    validation_results$valid <- FALSE
  }
  
  return(validation_results)
}

# Function to get job status (for monitoring)
msigportal.getRefittingStatus <- function(args, config = list()) {
  
  job_id <- args$jobId
  output_dir <- args$outputDir %||% "/tmp/refitting_output"
  
  if (is.null(job_id)) {
    return(list(
      status = "error",
      message = "Missing jobId parameter"
    ))
  }
  
  # Look for output files
  expected_output <- file.path(output_dir, paste0(job_id, "_H_Burden_est.csv"))
  
  if (file.exists(expected_output)) {
    # Read the results file to get summary info
    tryCatch({
      results_data <- read.csv(expected_output, stringsAsFactors = FALSE)
      
      return(list(
        status = "completed",
        message = "SBS refitting analysis completed successfully",
        outputFile = expected_output,
        fileSize = file.size(expected_output),
        lastModified = as.character(file.mtime(expected_output)),
        totalRows = nrow(results_data),
        uniqueSamples = length(unique(results_data$SAMPLE_ID)),
        uniqueSignatures = length(unique(results_data$Signature))
      ))
    }, error = function(e) {
      return(list(
        status = "error",
        message = paste("Error reading results file:", e$message)
      ))
    })
  } else {
    return(list(
      status = "running",
      message = "SBS refitting analysis in progress",
      expectedOutput = expected_output
    ))
  }
}

# Utility function for null coalescing
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}
