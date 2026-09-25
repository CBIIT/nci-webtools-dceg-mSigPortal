known_error <- function(message) {
  cond <- structure(
    class = c("known_error", "error", "condition"),
    list(message = message)
  )
  stop(cond)
}

# Guards downstream matrix code, which fails with an opaque "subscript out of bounds" on an empty set
assert_signature_set <- function(data, signature_set, profile = NULL, study = NULL) {
  if (nrow(data) > 0) {
    return(invisible(NULL))
  }
  profile_text <- if (!is.null(profile) && nzchar(profile)) paste0(profile, " ") else ""
  context <- if (!is.null(study) && nzchar(study)) paste0(" for study ", study) else ""
  known_error(paste0(
    "Signature set ", signature_set, " has no ", profile_text, "signatures", context,
    ". Please select another signature set."
  ))
}

# Load RData from S3 and return with custom variable name
s3load_as <- function(object, bucket) {
  temp_env <- new.env()
  s3load(object, bucket, envir = temp_env)
  # Return the first (and typically only) object loaded
  return(get(ls(temp_env)[1], envir = temp_env))
}

# Reads the refset data frame from an RData file; some study files store theirs as `signatures_refsets`
read_signature_refsets <- function(object, bucket) {
  env <- new.env()
  s3load(object, bucket, envir = env)
  available <- ls(env)
  preferred <- intersect(c("signature_refsets", "signatures_refsets"), available)
  get(if (length(preferred)) preferred[1] else available[1], envir = env)
}

# Reference signature sets combined with the selected study's de novo sets, if it has any
load_signature_refsets <- function(config, study = NULL, strategy = NULL) {
  reference <- read_signature_refsets(paste0(config$prefix, "Signature/signature_refsets.RData"), config$bucket)

  if (is.null(study) || is.null(strategy) || !nzchar(study) || !nzchar(strategy)) {
    return(reference)
  }

  study_file <- paste0(config$prefix, "Exposure/Study_Signatures/", study, "_", strategy, "_signature_refsets.RData")
  if (!aws.s3::object_exists(study_file, config$bucket)) {
    return(reference)
  }

  denovo <- read_signature_refsets(study_file, config$bucket)
  # a study set name that also exists in the reference file wins, so the study's own data is used
  reference <- reference %>% filter(!Signature_set_name %in% unique(denovo$Signature_set_name))

  bind_rows(reference, denovo[, intersect(colnames(denovo), colnames(reference))])
}
