known_error <- function(message) {
  cond <- structure(
    class = c("known_error", "error", "condition"),
    list(message = message)
  )
  stop(cond)
}

# Load RData from S3 and return with custom variable name
s3load_as <- function(object, bucket) {
    temp_env <- new.env()
    s3load(object, bucket, envir = temp_env)
    # Return the first (and typically only) object loaded
    return(get(ls(temp_env)[1], envir = temp_env))
}
