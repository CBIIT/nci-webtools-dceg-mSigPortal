source("renv/activate.R")

# Apply CRAN repository options only for macOS or Windows
if (.Platform$OS.type == "windows" || grepl("darwin", R.version$os)) {
    options(repos = c(CRAN = "https://packagemanager.posit.co/cran/latest"))
}
