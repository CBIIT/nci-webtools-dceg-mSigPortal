# converts relative path into static server fetchable path
# e.g. ../tmp/6c817baf-8d16-4fce-9750-f55cf60eade2/results/univariable/association_result.svg to /results/univariable/association_result.svg
getResultsPath <- function(path) {
  return(stringr::str_match(path, '\\/results\\/.*')[1])
}