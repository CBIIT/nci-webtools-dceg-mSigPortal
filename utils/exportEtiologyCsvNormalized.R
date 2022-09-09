library(dplyr, warn.conflicts = FALSE)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript exportEtiologyCsvNormalized.R inputFile.json outputFile.csv")
}

inputFile <- args[1]
processorFunction <- args[2]
outputFile <- args[3]

loadJson <- function(filepath) {
  dataset <- read_json(filepath, simplifyVector = T) %>%
    rowwise() %>%
    mutate(across(where(~ class(.x) == "list"), ~ paste(., collapse = " ")))
  dataset
}

CancerSpecificSignatures_2022 <- function(x) {
  x %>%
    mutate(
      category = processorFunction,
      etiologyDisplay = "Proposed Etiologies",
    ) %>%
    rename(
      etiology = Etiology,
      signature = "Signature Name",
      study = Study,
      studyUrl = Study_URL,
      genomeBuild = "Genome Build",
      signatureSource = "Signature Source",
      sourceUrl = Source_URL,
      cohort = "Observed Cohort",
      cosmic_v3_2 = "Identified in COSMICv3.2",
      refSig_v1 = "Identified in RefSigv1",
      descriptionStrandBias = "Replication Strand Biases",
      descriptionBaseContext = "Base Context",
      description = Description,
      note = Note,
      cisMutation = "Mutation Identified in Cis"
    ) %>%
    mutate(json = toJSON(unbox(across()))) %>%
    select(category, etiology, signature, json)
}

CancerSpecificSignatures <- function(x) {
  x %>%
    mutate(
      category = processorFunction,
      etiologyDisplay = "Tissue Types",
    ) %>%
    rename(
      etiology = Etiology,
      signature = "Tissue Specific Signature",
      referenceSignature = "Ref Signature",
      refSigProportion = "RefSig Proportion",
      study = Study,
      studyUrl = Study_URL,
      source = Source,
      sourceUrl = Source_URL,
      description = Description,
    ) %>%
    mutate(json = toJSON(unbox(across()))) %>%
    select(category, etiology, signature, json)
}

Cosmic <- function(x) {
  x %>%
    mutate(
      category = processorFunction,
      etiologyDisplay = "Proposed Etiologies"
    ) %>%
    rename(
      etiology = Etiology,
      signature = "Signature Name",
      study = Study,
      studyUrl = Study_URL,
      genomeBuild = "Genome Build",
      signatureSource = "Signature Source",
      sourceUrl = Source_URL,
      description = Description,
      tissueDistribution = Tissue_Distribution
    ) %>%
    mutate(json = toJSON(unbox(across()))) %>%
    select(category, etiology, signature, json)
}

EnviromentalMutagenesis <- function(x) {
  x %>%
    mutate(
      category = processorFunction,
      etiologyDisplay = "Proposed Mutagens"
    ) %>%
    rename(
      etiology = Etiology,
      signature = Signature,
      mutagen = Mutagen,
      treatment = Treatment,
      study = Study,
      studyUrl = Study_URL,
      source = Source,
      sourceUrl = Source_URL,
    ) %>%
    mutate(json = toJSON(unbox(across()))) %>%
    select(category, etiology, signature, json)
}

GeneEdits <- function(x) {
  x %>%
    mutate(
      category = processorFunction,
      etiologyDisplay = "Genes"
    ) %>%
    rename(
      etiology = Etiology,
      signature = Signature,
      cellLine = "Cell Line",
      study = Study,
      studyUrl = Study_URL,
      source = Source,
      sourceUrl = Source_URL,
    ) %>%
    mutate(json = toJSON(unbox(across()))) %>%
    select(category, etiology, signature, json)
}

CancerTherapies <- function(x) {
  x %>%
    mutate(
      category = processorFunction,
      etiologyDisplay = "Treatments"
    ) %>%
    rename(
      etiology = Treatments,
      signature = Signature,
      signatureExtractionMethod = Siganture_Extraction_Method,
      tumorType = Tumor_Type,
      study = Study,
      studyUrl = Study_URL,
    ) %>%
    mutate(json = toJSON(unbox(across()))) %>%
    select(category, etiology, signature, json)
}

Others <- function(x) {
  x %>%
    mutate(
      category = processorFunction,
      etiologyDisplay = "Proposed Etiologies"
    ) %>%
    rename(
      etiology = Etiology,
      signature = "Signature Name",
      signatureSource = "Signature Source",
      sourceUrl = URL,
      description = Description
    ) %>%
    mutate(json = toJSON(unbox(across()))) %>%
    select(category, etiology, signature, json)
}

dataset <- loadJson(inputFile)
processedDataset <- do.call(processorFunction, list(x = dataset))
write.csv(processedDataset, file = outputFile, na = "", row.names = F)