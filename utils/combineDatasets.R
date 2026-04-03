library(dplyr, warn.conflicts = FALSE)
library(vroom)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
    stop("Usage: Rscript combineDatasets.R 'inputFileGlobPattern' processorFunction outputFile.csv")
}

inputFiles <- unlist(lapply(strsplit(args[1], " ")[[1]], Sys.glob))
processorFunction <- args[2]
outputFile <- args[3]

extract <- function(str, separator, index) {
    unlist(strsplit(str, separator))[index]
}

regex_extract <- function(str, pattern) {
    regmatches(str, regexpr(pattern, str))
}

regex_group <- function(str, pattern) {
    regmatches(str, regexec(pattern, str))[[1]][2]
}

combineAssociationFiles <- function(x) {
    if (!"icgc_specimen_id" %in% names(x)) x$icgc_specimen_id <- NA_character_
    if (!"icgc_donor_id" %in% names(x)) x$icgc_donor_id <- NA_character_
    x %>%
        mutate(
            Study = sub("_.*", "", basename(filepath)),
            Dataset = sub("^[^_]+_([^_]+)_.*", "\\1", basename(filepath)),
            .before = Cancer_Type
        ) %>%
        rename(
            study = Study,
            strategy = Dataset,
            cancer = Cancer_Type,
            sample = Sample,
            icgcSpecimenId = icgc_specimen_id,
            icgcDonorId = icgc_donor_id,
            dataSource = data_source,
            dataType = data_type,
            variableName = variable_name,
            variableValue = variable_value,
            variableValueType = variable_value_type
        ) %>%
        select(study, strategy, cancer, sample, icgcSpecimenId, icgcDonorId, dataSource, dataType, variableName, variableValue, variableValueType)
}

combineExposureFiles <- function(x) {
    x %>%
        mutate(study = sub("_.*", "", basename(filepath))) %>%
        rename(
            strategy = Dataset,
            cancer = Cancer_Type,
            organ = Organ,
            sample = Sample,
            signatureSetName = Signature_set_name,
            signatureName = Signature_name,
            exposure = Exposure
        ) %>%
        select(study, strategy, cancer, organ, sample, signatureSetName, signatureName, exposure)
}

combineSeqmatrixFiles <- function(x) {
    x %>%
        mutate(
            profile = regex_extract(Profile, "^[A-Z]+"),
            matrix = regex_extract(Profile, "[0-9]+$"),
            study = sub("_.*", "", basename(filepath))
        ) %>%
        rename(
            cancer = Cancer_Type,
            sample = Sample,
            strategy = Dataset,
            mutationType = MutationType,
            mutations = Mutations
        ) %>%
        select(study, cancer, sample, strategy, profile, matrix, mutationType, mutations)
}

combineSignatureFiles <- function(x) {
    x %>%
        mutate(
            profile = regex_extract(Profile, "^[A-Z]+"),
            matrix = regex_extract(Profile, "[0-9]+$"),
            .before = Profile
        ) %>%
        rename(
            source = Source,
            signatureSetName = Signature_set_name,
            strategy = Dataset,
            strandInfo = Strand_info,
            strand = Strand,
            signatureName = Signature_name,
            mutationType = MutationType,
            contribution = Contribution
        ) %>%
        select(
            -Profile
        )
}

signatureSummary <- function(x) {
    x %>%
        mutate(
            profile = regex_extract(Profile, "^[A-Z]+"),
            matrix = regex_extract(Profile, "[0-9]+$"),
            .before = Profile
        ) %>%
        rename(
            species = Species,
            signatureSetName = Signature_set_name,
            count = N
        ) %>%
        select(
            -Profile
        )
}

combinePatternFiles <- function(x) {
    x %>%
        rowwise() %>%
        mutate(
            study = regex_group(Study, "^(.+)@"),
            cancer = regex_group(Study, "@(.+)$"),
            .before = Study
        ) %>%
        rename(
            sample = Sample,
            total = Total,
            pattern = Pattern,
            n0 = N0,
            n1 = N1,
            n2 = N2
        ) %>%
        select(
            -Study
        )
}

combineEtiology <- function(x) {
    x %>%
        rename(
            study = Study,
            strategy = Dataset,
            cancer = Cancer_Type,
            organ = Organ,
            sample = Sample,
            signatureSetName = Signature_set_name,
            mutations = Total_Mutations,
            cosineSimilarity = Cosine_similarity,
            sampleSize = Sample_size,
            signatureName = Signature_name,
            exposure = Exposure,
            burden = Burden,
            signatureSize = Signature_Size
        ) %>%
        select(
            -Study_Name,
            -Sample_Names
        )
}

combineEtiologyOrgan <- function(x) {
    x %>% rename(
        signature = Signature,
        cohort = Cohort,
        organ = Organ,
        prevalence = Prevalence,
        organSpecificSignature = "Organ-Specific Signature",
        contribution = Contribution
    )
}

combineEtiologySignature <- function(x) {
    x %>%
        rowwise() %>%
        mutate(
            etiology = stringr::str_split(Etiology, "\n")[[1]][1],
            author = stringr::str_split(Etiology, "\n")[[1]][2],
            profile = regex_extract(Profile, "^[A-Z]+"),
            matrix = regex_extract(Profile, "[0-9]+$"),
            .before = Profile
        ) %>%
        rename(
            signatureSetName = Signature_set_name,
            signatureName = Signature_name,
            mutationType = MutationType,
            contribution = Contribution
        ) %>%
        select(
            -Etiology,
            -Profile
        )
}

combineRefgenome <- function(x) {
    x %>%
        rowwise() %>%
        mutate(
            chr = regex_group(chr, "chr(.*)"),
        ) %>%
        rename(
            start = start2,
            end = end2,
        )
}

first <- TRUE
for (f in inputFiles) {
    env <- new.env()
    load(f, envir = env)
    dataset <- get(ls(env)[1], envir = env)
    dataset$filepath <- f
    processed <- do.call(processorFunction, list(x = dataset)) %>% select(-any_of("filepath"))
    vroom_write(processed, file = outputFile, delim = ",", na = "", append = !first)
    first <- FALSE
    rm(env, dataset, processed)
    gc()
}
