library(dplyr, warn.conflicts = FALSE)
library(vroom)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
    stop("Usage: Rscript combineDatasets.R 'inputFileGlobPattern' processorFunction outputFile.csv")
}

inputFiles <- Sys.glob(args[1])
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
    x %>%
        rowwise() %>%
        mutate(
            Study = extract(basename(filepath), "_", 1),
            Dataset = extract(basename(filepath), "_", 2),
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
        )
}

combineExposureFiles <- function(x) {
    x %>% rename(
        study = Study,
        strategy = Dataset,
        cancer = Cancer_Type,
        organ = Organ,
        sample = Sample,
        signatureSetName = Signature_set_name,
        signatureName = Signature_name,
        exposure = Exposure
    )
}

combineSeqmatrixFiles <- function(x) {
    x %>%
        mutate(
            profile = regex_extract(Profile, "^[A-Z]+"),
            matrix = regex_extract(Profile, "[0-9]+$"),
            .before = Profile
        ) %>%
        rename(
            study = Study,
            cancer = Cancer_Type,
            sample = Sample,
            strategy = Dataset,
            mutationType = MutationType,
            mutations = Mutations
        ) %>%
        select(
            -Profile
        )
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

datasets <- sapply(inputFiles, function(f) get(load(f)), simplify = F)
combinedDatasets <- bind_rows(datasets, .id = "filepath")
processedDatasets <- do.call(processorFunction, list(x = combinedDatasets)) %>% select(-filepath)
vroom_write(processedDatasets, file = outputFile, delim = ",", na = "")
