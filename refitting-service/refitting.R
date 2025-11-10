# refitting.R
# Run SBS refitting end-to-end given three input files + shared reference files.

# ---- Simple test function for r-wrapper ----
sum <- function(a, b) {
  return(a + b)
}

# ---- libraries used by the function ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(Biostrings)
  library(SATS)
  library(BSgenome)  # for getSeq()
})

# Safety: if someone previously masked base::sum, unmask it
if (exists("sum") && !identical(sum, base::sum)) {
  rm(sum)
}

# ---- main function ----
run_sbs_refitting <- function(maf_file,
                              genomic_file,
                              clinical_file,
                              output_dir,
                              common_files_dir,
                              genome = c("hg19","hg38"),
                              save_csv = TRUE,
                              out_file = "H_Burden_est.csv",
                              match_on_oncotree = FALSE) {

  # Add logging to track function execution
  cat("=== SBS Refitting Function Started ===\n")
  cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Input parameters:\n")
  cat("  maf_file:", maf_file, "\n")
  cat("  genomic_file:", genomic_file, "\n") 
  cat("  clinical_file:", clinical_file, "\n")
  cat("  output_dir:", output_dir, "\n")
  cat("  common_files_dir:", common_files_dir, "\n")
  cat("  genome:", genome, "\n")
  cat("  out_file:", out_file, "\n")
  cat("  match_on_oncotree:", match_on_oncotree, "\n")

  genome <- match.arg(genome)
  cat("Validated genome:", genome, "\n")

  # Check inputs exist
  cat("Checking input files exist...\n")
  user_files <- c(maf_file, genomic_file, clinical_file)
  missing_user <- user_files[!file.exists(user_files)]
  if (length(missing_user)) {
    cat("ERROR: Missing user files:\n", paste(missing_user, collapse = "\n"), "\n")
    stop("Missing user files:\n", paste(missing_user, collapse = "\n"))
  }
  cat("All user files found successfully\n")

  if (!dir.exists(output_dir)) {
    cat("ERROR: Output directory does not exist:", output_dir, "\n")
    stop("Output directory does not exist: ", output_dir)
  }
  cat("Output directory exists:", output_dir, "\n")
  
  if (missing(common_files_dir)) {
    cat("ERROR: common_files_dir must be provided\n")
    stop("common_files_dir must be provided (folder with the two shared CSVs).")
  }
  cat("Common files directory:", common_files_dir, "\n")

  # Shared reference files
  cat("Checking shared reference files...\n")
  ref_signature_file     <- file.path(common_files_dir, "Alex_Sigs_TMB_check_V3.4_SBS2_13 together.csv")
  cancer_dictionary_file <- file.path(common_files_dir, "0 Cancer_Dictionary_BZ.csv")
  common_needed <- c(ref_signature_file, cancer_dictionary_file)
  missing_common <- common_needed[!file.exists(common_needed)]
  if (length(missing_common)) {
    cat("ERROR: Missing common reference files:\n", paste(missing_common, collapse = "\n"), "\n")
    stop("Missing common reference files:\n", paste(missing_common, collapse = "\n"))
  }
  cat("All reference files found successfully\n")

  # Pick genome (package must be installed, but we don't attach it globally)
  Hsapiens <- switch(
    genome,
    hg19 = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
    hg38 = BSgenome.Hsapiens.UCSC.hg38::Hsapiens
  )

  # Load RefTMB dataset from SATS (dataset, not an exported object)
  utils::data("RefTMB", package = "SATS", envir = environment())
  if (!exists("RefTMB", envir = environment())) {
    stop("Failed to load RefTMB from SATS. Check SATS installation.")
  }

  # ---- read inputs ----
  SBS_TMBv3.4 <- utils::read.csv(ref_signature_file, row.names = 1, check.names = FALSE)
  Mut_category_order <- rownames(SBS_TMBv3.4)
  if (length(Mut_category_order) != 96) stop("Expected 96 SBS categories in reference signature file.")

  mutations2 <- utils::read.table(maf_file, header = TRUE, sep = "\t", quote = "",
                                  stringsAsFactors = FALSE, check.names = FALSE)
  genomic_information <- utils::read.table(genomic_file, header = TRUE, sep = "\t", quote = "",
                                           stringsAsFactors = FALSE, check.names = FALSE)
  clinical_sample <- utils::read.table(clinical_file, header = TRUE, sep = "\t", quote = "",
                                       stringsAsFactors = FALSE, check.names = FALSE)
  annoFile <- utils::read.csv(cancer_dictionary_file, check.names = FALSE)

  # Column checks
  need_mut <- c("Chromosome","Start_Position","Variant_Type","Reference_Allele","Tumor_Seq_Allele2","Tumor_Sample_Barcode")
  miss_mut <- setdiff(need_mut, names(mutations2))
  if (length(miss_mut)) stop("MAF missing columns: ", paste(miss_mut, collapse = ", "))

  need_cli <- c("SAMPLE_ID","SEQ_ASSAY_ID","CANCER_TYPE")
  miss_cli <- setdiff(need_cli, names(clinical_sample))
  if (length(miss_cli)) stop("Clinical file missing columns: ", paste(miss_cli, collapse = ", "))

  need_anno <- c("CANCER_TYPE","cancerAnalyzed")
  miss_anno <- setdiff(need_anno, names(annoFile))
  if (length(miss_anno)) stop("Annotation file missing columns: ", paste(miss_anno, collapse = ", "))

  # ---- helper: V from MAF (96 x N) ----
  make_V <- function(maf, order96, Hs) {
    maf <- tibble::as_tibble(maf)
    snv <- dplyr::filter(maf, Variant_Type == "SNP")
    if (nrow(snv) == 0) {
      return(matrix(0, nrow = length(order96), ncol = 0, dimnames = list(order96, NULL)))
    }

    gr5 <- GenomicRanges::GRanges(
      seqnames = paste0("chr", snv$Chromosome),
      IRanges::IRanges(start = snv$Start_Position - 1, end = snv$Start_Position - 1),
      strand = "+"
    )
    gr3 <- GenomicRanges::GRanges(
      seqnames = paste0("chr", snv$Chromosome),
      IRanges::IRanges(start = snv$Start_Position + 1, end = snv$Start_Position + 1),
      strand = "+"
    )

    s5 <- BSgenome::getSeq(Hs, gr5)
    s3 <- BSgenome::getSeq(Hs, gr3)

    ref_tri <- paste0(as.character(s5), snv$Reference_Allele, as.character(s3))
    mut_tri <- paste0(as.character(s5), snv$Tumor_Seq_Allele2, as.character(s3))

    flip <- snv$Reference_Allele %in% c("A","G")
    if (any(flip)) {
      ref_tri[flip] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(ref_tri[flip])))
      mut_tri[flip] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(mut_tri[flip])))
    }

    category <- paste0(substr(ref_tri,1,1),"[",substr(ref_tri,2,2),">",substr(mut_tri,2,2),"]",substr(ref_tri,3,3))
    counts <- tibble::tibble(Tumor_Sample_Barcode = snv$Tumor_Sample_Barcode, category = category) |>
      dplyr::count(Tumor_Sample_Barcode, category, name = "n")

    wide <- tidyr::pivot_wider(counts, names_from = "category", values_from = "n", values_fill = 0)

    V <- t(as.matrix(wide[, -1, drop = FALSE]))
    colnames(V) <- wide$Tumor_Sample_Barcode

    # add any missing of the 96 categories
    miss96 <- !order96 %in% rownames(V)
    if (any(miss96)) {
      add <- matrix(0, nrow = base::sum(miss96), ncol = ncol(V),
                    dimnames = list(order96[miss96], colnames(V)))
      V <- rbind(V, add)
    }

    V[match(order96, rownames(V)), , drop = FALSE]
  }

  # ---- V ----
  V <- make_V(mutations2, Mut_category_order, Hsapiens)

  # ---- Panel context & L ----
  Panel_context <- SATS::GeneratePanelSize(
    genomic_information = genomic_information,
    Class = "SBS",
    SBS_order = "COSMIC",
    ref.genome = genome
  )
  if (!identical(rownames(Panel_context), rownames(V)))
    stop("Panel_context/V rownames mismatch.")

  keep <- clinical_sample$SAMPLE_ID %in% colnames(V)
  if (!any(keep)) stop("No overlap between clinical SAMPLE_ID and V column names.")

  L <- SATS::GenerateLMatrix(
    Panel_context,
    data.frame(
      PATIENT_ID  = clinical_sample$SAMPLE_ID[keep],
      SEQ_ASSAY_ID = clinical_sample$SEQ_ASSAY_ID[keep]
    )
  )
  if (!identical(rownames(L), rownames(V)))
    stop("L/V rownames mismatch.")

  # ---- Refitting per sample ----
  tumor_ID <- colnames(V)
  out <- tibble::tibble(SAMPLE_ID = character(),
                        Signature  = character(),
                        Activity   = double(),
                        Burden     = double())

  for (sid in tumor_ID) {
    ct <- clinical_sample$CANCER_TYPE[match(sid, clinical_sample$SAMPLE_ID)]
    if (is.na(ct)) stop("No CANCER_TYPE for sample: ", sid)

    # map cancer type (optionally include ONCOTREE)
    if (isTRUE(match_on_oncotree) && "ONCOTREE_CODE" %in% names(clinical_sample) && "ONCOTREE_CODE" %in% names(annoFile)) {
      oc  <- clinical_sample$ONCOTREE_CODE[match(sid, clinical_sample$SAMPLE_ID)]
      idx <- which(annoFile$CANCER_TYPE == ct & annoFile$ONCOTREE_CODE == oc)
    } else {
      idx <- which(annoFile$CANCER_TYPE == ct)
    }
    if (!length(idx)) stop("No mapping in annotation file for sample ", sid, " (CANCER_TYPE=", ct, ").")

    CanType <- annoFile$cancerAnalyzed[idx[1]]

    tumor_sigs <- RefTMB$SBS_refSigs$COSMIC[RefTMB$SBS_refSigs$cancerType == CanType]
    if (!length(tumor_sigs)) stop("No COSMIC signatures for cancerType = ", CanType)

    W <- RefTMB$TMB_SBS_v3.4[, tumor_sigs, drop = FALSE]

    estH   <- SATS::EstimateSigActivity(V = V[, sid, drop = FALSE],
                                        L = L[, sid, drop = FALSE],
                                        W = W)
    Burden <- SATS::CalculateSignatureBurdens(L = L[, sid, drop = FALSE],
                                              W = W, H = estH$H)

    out <- dplyr::bind_rows(
      out,
      tibble::tibble(
        SAMPLE_ID = sid,
        Signature = rownames(estH$H),
        Activity  = round(unname(estH$H[, 1]), 3),
        Burden    = round(unname(Burden[, 1]), 2)
      )
    )
  }

  if (isTRUE(save_csv)) {
    utils::write.csv(out, file.path(output_dir, out_file), row.names = FALSE)
    message("Wrote: ", file.path(output_dir, out_file))
  }

  # Return a list so callers can inspect intermediates if needed
  list(H_Burden = out, V = V, L = L, Panel_context = Panel_context)
}

# ---- DBS refitting function ----
run_dbs_refitting <- function(maf_file,
                              genomic_file,
                              clinical_file,
                              output_dir,
                              common_files_dir,
                              genome = c("hg19","hg38"),
                              save_csv = TRUE,
                              out_file = "H_Burden_est_DBS.csv",
                              match_on_oncotree = FALSE) {

  # Add logging to track function execution
  cat("=== DBS Refitting Function Started ===\n")
  cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Input parameters:\n")
  cat("  maf_file:", maf_file, "\n")
  cat("  genomic_file:", genomic_file, "\n") 
  cat("  clinical_file:", clinical_file, "\n")
  cat("  output_dir:", output_dir, "\n")
  cat("  common_files_dir:", common_files_dir, "\n")
  cat("  genome:", genome, "\n")
  cat("  out_file:", out_file, "\n")
  cat("  match_on_oncotree:", match_on_oncotree, "\n")

  genome <- match.arg(genome)
  cat("Validated genome:", genome, "\n")

  # Check inputs exist
  cat("Checking input files exist...\n")
  user_files <- c(maf_file, genomic_file, clinical_file)
  missing_user <- user_files[!file.exists(user_files)]
  if (length(missing_user)) {
    cat("ERROR: Missing user files:\n", paste(missing_user, collapse = "\n"), "\n")
    stop("Missing user files:\n", paste(missing_user, collapse = "\n"))
  }
  cat("All user files found successfully\n")

  if (!dir.exists(output_dir)) {
    cat("ERROR: Output directory does not exist:", output_dir, "\n")
    stop("Output directory does not exist: ", output_dir)
  }
  cat("Output directory exists:", output_dir, "\n")
  
  if (missing(common_files_dir)) {
    cat("ERROR: common_files_dir must be provided\n")
    stop("common_files_dir must be provided (folder with the two shared CSVs).")
  }
  cat("Common files directory:", common_files_dir, "\n")

  # Shared reference files
  cat("Checking shared reference files...\n")
  ref_signature_file     <- file.path(common_files_dir, "DBS_78_TMB_v3.4_20sigs.csv")
  cancer_dictionary_file <- file.path(common_files_dir, "0 Cancer_Dictionary_BZ.csv")
  common_needed <- c(ref_signature_file, cancer_dictionary_file)
  missing_common <- common_needed[!file.exists(common_needed)]
  if (length(missing_common)) {
    cat("ERROR: Missing common reference files:\n", paste(missing_common, collapse = "\n"), "\n")
    stop("Missing common reference files:\n", paste(missing_common, collapse = "\n"))
  }
  cat("All reference files found successfully\n")

  # Pick genome (package must be installed, but we don't attach it globally)
  Hsapiens <- switch(
    genome,
    hg19 = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
    hg38 = BSgenome.Hsapiens.UCSC.hg38::Hsapiens
  )

  # Load RefTMB dataset from SATS (dataset, not an exported object)
  utils::data("RefTMB", package = "SATS", envir = environment())
  if (!exists("RefTMB", envir = environment())) {
    stop("Failed to load RefTMB from SATS. Check SATS installation.")
  }

  # ---- read inputs ----
  DBS_TMBv3.4 <- utils::read.csv(ref_signature_file, row.names = 1, check.names = FALSE)
  Mut_category_order <- rownames(DBS_TMBv3.4)
  if (length(Mut_category_order) != 78) stop("Expected 78 DBS categories in reference signature file.")

  mutations2 <- utils::read.table(maf_file, header = TRUE, sep = "\t", quote = "",
                                  stringsAsFactors = FALSE, check.names = FALSE)
  genomic_information <- utils::read.table(genomic_file, header = TRUE, sep = "\t", quote = "",
                                           stringsAsFactors = FALSE, check.names = FALSE)
  clinical_sample <- utils::read.table(clinical_file, header = TRUE, sep = "\t", quote = "",
                                       stringsAsFactors = FALSE, check.names = FALSE)
  annoFile <- utils::read.csv(cancer_dictionary_file, check.names = FALSE)

  # Column checks
  need_mut <- c("Chromosome","Start_Position","End_Position","Variant_Type","Reference_Allele","Tumor_Seq_Allele2","Tumor_Sample_Barcode")
  miss_mut <- setdiff(need_mut, names(mutations2))
  if (length(miss_mut)) stop("MAF missing columns: ", paste(miss_mut, collapse = ", "))

  need_cli <- c("SAMPLE_ID","SEQ_ASSAY_ID","CANCER_TYPE")
  miss_cli <- setdiff(need_cli, names(clinical_sample))
  if (length(miss_cli)) stop("Clinical file missing columns: ", paste(miss_cli, collapse = ", "))

  need_anno <- c("CANCER_TYPE","cancerAnalyzed")
  miss_anno <- setdiff(need_anno, names(annoFile))
  if (length(miss_anno)) stop("Annotation file missing columns: ", paste(miss_anno, collapse = ", "))

  # ---- helper: V from MAF (78 x N) for DBS ----
  make_DBS_V <- function(maf, order78, Hs) {
    maf <- tibble::as_tibble(maf)
    dnv <- dplyr::filter(maf, Variant_Type == "DNP")
    if (nrow(dnv) == 0) {
      return(matrix(0, nrow = length(order78), ncol = 0, dimnames = list(order78, NULL)))
    }

    dnv_5_GRanges <- GenomicRanges::GRanges(
      seqnames = paste0("chr", dnv$Chromosome),
      IRanges::IRanges(start = dnv$Start_Position, end = dnv$Start_Position),
      strand = "+"
    )
    
    dnv_3_GRanges <- GenomicRanges::GRanges(
      seqnames = paste0("chr", dnv$Chromosome),
      IRanges::IRanges(start = dnv$End_Position, end = dnv$End_Position),
      strand = "+"
    )
    
    DNV_5 <- BSgenome::getSeq(Hs, dnv_5_GRanges)
    DNV_3 <- BSgenome::getSeq(Hs, dnv_3_GRanges)
    
    dnv$ref_di <- paste0(as.character(DNV_5), as.character(DNV_3))
    dnv$mut_di <- dnv$Tumor_Seq_Allele2
    
    dnv$ref_di_10 <- dnv$ref_di
    dnv$mut_di_10 <- dnv$mut_di
    
    # Misspecified sequence check
    Misspecified <- mapply(function(s1, s2) {
      chars1 <- strsplit(s1, "")[[1]]
      chars2 <- strsplit(s2, "")[[1]]
      base::sum(chars1 == chars2)
    }, dnv$ref_di_10, dnv$mut_di_10)
    
    dnv <- dnv[Misspecified == 0, ]
    
    # Flip to canonical representation
    flip_idx <- !dnv$ref_di %in% c("AC","AT","CC","CG","CT","GC","TA","TC","TG","TT")
    if (any(flip_idx)) {
      dnv$ref_di_10[flip_idx] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(dnv$ref_di_10[flip_idx])))
      dnv$mut_di_10[flip_idx] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(dnv$mut_di_10[flip_idx])))
    }
    
    # Palindrome sequences handling
    flip_idx_AT <- (dnv$ref_di %in% "AT") & (!dnv$mut_di_10 %in% c("CA", "CC", "CG", "GA", "GC", "TA"))
    if (any(flip_idx_AT)) {
      dnv$mut_di_10[flip_idx_AT] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(dnv$mut_di_10[flip_idx_AT])))
    }
    
    flip_idx_TA <- (dnv$ref_di %in% "TA") & (!dnv$mut_di_10 %in% c("AT", "CG", "CT", "GC", "GG", "GT"))
    if (any(flip_idx_TA)) {
      dnv$mut_di_10[flip_idx_TA] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(dnv$mut_di_10[flip_idx_TA])))
    }
    
    flip_idx_CG <- (dnv$ref_di %in% "CG") & (!dnv$mut_di_10 %in% c("AT", "GC", "GT", "TA", "TC", "TT"))
    if (any(flip_idx_CG)) {
      dnv$mut_di_10[flip_idx_CG] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(dnv$mut_di_10[flip_idx_CG])))
    }
    
    flip_idx_GC <- (dnv$ref_di %in% "GC") & (!dnv$mut_di_10 %in% c("AA", "AG", "AT", "CA", "CG", "TA"))
    if (any(flip_idx_GC)) {
      dnv$mut_di_10[flip_idx_GC] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(dnv$mut_di_10[flip_idx_GC])))
    }
    
    dnv$category <- paste0(dnv$ref_di_10, dnv$mut_di_10)
    
    counts <- tibble::tibble(Tumor_Sample_Barcode = dnv$Tumor_Sample_Barcode, category = dnv$category) |>
      dplyr::count(Tumor_Sample_Barcode, category, name = "n")

    wide <- tidyr::pivot_wider(counts, names_from = "category", values_from = "n", values_fill = 0)

    V <- t(as.matrix(wide[, -1, drop = FALSE]))
    colnames(V) <- wide$Tumor_Sample_Barcode

    # add any missing of the 78 categories
    miss78 <- !order78 %in% rownames(V)
    if (any(miss78)) {
      add <- matrix(0, nrow = base::sum(miss78), ncol = ncol(V),
                    dimnames = list(order78[miss78], colnames(V)))
      V <- rbind(V, add)
    }

    V[match(order78, rownames(V)), , drop = FALSE]
  }

  # ---- V ----
  V <- make_DBS_V(mutations2, Mut_category_order, Hsapiens)

  # ---- Panel context & L ----
  Panel_context <- SATS::GeneratePanelSize(
    genomic_information = genomic_information,
    Class = "DBS",
    SBS_order = "COSMIC",
    ref.genome = genome
  )
  if (!identical(rownames(Panel_context), rownames(V)))
    stop("Panel_context/V rownames mismatch.")

  keep <- clinical_sample$SAMPLE_ID %in% colnames(V)
  if (!any(keep)) stop("No overlap between clinical SAMPLE_ID and V column names.")

  L <- SATS::GenerateLMatrix(
    Panel_context,
    data.frame(
      PATIENT_ID  = clinical_sample$SAMPLE_ID[keep],
      SEQ_ASSAY_ID = clinical_sample$SEQ_ASSAY_ID[keep]
    )
  )
  if (!identical(rownames(L), rownames(V)))
    stop("L/V rownames mismatch.")

  # ---- Refitting per sample ----
  tumor_ID <- colnames(V)
  out <- tibble::tibble(SAMPLE_ID = character(),
                        Signature  = character(),
                        Activity   = double(),
                        Burden     = double())

  for (sid in tumor_ID) {
    ct <- clinical_sample$CANCER_TYPE[match(sid, clinical_sample$SAMPLE_ID)]
    if (is.na(ct)) stop("No CANCER_TYPE for sample: ", sid)

    # map cancer type (optionally include ONCOTREE)
    if (isTRUE(match_on_oncotree) && "ONCOTREE_CODE" %in% names(clinical_sample) && "ONCOTREE_CODE" %in% names(annoFile)) {
      oc  <- clinical_sample$ONCOTREE_CODE[match(sid, clinical_sample$SAMPLE_ID)]
      idx <- which(annoFile$CANCER_TYPE == ct & annoFile$ONCOTREE_CODE == oc)
    } else {
      idx <- which(annoFile$CANCER_TYPE == ct)
    }
    if (!length(idx)) stop("No mapping in annotation file for sample ", sid, " (CANCER_TYPE=", ct, ").")

    CanType <- annoFile$cancerAnalyzed[idx[1]]

    tumor_sigs <- RefTMB$DBS_refSigs$COSMIC[RefTMB$DBS_refSigs$cancerType == CanType]
    if (!length(tumor_sigs)) stop("No COSMIC signatures for cancerType = ", CanType)

    W <- RefTMB$TMB_DBS_v3.4[, tumor_sigs, drop = FALSE]

    estH   <- SATS::EstimateSigActivity(V = V[, sid, drop = FALSE],
                                        L = L[, sid, drop = FALSE],
                                        W = W)
    Burden <- SATS::CalculateSignatureBurdens(L = L[, sid, drop = FALSE],
                                              W = W, H = estH$H)

    out <- dplyr::bind_rows(
      out,
      tibble::tibble(
        SAMPLE_ID = sid,
        Signature = rownames(estH$H),
        Activity  = round(unname(estH$H[, 1]), 3),
        Burden    = round(unname(Burden[, 1]), 2)
      )
    )
  }

  if (isTRUE(save_csv)) {
    utils::write.csv(out, file.path(output_dir, out_file), row.names = FALSE)
    message("Wrote: ", file.path(output_dir, out_file))
  }

  # Return a list so callers can inspect intermediates if needed
  list(H_Burden = out, V = V, L = L, Panel_context = Panel_context)
}
