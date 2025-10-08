# Load required libraries
library(tidyverse)
library(GenomicRanges)
library(Biostrings)
library(SATS)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)

sum <- function(a,b) {
   return(a + b)
 
 }

run_sbs_refitting <- function(maf_file,
                              genomic_file,
                              clinical_file,
                              output_dir,
                              common_files_dir = "/app/refitting-service/data",
                              genome = c("hg19","hg38"),
                              save_csv = TRUE,
                              out_file = "H_Burden_est.csv",
                              match_on_oncotree = FALSE) {

  genome <- match.arg(genome)
  
  # Define paths to common reference files
  ref_signature_file <- file.path(common_files_dir, "Alex_Sigs_TMB_check_V3.4_SBS2_13 together.csv")
  cancer_dictionary_file <- file.path(common_files_dir, "0 Cancer_Dictionary_BZ.csv")
  
  # Check that user-uploaded files exist
  user_files <- c(maf_file, genomic_file, clinical_file)
  missing_user_files <- user_files[!file.exists(user_files)]
  if (length(missing_user_files)) stop("Missing user files: ", paste(missing_user_files, collapse = ", "))
  
  # Check that common reference files exist
  common_files <- c(ref_signature_file, cancer_dictionary_file)
  missing_common_files <- common_files[!file.exists(common_files)]
  if (length(missing_common_files)) stop("Missing common reference files: ", paste(missing_common_files, collapse = ", "))
  
  # Check that output directory exists
  if (!dir.exists(output_dir)) stop("Output directory does not exist: ", output_dir)

  # Get appropriate genome reference
  if (genome == "hg19") {
    Hsapiens <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
  } else {
    Hsapiens <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
  }

  utils::data("RefTMB", package = "SATS", envir = environment())

  # --- read inputs ---
  SBS_TMBv3.4 <- utils::read.csv(ref_signature_file, row.names = 1, check.names = FALSE)
  Mut_category_order <- row.names(SBS_TMBv3.4)
  if (length(Mut_category_order) != 96) stop("Expected 96 SBS categories in reference.")

  mutations2 <- utils::read.table(maf_file, header = TRUE, sep = "\t",
                                  quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  genomic_information <- utils::read.table(genomic_file, header = TRUE, sep = "\t",
                                           quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  clinical_sample <- utils::read.table(clinical_file, header = TRUE, sep = "\t",
                                       quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  annoFile <- utils::read.csv(cancer_dictionary_file, check.names = FALSE)

  # --- column checks ---
  need_mut <- c("Chromosome","Start_Position","Variant_Type","Reference_Allele","Tumor_Seq_Allele2","Tumor_Sample_Barcode")
  if (anyNA(match(need_mut, names(mutations2)))) stop("MAF missing: ", paste(setdiff(need_mut, names(mutations2)), collapse=", "))
  need_cli <- c("SAMPLE_ID","SEQ_ASSAY_ID","CANCER_TYPE")
  if (anyNA(match(need_cli, names(clinical_sample)))) stop("Clinical file missing: ", paste(setdiff(need_cli, names(clinical_sample)), collapse=", "))
  need_anno <- c("CANCER_TYPE","cancerAnalyzed")
  if (anyNA(match(need_anno, names(annoFile)))) stop("Annotation file missing: ", paste(setdiff(need_anno, names(annoFile)), collapse=", "))

  # --- helper: V from MAF ---
  make_V <- function(maf, order96, Hs) {
    maf <- tibble::as_tibble(maf)
    snv <- dplyr::filter(maf, .data$Variant_Type == "SNP")
    if (nrow(snv) == 0) {
      V0 <- matrix(0, nrow = length(order96), ncol = 0, dimnames = list(order96, NULL))
      return(V0)
    }

    gr5 <- GenomicRanges::GRanges(seqnames = paste0("chr", snv$Chromosome),
                                  IRanges::IRanges(start = snv$Start_Position - 1,
                                                   end   = snv$Start_Position - 1),
                                  strand = "+")
    gr3 <- GenomicRanges::GRanges(seqnames = paste0("chr", snv$Chromosome),
                                  IRanges::IRanges(start = snv$Start_Position + 1,
                                                   end   = snv$Start_Position + 1),
                                  strand = "+")

    s5 <- BSgenome::getSeq(Hs, gr5)
    s3 <- BSgenome::getSeq(Hs, gr3)

    ref_tri <- paste0(as.character(s5), snv$Reference_Allele, as.character(s3))
    mut_tri <- paste0(as.character(s5), snv$Tumor_Seq_Allele2, as.character(s3))

    flip <- snv$Reference_Allele %in% c("A","G")
    if (any(flip)) {
      ref_tri[flip] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(ref_tri[flip])))
      mut_tri[flip] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(mut_tri[flip])))
    }

    cat <- paste0(substr(ref_tri,1,1),"[",substr(ref_tri,2,2),">",substr(mut_tri,2,2),"]",substr(ref_tri,3,3))
    snv$category <- cat

    counts <- snv %>% dplyr::group_by(.data$Tumor_Sample_Barcode) %>% dplyr::count(.data$category, name = "n")
    wide <- tidyr::pivot_wider(counts, names_from = .data$category, values_from = .data$n, values_fill = 0)

    V <- t(as.matrix(wide[, -1, drop = FALSE]))
    colnames(V) <- wide$Tumor_Sample_Barcode

    # add any missing 96
    miss96 <- !order96 %in% rownames(V)
    if (any(miss96)) {
      add <- matrix(0, nrow = sum(miss96), ncol = ncol(V), dimnames = list(order96[miss96], colnames(V)))
      V <- rbind(V, add)
    }
    V[match(order96, rownames(V)), , drop = FALSE]
  }

  V <- make_V(mutations2, Mut_category_order, Hsapiens)

  # --- Panel/L ---
  Panel_context <- SATS::GeneratePanelSize(genomic_information, Class = "SBS", SBS_order = "COSMIC", ref.genome = genome)
  if (!identical(rownames(Panel_context), rownames(V))) stop("Panel_context/V rownames mismatch.")

  keep <- clinical_sample$SAMPLE_ID %in% colnames(V)
  if (!any(keep)) stop("No overlap between clinical SAMPLE_ID and V colnames.")

  L <- SATS::GenerateLMatrix(
    Panel_context,
    data.frame(PATIENT_ID = clinical_sample$SAMPLE_ID[keep],
               SEQ_ASSAY_ID = clinical_sample$SEQ_ASSAY_ID[keep])
  )
  if (!identical(rownames(L), rownames(V))) stop("L/V rownames mismatch.")

  # --- Refitting ---
  tumor_ID <- colnames(V)
  out <- tibble::tibble(SAMPLE_ID = character(), Signature = character(), Activity = double(), Burden = double())

  for (sid in tumor_ID) {
    ct <- clinical_sample$CANCER_TYPE[match(sid, clinical_sample$SAMPLE_ID)]
    if (is.na(ct)) stop("No CANCER_TYPE for sample: ", sid)

    if (isTRUE(match_on_oncotree) && "ONCOTREE_CODE" %in% names(clinical_sample) && "ONCOTREE_CODE" %in% names(annoFile)) {
      oc <- clinical_sample$ONCOTREE_CODE[match(sid, clinical_sample$SAMPLE_ID)]
      idx <- which(annoFile$CANCER_TYPE == ct & annoFile$ONCOTREE_CODE == oc)
    } else {
      idx <- which(annoFile$CANCER_TYPE == ct)
    }
    if (!length(idx)) stop("No mapping in annotation file for sample ", sid, " (CANCER_TYPE=", ct, ").")

    CanType <- annoFile$cancerAnalyzed[idx[1]]
    tumor_sigs <- SATS::RefTMB$SBS_refSigs$COSMIC[SATS::RefTMB$SBS_refSigs$cancerType == CanType]
    if (!length(tumor_sigs)) stop("No COSMIC signatures for cancerType=", CanType)

    W <- SATS::RefTMB$TMB_SBS_v3.4[, tumor_sigs, drop = FALSE]

    estH <- SATS::EstimateSigActivity(V = V[, sid, drop = FALSE], L = L[, sid, drop = FALSE], W = W)
    Burden <- SATS::CalculateSignatureBurdens(L = L[, sid, drop = FALSE], W = W, H = estH$H)

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

  if (save_csv) {
    utils::write.csv(out, file.path(output_dir, out_file), row.names = FALSE)
    message("Wrote: ", file.path(output_dir, out_file))
  }

  # return intermediates too
  list(H_Burden = out, V = V, L = L, Panel_context = Panel_context)
}
