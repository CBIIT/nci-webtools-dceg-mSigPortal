sum <- function(a,b) {
   return(a + b)
 }

 run_sbs_refitting <- function(path, genome = c("hg19","hg38"),
                              save_csv = TRUE, out_file = "H_Burden_est.csv") {
  genome <- match.arg(genome)
  stopifnot(dir.exists(path))

  # Load packages used directly in code
  suppressPackageStartupMessages({
    library(tidyverse)
    library(GenomicRanges)
    library(Biostrings)
    library(SATS)
    if (genome == "hg19") {
      library(BSgenome.Hsapiens.UCSC.hg19); Hsapiens <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
    } else {
      library(BSgenome.Hsapiens.UCSC.hg38); Hsapiens <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
    }
  })
  data("RefTMB", package = "SATS")

  # Input files
  req <- c("Alex_Sigs_TMB_check_V3.4_SBS2_13 together.csv",
           "SBS_MAF_two_samples.txt",
           "Genomic_information_sample.txt",
           "Clinical_sample.txt",
           "0 Cancer_Dictionary_BZ.csv")
  miss <- req[!file.exists(file.path(path, req))]
  if (length(miss)) stop("Missing files: ", paste(miss, collapse=", "))

  # Read data
  SBS_TMBv3.4 <- read.csv(file.path(path, req[1]), row.names = 1, check.names = FALSE)
  Mut_category_order <- row.names(SBS_TMBv3.4)
  if (length(Mut_category_order) != 96) stop("Expected 96 SBS categories.")

  mutations2 <- read.table(file.path(path, req[2]), header=TRUE, sep="\t", quote="",
                           stringsAsFactors=FALSE, check.names=FALSE)
  genomic_information <- read.table(file.path(path, req[3]), header=TRUE, sep="\t", quote="",
                                    stringsAsFactors=FALSE, check.names=FALSE)
  clinical_sample <- read.table(file.path(path, req[4]), header=TRUE, sep="\t", quote="",
                                stringsAsFactors=FALSE, check.names=FALSE)
  annoFile <- read.csv(file.path(path, req[5]), check.names=FALSE)

  # Helper: build V from MAF
  SBS_V_from_MAF <- function(MAF_with_min_info, Mut_category_order, Hsapiens) {
    MAF_with_min_info <- tibble::as_tibble(MAF_with_min_info)
    SNV <- dplyr::filter(MAF_with_min_info, Variant_Type == "SNP")

    SNV_5_GRanges <- GenomicRanges::GRanges(
      seqnames = paste0("chr", SNV$Chromosome),
      IRanges::IRanges(start = SNV$Start_Position-1, end = SNV$Start_Position-1),
      strand = "+"
    )
    SNV_3_GRanges <- GenomicRanges::GRanges(
      seqnames = paste0("chr", SNV$Chromosome),
      IRanges::IRanges(start = SNV$Start_Position+1, end = SNV$Start_Position+1),
      strand = "+"
    )

    SNV_5 <- BSgenome::getSeq(Hsapiens, SNV_5_GRanges)
    SNV_3 <- BSgenome::getSeq(Hsapiens, SNV_3_GRanges)

    SNV$ref_tri <- paste0(as.character(SNV_5), SNV$Reference_Allele, as.character(SNV_3))
    SNV$mut_tri <- paste0(as.character(SNV_5), SNV$Tumor_Seq_Allele2, as.character(SNV_3))

    flip <- SNV$Reference_Allele %in% c("A","G")
    ref_ct <- SNV$ref_tri; mut_ct <- SNV$mut_tri
    if (any(flip)) {
      ref_ct[flip] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(ref_ct[flip])))
      mut_ct[flip] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(mut_ct[flip])))
    }
    SNV$category <- paste0(substr(ref_ct,1,1),"[",substr(ref_ct,2,2),">",substr(mut_ct,2,2),"]",substr(ref_ct,3,3))

    MutationCounts <- SNV %>% dplyr::group_by(Tumor_Sample_Barcode) %>% dplyr::count(category)
    MutationCounts96 <- MutationCounts %>% tidyr::spread(category, n) %>% replace(is.na(.), 0)

    V <- t(as.matrix(MutationCounts96[, -1]))
    colnames(V) <- dplyr::pull(MutationCounts96, Tumor_Sample_Barcode)

    add_chk <- !Mut_category_order %in% rownames(V)
    if (any(add_chk)) {
      add_V <- matrix(0, nrow = sum(add_chk), ncol = ncol(V),
                      dimnames = list(Mut_category_order[add_chk], colnames(V)))
      V <- rbind(V, add_V)
    }
    if (nrow(V) != length(Mut_category_order)) stop("MAF categories don't match 96-category order.")
    V[match(Mut_category_order, rownames(V)), , drop = FALSE]
  }

  # Build V
  V <- SBS_V_from_MAF(mutations2, Mut_category_order, Hsapiens)

  # Panel context & L
  Panel_context <- GeneratePanelSize(genomic_information, Class="SBS", SBS_order="COSMIC", ref.genome = genome)
  if (!identical(rownames(Panel_context), rownames(V))) stop("Panel_context/V rowname mismatch.")
  keep <- clinical_sample$SAMPLE_ID %in% colnames(V)
  if (!any(keep)) stop("No overlap between clinical SAMPLE_ID and V colnames.")
  L <- GenerateLMatrix(Panel_context,
                       data.frame(PATIENT_ID = clinical_sample$SAMPLE_ID[keep],
                                  SEQ_ASSAY_ID = clinical_sample$SEQ_ASSAY_ID[keep]))
  if (!identical(rownames(L), rownames(V))) stop("L/V rowname mismatch.")

  # Refitting per sample
  tumor_ID <- colnames(V)
  out <- tibble::tibble(SAMPLE_ID=character(), Signature=character(), Activity=double(), Burden=double())

  for (sid in tumor_ID) {
    ct <- clinical_sample$CANCER_TYPE[match(sid, clinical_sample$SAMPLE_ID)]
    idx <- which(annoFile$CANCER_TYPE == ct)
    if (!length(idx)) stop("CANCER_TYPE not mapped in annotation file: ", ct)
    CanType <- annoFile$cancerAnalyzed[idx[1]]

    tumor_sigs <- RefTMB$SBS_refSigs$COSMIC[RefTMB$SBS_refSigs$cancerType == CanType]
    if (!length(tumor_sigs)) stop("No COSMIC signatures for cancerType: ", CanType)
    W <- RefTMB$TMB_SBS_v3.4[, tumor_sigs, drop = FALSE]

    estH <- EstimateSigActivity(V = V[, sid, drop = FALSE], L = L[, sid, drop = FALSE], W = W)
    Burden <- CalculateSignatureBurdens(L = L[, sid, drop = FALSE], W = W, H = estH$H)

    out <- dplyr::bind_rows(out, tibble::tibble(
      SAMPLE_ID = sid,
      Signature = rownames(estH$H),
      Activity  = round(unname(estH$H[,1]), 3),
      Burden    = round(unname(Burden[,1]), 2)
    ))
  }

  if (save_csv) {
    utils::write.csv(out, file.path(path, out_file), row.names = FALSE)
    message("Wrote: ", file.path(path, out_file))
  }
  out
}
