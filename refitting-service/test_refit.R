# test_refit.R
# run 
# >> Rscript --vanilla test_refit.R

source("refitting.R")

maf_file      <- "/Users/nguyenthun/Documents/mSigportal/Refitting/SBS_refitting_v2/SBS_MAF_two_samples.txt"
genomic_file  <- "/Users/nguyenthun/Documents/mSigportal/Refitting/SBS_refitting_v2/Genomic_information_sample.txt"
clinical_file <- "/Users/nguyenthun/Documents/mSigportal/Refitting/SBS_refitting_v2/Clinical_sample.txt"

# Shared reference files under ./data/
common_dir <- file.path(getwd(), "data")

# ✅ Output goes to ./data/output
out_dir <- file.path(getwd(), "data", "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

res <- run_sbs_refitting(
  maf_file         = maf_file,
  genomic_file     = genomic_file,
  clinical_file    = clinical_file,
  output_dir       = out_dir,
  common_files_dir = common_dir,
  genome           = "hg19",
  save_csv         = TRUE,
  out_file         = "H_Burden_est.csv",
  match_on_oncotree = FALSE
)

print(head(res$H_Burden))
cat("CSV written to:", file.path(out_dir, "H_Burden_est.csv"), "\n")
