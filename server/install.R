# parallelize installation using all available cores
options(Ncpus = parallel::detectCores())
Sys.setenv(MAKEFLAGS = paste0("-j", parallel::detectCores()))

install.packages("remotes")
# Keep installed dependencies unless an explicit version below replaces them.
options(remotes.upgrade = "never")

remotes::install_version("tidyverse", version = "2.0.0")
remotes::install_version("jsonlite", version = "2.0.0")
remotes::install_version("BiocManager", version = "1.30.27")
remotes::install_version("aws.ec2metadata", version = "0.2.0")
remotes::install_version("cowplot", version = "1.2.0")
remotes::install_version("svglite", version = "2.2.2")
remotes::install_version("conflicted", version = "1.2.0")
remotes::install_version("coop", version = "0.6-3")
remotes::install_version("dplyr", version = "1.1.4")
remotes::install_version("entropy", version = "1.3.2")
remotes::install_version("factoextra", version = "1.0.7")
remotes::install_version("ggdendro", version = "0.2.0")
remotes::install_version("ggforce", version = "0.5.0")
remotes::install_version("ggplot2", version = "4.0.4")
remotes::install_version("ggpubr", version = "0.6.2")
remotes::install_version("ggrepel", version = "0.9.8")
remotes::install_version("ggridges", version = "0.5.7")
remotes::install_version("ggsci", version = "4.1.0")
remotes::install_version("ggstatsplot", version = "0.13.3")
remotes::install_version("ggtext", version = "0.1.2")
remotes::install_version("glue", version = "1.8.1")
remotes::install_version("hrbrthemes", version = "0.9.3")
remotes::install_version("janitor", version = "2.2.1")
remotes::install_version("scales", version = "1.4.0")
remotes::install_version("statsExpressions", version = "1.7.1")
remotes::install_version("aws.s3", version = "0.3.22")
remotes::install_version("broom", version = "1.0.13")
remotes::install_version("purrr", version = "1.2.2")
remotes::install_version("plotly", version = "4.11.0")
remotes::install_version("stringi", version = "1.8.9")
remotes::install_version("stringr", version = "1.6.0")
remotes::install_version("igraph", version = "2.2.1")
remotes::install_version("networkD3", version = "0.4.1")

remotes::install_local("./r-packages/skit-0.0.2.tar.gz")

remotes::install_github("YuLab-SMU/treeio")
remotes::install_github("YuLab-SMU/ggtree")

# Resolve Bioconductor dependencies together for the installed R version.
BiocManager::install(c(
  "BSgenome.Hsapiens.UCSC.hg19",
  "BSgenome.Hsapiens.UCSC.hg38",
  "org.Hs.eg.db",
  "ChIPseeker",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "BSgenome",
  "Biostrings",
  "GenomicRanges",
  "TreeAndLeaf"
), update = FALSE, ask = FALSE)
