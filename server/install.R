install.packages("devtools")
devtools::install_version("tidyverse", version = "2.0.0")
devtools::install_version("jsonlite", version = "2.0.0")
devtools::install_version("BiocManager", version = "1.30.27")
devtools::install_version("PMCMRplus", version = "1.9.12")
devtools::install_version("aws.ec2metadata", version = "0.2.0")
devtools::install_version("cowplot", version = "1.2.0")
devtools::install_version("ggExtra", version = "0.11.0")
devtools::install_version("ggside", version = "0.4.0")
devtools::install_version("svglite", version = "2.2.2")
devtools::install_version("renv", version = "1.1.5")
devtools::install_version("conflicted", version = "1.2.0")
devtools::install_version("coop", version = "0.6-3")
devtools::install_version("dplyr", version = "1.1.4")
devtools::install_version("entropy", version = "1.3.2")
devtools::install_version("factoextra", version = "1.0.7")
devtools::install_version("ggdendro", version = "0.2.0")
devtools::install_version("ggforce", version = "0.5.0")
devtools::install_version("ggplot2", version = "4.0.1")
devtools::install_version("ggpubr", version = "0.6.2")
devtools::install_version("ggrepel", version = "0.9.6")
devtools::install_version("ggridges", version = "0.5.7")
devtools::install_version("ggsci", version = "4.1.0")
devtools::install_version("ggstatsplot", version = "0.13.3")
devtools::install_version("ggtext", version = "0.1.2")
devtools::install_version("glue", version = "1.8.0")
devtools::install_version("hrbrthemes", version = "0.8.7")
devtools::install_version("janitor", version = "2.2.1")
devtools::install_version("scales", version = "1.4.0")
devtools::install_version("statsExpressions", version = "1.7.1")
devtools::install_version("aws.s3", version = "0.3.22")
devtools::install_version("broom", version = "1.0.10")
devtools::install_version("purrr", version = "1.2.0")
devtools::install_version("plotly", version = "4.11.0")
devtools::install_version("stringi", version = "1.8.7")
devtools::install_version("stringr", version = "1.6.0")
devtools::install_version("igraph", version = "2.2.1")
devtools::install_version("networkD3", version = "0.4.1")

devtools::install_local("./r-packages/skit-0.0.2.tar.gz")

BiocManager::install(c(
  "BSgenome.Hsapiens.UCSC.hg19",
  "BSgenome.Hsapiens.UCSC.hg38",
  "ChIPseeker",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "BSgenome",
  "Biostrings",
  "GenomicRanges",
  "TreeAndLeaf"
))
