
Sys.setenv(
  "R_REMOTES_UPGRADE" = "always",
  "R_REMOTES_NO_ERRORS_FROM_WARNINGS" = TRUE
)

# warning from failing to install are turned into errors, halting the entire install
# options(warn = 2)

tryCatch({
  packages = c(
        'aws.ec2metadata',
        'aws.s3',
        'BiocManager',
        'car',
        'coop',
        'cowplot',
        'entropy',
        'factoextra',
        'FactoMineR',
        'ggplot2',
        'ggdendro',
        'ggExtra',
        'ggforce',
        'ggpubr',
        'ggrepel',
        'ggridges',
        'ggsci',
        'ggside',
        'ggstatsplot',
        'ggtext',
        'hrbrthemes',
        'janitor',
        'scales',
        'svglite',
        'tidyverse',
        'PMCMRplus'
    )

  install.packages(packages, repos = 'https://cloud.r-project.org/')

  # verify packages are installed
  # lapply(packages, require, character.only = TRUE)

  BiocManager::install(
    c(
        'BSgenome.Hsapiens.UCSC.hg19',
        'BSgenome.Hsapiens.UCSC.hg38',
        'ChIPseeker',
        'org.Hs.eg.db',
        'TxDb.Hsapiens.UCSC.hg19.knownGene',
        'TxDb.Hsapiens.UCSC.hg38.knownGene'
    ))
}, warning = function(e) {
  stop(e)
})