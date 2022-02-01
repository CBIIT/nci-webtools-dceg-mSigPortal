
Sys.setenv(
  "R_REMOTES_UPGRADE" = "always",
  "R_REMOTES_NO_ERRORS_FROM_WARNINGS" = TRUE
)

install.packages(
  c(
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
  ),
  repos = 'https://cloud.r-project.org/'
)

BiocManager::install(
  c(
    'BSgenome.Hsapiens.UCSC.hg19',
    'BSgenome.Hsapiens.UCSC.hg38',
    'ChIPseeker',
    'org.Hs.eg.db',
    'TxDb.Hsapiens.UCSC.hg19.knownGene',
    'TxDb.Hsapiens.UCSC.hg38.knownGene'
  )
)