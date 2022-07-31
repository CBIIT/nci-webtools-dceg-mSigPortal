#!/bin/bash
set -ex

Rscript exportJsonAsCsv.R Etiology/Etiology_cancer_specific_signatures.json Category=CancerSpecificSignatures Etiology/Etiology_cancer_specific_signatures.csv
Rscript exportJsonAsCsv.R Etiology/Etiology_cancer_therapies.json Category=CancerTherapies Etiology/Etiology_cancer_therapies.csv
Rscript exportJsonAsCsv.R Etiology/Etiology_cosmic.json Category=Cosmic Etiology/Etiology_cosmic.csv
Rscript exportJsonAsCsv.R Etiology/Etiology_enviromental_mutagenesis.json Category=EnviromentalMutagenesis Etiology/Etiology_enviromental_mutagenesis.csv
Rscript exportJsonAsCsv.R Etiology/Etiology_gene_edits.json Category=GeneEdits Etiology/Etiology_gene_edits.csv
Rscript exportJsonAsCsv.R Etiology/Etiology_others.json Category=Other Etiology/Etiology_others.csv