#!/bin/bash
set -ex

Rscript combineDatasets.R 'Association/*_*_vardata.RData' combineAssociationFiles Association/data.csv
Rscript combineDatasets.R 'Exposure/exposure_refdata.RData' identity Exposure/data.csv
Rscript combineDatasets.R 'Seqmatrix/*_*_seqmatrix_refdata.RData' identity Seqmatrix/data.csv
Rscript combineDatasets.R 'Signature/signature_refsets.RData' identity Signature/data.csv