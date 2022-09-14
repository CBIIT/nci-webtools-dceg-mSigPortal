#!/bin/bash
set -ex

Rscript combineDatasets.R 'Association/*_*_*_vardata.RData' combineAssociationFiles Association/data.csv
Rscript combineDatasets.R 'Exposure/exposure_refdata.RData' combineExposureFiles Exposure/data.csv
Rscript combineDatasets.R 'Seqmatrix/*_*_seqmatrix_refdata.RData' combineSeqmatrixFiles Seqmatrix/data.csv
Rscript combineDatasets.R 'Signature/signature_refsets.RData' combineSignatureFiles Signature/data.csv
Rscript combineDatasets.R 'Others/content_data_all.RData' combinePatternFiles Others/pattern.csv
Rscript combineDatasets.R 'Etiology/aetiology_exposure.RData' combineEtiology Etiology/aetiology_exposure.csv
Rscript combineDatasets.R 'Etiology/aetiology_organ_specific_signature.RData' combineEtiologyOrgan Etiology/aetiology_organ_specific_signature.csv