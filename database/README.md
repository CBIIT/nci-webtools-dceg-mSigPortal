# mSigPortal Database

Notes:
1. Sequencing Strategy/Dataset are used interchangably (eg: WES, WGS, etc)

## Association

### Files
- *Study*\_vardata.RData
- *Study*\_*Strategy*\_*Cancer*\_vardata.RData (partition of *Study*\_vardata.RData on Cancer_Type column)

### Columns
1. Cancer\_Type
2. Sample
3. icgc\_specimen\_id
4. icgc\_donor\_id
5. data\_source
6. data\_type
7. variable\_name
8. variable\_value
9. variable\_value\_type

### Sample Rows
|Cancer_Type |Sample |icgc_specimen_id |icgc_donor_id |data_source   |data_type          |variable_name            |variable_value                |variable_value_type |
|:-----------|:------|:----------------|:-------------|:-------------|:------------------|:------------------------|:-----------------------------|:-------------------|
|Breast-DCIS |SP2143 |SP2143           |DO1000        |clinical data |clinical variables |histology_tier1          |ECTODERM                      |character           |
|Breast-DCIS |SP2143 |SP2143           |DO1000        |clinical data |clinical variables |histology_tier2          |Breast                        |character           |
|Breast-DCIS |SP2143 |SP2143           |DO1000        |clinical data |clinical variables |histology_tier3          |In situ adenocarcinoma        |character           |
|Breast-DCIS |SP2143 |SP2143           |DO1000        |clinical data |clinical variables |histology_tier4          |Duct micropapillary carcinoma |character           |
|Breast-DCIS |SP2143 |SP2143           |DO1000        |clinical data |clinical variables |tumour_histological_type |Duct micropapillary carcinoma |character           |
|Breast-DCIS |SP2143 |SP2143           |DO1000        |clinical data |clinical variables |tumour_grade             |G2                            |character           |

#### Partitioned JSON Paths
- /Database/Association/*Study*/*Strategy*/*Cancer\_Type*/*Sample*/data.json


## Exposure

### Files
- exposure_refdata.RData
- *Study*\_*Dataset*\_exposure_refdata.RData (partition of exposure_refdata.RData on Study, Dataset column)

### Columns
1. Study
2. Dataset
3. Cancer_Type
4. Organ
5. Sample
6. Signature_set_name
7. Signature_name
8. Exposure

### Sample Rows
|Study     |Dataset |Cancer_Type |Organ  |Sample   |Signature_set_name                            |Signature_name         | Exposure|
|:---------|:-------|:-----------|:------|:--------|:---------------------------------------------|:----------------------|--------:|
|Breast560 |WGS     |Breast      |Breast |PD10010a |Organ-specific_Cancer_Signatures_GRCh37_SBS96 |Breast_A (Breast_MMR1) |   0.0000|
|Breast560 |WGS     |Breast      |Breast |PD10010a |Organ-specific_Cancer_Signatures_GRCh37_SBS96 |Breast_B (Breast_2)    |   0.0000|
|Breast560 |WGS     |Breast      |Breast |PD10010a |Organ-specific_Cancer_Signatures_GRCh37_SBS96 |Breast_C (Breast_13)   |   0.0000|
|Breast560 |WGS     |Breast      |Breast |PD10010a |Organ-specific_Cancer_Signatures_GRCh37_SBS96 |Breast_D (Breast_MMR2) |   0.0000|
|Breast560 |WGS     |Breast      |Breast |PD10010a |Organ-specific_Cancer_Signatures_GRCh37_SBS96 |Breast_E (Breast_8)    |   0.0000|
|Breast560 |WGS     |Breast      |Breast |PD10010a |Organ-specific_Cancer_Signatures_GRCh37_SBS96 |Breast_F (Breast_18)   | 299.5687|

#### Partitioned JSON Paths
- /Database/Exposure/*Study*/*Dataset*/*Cancer\_Type*/*Signature_set_name*/data.json


## Seqmatrix

### Files
- *Study*\_*Strategy*\_*Cancer*\_seqmatrix_refdata.RData 
- *Study*\_*Strategy*\_*Cancer*\_seqmatrix_refdata_info.RData (unused)

### Columns
1. Study
2. Cancer_Type
3. Sample
4. Dataset
5. Profile
6. MutationType
7. Mutations

### Sample Rows
|Study |Cancer_Type     |Sample   |Dataset |Profile |MutationType | Mutations|
|:-----|:---------------|:--------|:-------|:-------|:------------|---------:|
|PCAWG |Biliary-AdenoCA |SP117655 |WGS     |SBS96   |A[C>A]A      |       269|
|PCAWG |Biliary-AdenoCA |SP117556 |WGS     |SBS96   |A[C>A]A      |       114|
|PCAWG |Biliary-AdenoCA |SP117627 |WGS     |SBS96   |A[C>A]A      |       105|
|PCAWG |Biliary-AdenoCA |SP117775 |WGS     |SBS96   |A[C>A]A      |       217|
|PCAWG |Biliary-AdenoCA |SP117332 |WGS     |SBS96   |A[C>A]A      |        52|
|PCAWG |Biliary-AdenoCA |SP117712 |WGS     |SBS96   |A[C>A]A      |       192|


### Partitioned JSON Paths
- /Database/Seqmatrix/*Study*/*Cancer\_Type*/*Sample*/*Dataset*/*Profile*/data.json


## Signature

### Files
- signature_refsets.RData

### Columns
1. Source
2. Profile
3. Signature_set_name
4. Dataset
5. Strand_info
6. Strand
7. Signature_name
8. MutationType
9. Contribution

### Sample Rows
|Source               |Profile |Signature_set_name                      |Dataset |Strand_info |Strand |Signature_name |MutationType | Contribution|
|:--------------------|:-------|:---------------------------------------|:-------|:-----------|:------|:--------------|:------------|------------:|
|Published_signatures |DBS78   |Other_Published_Signatures_GRCh37_DBS78 |WGS     |N           |NA     |DBS_BPA_WGS    |AC>CA        |    0.0192762|
|Published_signatures |DBS78   |Other_Published_Signatures_GRCh37_DBS78 |WGS     |N           |NA     |DBS_BPA_WGS    |AC>CG        |    0.0050342|
|Published_signatures |DBS78   |Other_Published_Signatures_GRCh37_DBS78 |WGS     |N           |NA     |DBS_BPA_WGS    |AC>CT        |    0.0155328|
|Published_signatures |DBS78   |Other_Published_Signatures_GRCh37_DBS78 |WGS     |N           |NA     |DBS_BPA_WGS    |AC>GA        |    0.0049635|
|Published_signatures |DBS78   |Other_Published_Signatures_GRCh37_DBS78 |WGS     |N           |NA     |DBS_BPA_WGS    |AC>GG        |    0.0083474|
|Published_signatures |DBS78   |Other_Published_Signatures_GRCh37_DBS78 |WGS     |N           |NA     |DBS_BPA_WGS    |AC>GT        |    0.0275530|

### Partitioned JSON Paths
- /Database/Signature/*Profile*/data.json
- /Database/Signature/*Profile*/*Signature_set_name*/data.json

