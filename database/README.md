# mSigPortal Database

## Getting Started

### Importing Data

1. Ensure databse/config.json exists and specifies a Postgres database configuration and s3 bucket/key prefix
2. If starting from a clean database, run the createDatabase.js script

```
node createDatabase.js --schema schema.js
```

3. Sync the msigportal Database/ folder from the s3 folder to a local folder (eg: data/)

```
aws s3 sync s3://bucket-name/Database/ path-to-local-folder/Database/
```

4. Copy over the scripts in the utils/ folder to the Database folder
5. Run the following scripts to convert RData into CSV for Postgres Database import:

```
sh exportCombinedDatasets.sh
sh exportEtiologyNormalized.sh
sh exportPublications.sh
```

6. To import to postgres database, Execute the startDatabaseImport.js script with the following arguments:

```
node startDatabaseImport.js --schema schema.js --sources sources.js --provider local path/to/data/folder
```

- Optionally, we can execute the database import remotely in S3 (sync your new Database files there first):

```
aws s3 sync path-to-local-folder/Database/ s3://bucket-name/msigportal/Database/
node startDatabaseImport.js --schema schema.js --sources sources.js --provider s3 s3://bucket-name/msigportal/Database
```

7. Copy Database folder to our EFS volume to update Etiology svg thumbnails

Notes:

1. Sequencing Strategy/Dataset are used interchangably (eg: WES, WGS, etc)

### Upload to EFS

1. SSH to the EFS host for the account's EFS volume (usually the docker host)
2. Sync s3 bucket to the EFS folder
```
aws s3 sync s3://bucket-name/msigportal/Database/ /path/to/analysistools_efs/msigportal/Database/ --delete --dryrun (verify paths with dryrun tag first)
```
3. After sync is complete, update the permissions of the folder
```
chmod -R 755 /paht/to/analysistools_efs/msigportal/
chown -R 1000:1000 /path/to/analysistools_efs/msigportal/
```

## Association

### Files

- _Study_\_vardata.RData
- _Study_\__Strategy_\__Cancer_\_vardata.RData (partition of _Study_\_vardata.RData on Cancer_Type column)

### Columns

1. Cancer_Type
2. Sample
3. icgc_specimen_id
4. icgc_donor_id
5. data_source
6. data_type
7. variable_name
8. variable_value
9. variable_value_type

### Sample Rows

| Cancer_Type | Sample | icgc_specimen_id | icgc_donor_id | data_source   | data_type          | variable_name            | variable_value                | variable_value_type |
| :---------- | :----- | :--------------- | :------------ | :------------ | :----------------- | :----------------------- | :---------------------------- | :------------------ |
| Breast-DCIS | SP2143 | SP2143           | DO1000        | clinical data | clinical variables | histology_tier1          | ECTODERM                      | character           |
| Breast-DCIS | SP2143 | SP2143           | DO1000        | clinical data | clinical variables | histology_tier2          | Breast                        | character           |
| Breast-DCIS | SP2143 | SP2143           | DO1000        | clinical data | clinical variables | histology_tier3          | In situ adenocarcinoma        | character           |
| Breast-DCIS | SP2143 | SP2143           | DO1000        | clinical data | clinical variables | histology_tier4          | Duct micropapillary carcinoma | character           |
| Breast-DCIS | SP2143 | SP2143           | DO1000        | clinical data | clinical variables | tumour_histological_type | Duct micropapillary carcinoma | character           |
| Breast-DCIS | SP2143 | SP2143           | DO1000        | clinical data | clinical variables | tumour_grade             | G2                            | character           |

## Exposure

### Files

- exposure_refdata.RData
- _Study_\__Dataset_\_exposure_refdata.RData (partition of exposure_refdata.RData on Study, Dataset column)

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

| Study     | Dataset | Cancer_Type | Organ  | Sample   | Signature_set_name                            | Signature_name         | Exposure |
| :-------- | :------ | :---------- | :----- | :------- | :-------------------------------------------- | :--------------------- | -------: |
| Breast560 | WGS     | Breast      | Breast | PD10010a | Organ-specific_Cancer_Signatures_GRCh37_SBS96 | Breast_A (Breast_MMR1) |   0.0000 |
| Breast560 | WGS     | Breast      | Breast | PD10010a | Organ-specific_Cancer_Signatures_GRCh37_SBS96 | Breast_B (Breast_2)    |   0.0000 |
| Breast560 | WGS     | Breast      | Breast | PD10010a | Organ-specific_Cancer_Signatures_GRCh37_SBS96 | Breast_C (Breast_13)   |   0.0000 |
| Breast560 | WGS     | Breast      | Breast | PD10010a | Organ-specific_Cancer_Signatures_GRCh37_SBS96 | Breast_D (Breast_MMR2) |   0.0000 |
| Breast560 | WGS     | Breast      | Breast | PD10010a | Organ-specific_Cancer_Signatures_GRCh37_SBS96 | Breast_E (Breast_8)    |   0.0000 |
| Breast560 | WGS     | Breast      | Breast | PD10010a | Organ-specific_Cancer_Signatures_GRCh37_SBS96 | Breast_F (Breast_18)   | 299.5687 |

## Seqmatrix

### Files

- _Study_\__Strategy_\__Cancer_\_seqmatrix_refdata.RData
- _Study_\__Strategy_\__Cancer_\_seqmatrix_refdata_info.RData (unused)

### Columns

1. Study
2. Cancer_Type
3. Sample
4. Dataset
5. Profile
6. MutationType
7. Mutations

### Sample Rows

| Study | Cancer_Type     | Sample   | Dataset | Profile | MutationType | Mutations |
| :---- | :-------------- | :------- | :------ | :------ | :----------- | --------: |
| PCAWG | Biliary-AdenoCA | SP117655 | WGS     | SBS96   | A[C>A]A      |       269 |
| PCAWG | Biliary-AdenoCA | SP117556 | WGS     | SBS96   | A[C>A]A      |       114 |
| PCAWG | Biliary-AdenoCA | SP117627 | WGS     | SBS96   | A[C>A]A      |       105 |
| PCAWG | Biliary-AdenoCA | SP117775 | WGS     | SBS96   | A[C>A]A      |       217 |
| PCAWG | Biliary-AdenoCA | SP117332 | WGS     | SBS96   | A[C>A]A      |        52 |
| PCAWG | Biliary-AdenoCA | SP117712 | WGS     | SBS96   | A[C>A]A      |       192 |

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

| Source               | Profile | Signature_set_name                      | Dataset | Strand_info | Strand | Signature_name | MutationType | Contribution |
| :------------------- | :------ | :-------------------------------------- | :------ | :---------- | :----- | :------------- | :----------- | -----------: |
| Published_signatures | DBS78   | Other_Published_Signatures_GRCh37_DBS78 | WGS     | N           | NA     | DBS_BPA_WGS    | AC>CA        |    0.0192762 |
| Published_signatures | DBS78   | Other_Published_Signatures_GRCh37_DBS78 | WGS     | N           | NA     | DBS_BPA_WGS    | AC>CG        |    0.0050342 |
| Published_signatures | DBS78   | Other_Published_Signatures_GRCh37_DBS78 | WGS     | N           | NA     | DBS_BPA_WGS    | AC>CT        |    0.0155328 |
| Published_signatures | DBS78   | Other_Published_Signatures_GRCh37_DBS78 | WGS     | N           | NA     | DBS_BPA_WGS    | AC>GA        |    0.0049635 |
| Published_signatures | DBS78   | Other_Published_Signatures_GRCh37_DBS78 | WGS     | N           | NA     | DBS_BPA_WGS    | AC>GG        |    0.0083474 |
| Published_signatures | DBS78   | Other_Published_Signatures_GRCh37_DBS78 | WGS     | N           | NA     | DBS_BPA_WGS    | AC>GT        |    0.0275530 |
