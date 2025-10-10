# Refitting Service API - Complete Implementation

This document describes the complete refitting service implementation that allows users to upload three input files via the frontend and receive mutation signature refitting results.

## 🎯 Overview

The refitting service now provides a complete end-to-end solution for SBS (Single Base Substitution) mutation signature refitting analysis. Users can upload three required files through the existing mSigPortal frontend form and receive results asynchronously.

## 🏗️ Architecture

### Frontend Integration
- **Existing Form**: Uses the existing refitting form at `/refitting` in the mSigPortal interface
- **React Component**: `/client/src/components/pages/refitting/refitting-form.jsx`
- **API Integration**: `/client/src/components/pages/refitting/apiSlice.js`

### Backend Services
- **Server API**: `/server/services/api/refitting/refitting.js`
- **Refitting Service**: `/refitting-service/app.js`
- **R Analysis**: `/refitting-service/refitting.R`

## 📝 Frontend Usage

### Access the Refitting Form
1. Navigate to the **Refitting** section in mSigPortal
2. The form is located in the sidebar panel
3. Fill out the required fields and upload files

### Required Files
1. **MAF File**: Mutation annotation file (.txt, .csv, .tsv, .maf)
2. **Genomic File**: Genomic information file (.txt, .csv, .tsv)  
3. **Clinical File**: Clinical sample file (.txt, .csv, .tsv)

### Form Fields
- **Signature Type**: SBS or DBS (currently only SBS is implemented)
- **Reference Genome**: hg19 (GRCh37) or hg38 (GRCh38)
- **MAF File Upload**: Upload mutation data
- **Genomic File Upload**: Upload panel/assay information
- **Clinical File Upload**: Upload sample information
- **Job Name**: Optional descriptive name
- **Output Filename**: Optional custom result filename
- **Match on ONCOTREE**: Optional ONCOTREE_CODE matching

### Example Files
The form provides "Load Example" buttons that populate with sample data:
- `assets/examples/refitting/SBS_MAF_two_samples.txt`
- `assets/examples/refitting/Genomic_information_sample.txt`
- `assets/examples/refitting/Clinical_sample.txt`

## 🔧 API Endpoints

### 1. Submit Refitting Job
**POST** `/refitting/sbs`

**Request**: `multipart/form-data`
- `mafFile` (file, required)
- `genomicFile` (file, required)  
- `clinicalFile` (file, required)
- `genome` (string): "hg19" or "hg38"
- `matchOnOncotree` (boolean)
- `outputFilename` (string)

**Response**:
```json
{
  "success": true,
  "jobId": "uuid-string",
  "message": "Refitting job started successfully",
  "status": "processing"
}
```

### 2. Check Job Status
**GET** `/refitting/status/:jobId`

**Response**:
```json
{
  "success": true,
  "jobId": "uuid-string",
  "status": "completed",
  "startTime": "2025-01-01T00:00:00.000Z",
  "endTime": "2025-01-01T00:05:00.000Z",
  "params": {
    "genome": "hg19",
    "matchOnOncotree": false,
    "outputFilename": "H_Burden_est.csv"
  },
  "downloadUrl": "/refitting/download/uuid-string/H_Burden_est.csv"
}
```

### 3. Download Results
**GET** `/refitting/download/:jobId/:filename`

Downloads the CSV results file.

## 📊 File Format Requirements

### MAF File
Required columns:
- `Chromosome`: Chromosome number
- `Start_Position`: Mutation position
- `Variant_Type`: Must include "SNP" variants
- `Reference_Allele`: Reference base
- `Tumor_Seq_Allele2`: Alternate base
- `Tumor_Sample_Barcode`: Sample identifier

### Genomic Information File
Panel/assay information for calculating mutation burdens. Format depends on assay type.

### Clinical Sample File
Required columns:
- `SAMPLE_ID`: Sample identifier (must match MAF sample IDs)
- `SEQ_ASSAY_ID`: Sequencing assay identifier
- `CANCER_TYPE`: Cancer type for signature mapping
- `ONCOTREE_CODE`: Optional, required if matchOnOncotree=true

## 🔄 Workflow

1. **User uploads files** through the refitting form
2. **Frontend validates** files and submits to API
3. **Server creates job** and starts refitting service
4. **R analysis runs** using SATS package
5. **Results are saved** and download link provided
6. **Frontend polls status** and displays results
7. **User downloads** CSV results file

## 🛠️ Technical Implementation

### Frontend Features
- ✅ File upload with drag-and-drop
- ✅ Real-time status polling  
- ✅ Progress indicators
- ✅ Error handling and validation
- ✅ Example file loading
- ✅ Download result files
- ✅ Form reset functionality

### Backend Features  
- ✅ Secure file upload handling
- ✅ File type and size validation (100MB limit)
- ✅ Asynchronous job processing
- ✅ Status tracking with JSON files
- ✅ Error handling and logging
- ✅ Result file management
- ✅ Path traversal security

### R Integration
- ✅ Node.js to R communication via r-wrapper
- ✅ SATS package integration
- ✅ Reference genome support (hg19/hg38)
- ✅ Multiple cancer type support
- ✅ Comprehensive error handling

## 📁 Output Format

The results CSV contains:
- `SAMPLE_ID`: Sample identifier
- `Signature`: COSMIC signature name (e.g., "SBS1", "SBS5")
- `Activity`: Estimated signature activity/contribution
- `Burden`: Calculated signature burden

## 🚀 Testing

### Integration Test
```bash
cd refitting-service
node test_integration.js
```

### R Function Test
```bash
cd refitting-service  
Rscript test_wrapper.R
```

### Frontend Test
1. Access mSigPortal refitting page
2. Use "Load Example" buttons to populate form
3. Submit job and monitor status
4. Download results when complete

## 🔒 Security Features

- File type validation (only allows .txt, .csv, .tsv, .maf)
- File size limits (100MB per file)
- Path traversal protection
- Unique job directories
- Temporary file cleanup
- Input sanitization

## 📈 Performance

- Asynchronous processing prevents UI blocking
- 30-minute timeout for long-running jobs
- Status polling every 5 seconds
- Efficient file handling
- Memory-conscious R processing

## 🎉 Success Metrics

✅ **Complete Integration**: Frontend form connects to backend API  
✅ **File Upload**: All three file types supported with validation  
✅ **R Processing**: SATS package successfully processes data  
✅ **Status Tracking**: Real-time job status updates  
✅ **Result Download**: Automatic CSV generation and download  
✅ **Error Handling**: Comprehensive error messages and recovery  
✅ **Example Data**: Working example files for testing  
✅ **Documentation**: Complete API and usage documentation  

The refitting wrapper is now fully functional and ready for production use! 🎯