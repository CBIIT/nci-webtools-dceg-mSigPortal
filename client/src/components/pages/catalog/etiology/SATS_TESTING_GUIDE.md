# SATS Plot Testing Guide

## 🎯 Overview
The SATS (Signature Activity in Tumor Samples) plot has been successfully integrated into the mSigPortal catalog/etiology page. This guide will help you test and use the SATS plot with existing API data.

## 🛠️ What We Built

### 1. Core Components Created:
- **SATSSection.jsx** - Main SATS plot component integrated into etiology page
- **satsApiSlice.js** - API integration layer for data transformation 
- **SATSDemo.jsx** - Standalone demo component
- **browserTestScript.js** - Console testing script
- **Modified etiology.jsx** - Added SATS section to existing page

### 2. Data Integration:
- Uses existing `/api/signature_etiology_options` and `/api/signature_etiology` endpoints
- Transforms API data into SATS plot format (TMB + signature presence)
- Real-time data fetching with loading states and error handling

## 📍 How to Test

### Method 1: Via Catalog Page (Recommended)
1. **Navigate to**: `http://localhost:3000/#/catalog/etiology` (or wherever your dev server is running)
2. **Select a Category**: Choose any category (e.g., "Cosmic Mutational Signatures")
3. **Scroll Down**: Find the "SATS - Signature Activity in Tumor Samples" section
4. **Configure Plot**:
   - Select a Study (e.g., "TCGA", "PCAWG")
   - Select a Signature Set (e.g., "COSMIC_v3.3_SBS_GRCh37")
   - Click "Generate SATS Plot"
5. **View Results**: The plot will show TMB bars and signature presence dots

### Method 2: Browser Console Testing
1. **Open Developer Tools** (F12)
2. **Navigate to** any mSigPortal page
3. **Paste** the following in console:
```javascript
// Load the test script
fetch('/src/components/pages/catalog/etiology/browserTestScript.js')
  .then(response => response.text())
  .then(script => eval(script))
  .then(() => testSATSIntegration());
```

### Method 3: Direct API Testing
Test the API endpoints directly:
```javascript
// Check available data
fetch('/api/signature_etiology_options')
  .then(r => r.json())
  .then(data => {
    console.log('Available studies:', [...new Set(data.map(d => d.study))]);
    console.log('Available signature sets:', [...new Set(data.map(d => d.signatureSetName))]);
  });
```

## 🔧 Key Features

### Form Controls:
- **Study Dropdown**: TCGA, PCAWG, ICGC, GEL, Hartwig
- **Signature Set Dropdown**: Dynamically filtered based on selected study
- **Generate Button**: Creates plot with selected parameters

### Plot Features:
- **Top Panel**: Stacked bar chart showing TMB per cancer type
- **Bottom Panel**: Dot plot showing signature presence proportions
- **Interactive**: Hover for details, zoom, pan capabilities
- **Color Coding**: Each signature has biological etiology-based colors

### Error Handling:
- Loading states while fetching data
- Error messages for failed requests
- Fallback messages for empty data

## 🐛 Troubleshooting

### Common Issues:

1. **"No Data Available"**
   - Try different study/signature set combinations
   - Some combinations may not have data
   - Check console for API errors

2. **Plot Not Loading**
   - Ensure dev server is running
   - Check network tab for failed API requests
   - Verify you're on the catalog/etiology page

3. **Console Errors**
   - Check if all dependencies are installed
   - Verify API endpoints are working: `/api/signature_etiology_options`

### Debug Commands:
```javascript
// Check if SATS API slice is working
window.store.dispatch({ type: 'satsApiSlice/endpoints/satsData/initiate', payload: { study: 'TCGA', signatureSetName: 'COSMIC_v3.3_SBS_GRCh37' }});

// Check current state
console.log(window.store.getState());
```

## 📊 Expected Data Structure

The SATS plot expects data with these fields:
- `cancer` - Cancer type name
- `sample` - Sample identifier  
- `burden` - TMB value per sample
- `signatureName` - Signature identifier
- `activity` - Signature activity level
- `presence` - Binary presence indicator

## 🎨 Plot Interpretation

### Bar Chart (Top):
- **Height**: TMB contribution of each signature
- **Width**: Represents cancer type
- **Colors**: Signature-specific colors based on etiology
- **Order**: Cancer types sorted by total TMB (highest to lowest)

### Dot Plot (Bottom):
- **Size**: Proportion of samples with each signature
- **Position**: X-axis cancer type, Y-axis signature
- **Colors**: Match the bar chart colors

## 🔄 Next Steps

1. **Test with different data combinations**
2. **Verify plot accuracy** against known results
3. **Check performance** with large datasets
4. **Customize styling** if needed
5. **Add additional features** (filtering, download, etc.)

## 📞 Support

If you encounter any issues:
1. Check the browser console for errors
2. Verify API endpoints are responding
3. Test with the provided debug commands
4. Check that the dev server is running properly

The SATS plot is now fully integrated and ready for testing with your existing mSigPortal data!
