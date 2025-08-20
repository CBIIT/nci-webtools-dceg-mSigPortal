# SATS Plot - Clean Production Code

## Overview
The SATS (Signature Activities in Tumor Samples) plot functionality has been cleaned up and is now ready for production use. All demo components and debug code have been removed.

## Core Components

### 1. **SATS Plot Generator**
**File**: `client/src/components/controls/plotly/SATS/satsSignaturePresence.jsx`

Main plotting function that generates the SATS visualization:
- **Input**: Data object with `tmbData` and `dotData` arrays
- **Output**: Plotly configuration object with `traces`, `layout`, and `config`
- **Features**:
  - Stacked bar chart (top panel) showing TMB contributions by signature
  - Dot plot (bottom panel) showing signature presence across cancer types
  - Automatic color mapping for different mutational signatures
  - Supports both new format (CancerType/SBS) and legacy format
  - Cancer types sorted by total count (highest to lowest)

### 2. **API Integration**
**File**: `client/src/components/pages/catalog/etiology/satsApiSlice.js`

RTK Query slice for fetching and transforming data:
- **`useSatsDataQuery`**: Fetches etiology data and generates SATS plot
- **`useSatsOptionsQuery`**: Fetches available studies and signature sets
- **`transformEtiologyDataToSATS`**: Converts API data to SATS format

### 3. **UI Component**
**File**: `client/src/components/pages/catalog/etiology/SATSSection.jsx`

React component with form controls and plot display:
- Study selection dropdown
- Signature set selection dropdown
- Plot generation and display
- Loading states and error handling

## Data Format

### Input Format
```javascript
{
  tmbData: [
    {
      CancerType: "Bladder Cancer",
      SBS: "SBS1(Deamination of 5meC)",
      Count: 3865.526338,
      Proportion: 0.08222075,
      N: 3608,
      Presence: 0.454268293,
      TMB_all: 7.985976405,
      Label: "3608",
      CancerType_num: 3
    },
    // ... more entries
  ],
  dotData: [...] // Same format as tmbData
}
```

### Plot Features
- **Stacked Bar Chart**: Shows total mutation burden with signature contributions
- **Dot Plot**: Circle size represents proportion of samples with signature
- **Color Coding**: Each signature has a distinct color
- **Sorting**: Cancer types ordered by total count (descending)
- **Interactive**: Hover tooltips with detailed information
- **Export**: Built-in export functionality to PNG/SVG

## Usage

### Basic Usage
```javascript
import SATSSignaturePresence from './path/to/satsSignaturePresence';

const plotConfig = SATSSignaturePresence({
  tmbData: yourTmbData,
  dotData: yourDotData
});

// Use with react-plotly.js
<Plot
  data={plotConfig.traces}
  layout={plotConfig.layout}
  config={plotConfig.config}
/>
```

### With React Hook Form Integration
```javascript
import SATSSection from './path/to/SATSSection';

// Simply include in your component
<SATSSection />
```

## Files Removed
- `SATSDemo.jsx` - Demo component
- `SATSSectionDebug.jsx` - Debug component  
- `sampleSatsData.json` - Sample data file
- Demo route from `app.jsx`
- All debug console.log statements

## Files Kept (Core Functionality)
- `satsSignaturePresence.jsx` - Core plotting function
- `satsApiSlice.js` - API integration
- `SATSSection.jsx` - UI component
- Integration in `etiology.jsx` catalog page

The code is now clean, production-ready, and contains only the essential SATS plot generation functionality.
