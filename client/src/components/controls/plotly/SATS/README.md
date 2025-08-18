# SATS Plotly Components

This folder contains Plotly.js components for creating SATS (Signature Activity in Tumor Samples) visualizations. These components replicate the R ggplot2 visualizations from the SATS analysis.

## Components

### 1. SATSSignaturePresence (`satsSignaturePresence.jsx`)
Creates a combined plot with:
- **Top panel**: Stacked bar chart showing TMB (Tumor Mutational Burden) per megabase for each cancer type
- **Bottom panel**: Dot plot showing signature presence across cancer types

**Data Format:**
```javascript
{
  tmbData: [
    { cancer: 'BRCA', signature: 'SBS1', tmb: 0.5, sampleCount: 100 },
    // ... more TMB data
  ],
  dotData: [
    { cancer: 'BRCA', signature: 'SBS1', presence: 0.8, sampleCount: 100 },
    // ... more presence data
  ]
}
```

### 2. SATSDotPlot (`satsDotPlot.jsx`)
Creates a standalone dot plot showing signature presence across cancer types.

**Data Format:**
```javascript
[
  { cancer: 'BRCA', signature: 'SBS1', presence: 0.8, sampleCount: 100 },
  // ... more data
]
```

### 3. SATSTMBBar (`satsTMBBar.jsx`)
Creates a standalone stacked bar chart showing TMB by cancer type.

**Data Format:**
```javascript
[
  { cancer: 'BRCA', signature: 'SBS1', tmb: 0.5 },
  // ... more data
]
```

## Usage

```javascript
import { SATSSignaturePresence, SATSDotPlot, SATSTMBBar } from '../plotly/SATS';

// For combined plot
const combinedPlot = SATSSignaturePresence(data, options);

// For individual plots
const dotPlot = SATSDotPlot(data, options);
const barChart = SATSTMBBar(data, options);
```

## Features

- **Color coding**: Signature-specific colors based on the R implementation
- **Interactive tooltips**: Hover information showing cancer type, signature name, TMB/presence values
- **Responsive design**: Plots adapt to container size
- **Export capability**: SVG export functionality
- **Signature annotations**: Full signature names with biological context

## Signature Color Scheme

The components use the same color scheme as the R implementation:
- SBS1 (Deamination): Blue (#1F77B4)
- SBS2/13 (APOBEC): Orange (#FF7F0E)
- SBS4 (Tobacco): Dark red (#A50021)
- POLE signatures: Purple variants
- MMR deficiency: Green variants
- Unknown signatures: Gray (#808080)
- And more...

## Data Requirements

- **cancer**: Cancer type abbreviation (e.g., 'BRCA', 'LUAD')
- **signature**: Signature identifier (e.g., 'SBS1', 'SBS2_13')
- **tmb**: Tumor Mutational Burden value (mutations per megabase)
- **presence**: Proportion of samples with the signature (0-1)
- **sampleCount**: Number of samples (optional, for tooltips)

## Customization

The components accept an optional `options` parameter for customization:
- Layout modifications
- Color overrides
- Size adjustments
- Title changes
