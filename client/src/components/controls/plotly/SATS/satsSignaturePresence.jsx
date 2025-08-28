import { groupBy } from 'lodash';

export default function SATSSignaturePresence(data, options = {}) {
  console.log('SATS function called with data:', data);
  console.log('Data type:', typeof data);
  console.log('Data keys:', data ? Object.keys(data) : 'null/undefined');
  
  if (!data) {
    console.error('SATS: No data provided');
    return { traces: [], layout: {}, config: {} };
  }
  // Define signature color mapping based on the R code
  const signatureColors = {
    'SBS1': '#1F77B4',
    'SBS_Flat': '#AEC7E8', 
    'SBS_Artefactes': '#37474F',
    'SBS2_13': '#FF7F0E',
    'SBS2/13': '#FF7F0E',
    'SBS84': '#FFBB78',
    'SBS7a': '#2CA02C',
    'SBS7b': '#98DF8A',
    'SBS4': '#A50021',
    'SBS25': '#ADFF2F',
    'SBS11': '#9467BD',
    'SBS32': '#F7B6D2',
    'SBS87': '#FF95A8',
    'SBS10a': '#860086',
    'SBS10b': '#BB00BB',
    'SBS10c': '#F100F1',
    'SBS6': '#94A323',
    'SBS14': '#F07E6E',
    'SBS15': '#9EDAE5',
    'SBS30': '#F9D23C',
    'SBS44': '#80CDC1',
    'SBS12': '#AAA488',
    'SBS19': '#F6E8C3',
    'SBS89': '#BE9C2E',
    'SBS94': '#00FFFF',
    'SBS97': '#FF3D3D',
    'SBS8': '#808080',
    'SBS21': '#808080',
    'SBS31': '#808080',
    'SBS37': '#808080',
    'SBS39': '#808080',
    'SBS92': '#808080',
    // Add patterns for complex signature names
    'Flat': '#AEC7E8',
    'Artefactes': '#37474F'
  };
  
  // Function to get color for a signature
  const getSignatureColor = (signature) => {
    // First try exact match
    if (signatureColors[signature]) {
      return signatureColors[signature];
    }
    
    // Extract base signature name from complex descriptions
    if (signature.includes('SBS1(')) return signatureColors['SBS1'];
    if (signature.includes('SBS2/13(') || signature.includes('SBS2_13(')) return signatureColors['SBS2/13'];
    if (signature.includes('SBS4(')) return signatureColors['SBS4'];
    if (signature.includes('SBS6(')) return signatureColors['SBS6'];
    if (signature.includes('SBS7a(')) return signatureColors['SBS7a'];
    if (signature.includes('SBS7b(')) return signatureColors['SBS7b'];
    if (signature.includes('SBS8(')) return signatureColors['SBS8'];
    if (signature.includes('SBS10a(')) return signatureColors['SBS10a'];
    if (signature.includes('SBS10b(')) return signatureColors['SBS10b'];
    if (signature.includes('SBS10c(')) return signatureColors['SBS10c'];
    if (signature.includes('SBS11(')) return signatureColors['SBS11'];
    if (signature.includes('SBS12(')) return signatureColors['SBS12'];
    if (signature.includes('SBS14(')) return signatureColors['SBS14'];
    if (signature.includes('SBS15(')) return signatureColors['SBS15'];
    if (signature.includes('SBS19(')) return signatureColors['SBS19'];
    if (signature.includes('SBS25(')) return signatureColors['SBS25'];
    if (signature.includes('SBS30(')) return signatureColors['SBS30'];
    if (signature.includes('SBS32(')) return signatureColors['SBS32'];
    if (signature.includes('SBS44(')) return signatureColors['SBS44'];
    if (signature.includes('SBS84(')) return signatureColors['SBS84'];
    if (signature.includes('SBS87(')) return signatureColors['SBS87'];
    if (signature.includes('SBS89(')) return signatureColors['SBS89'];
    if (signature.includes('SBS94(')) return signatureColors['SBS94'];
    if (signature.includes('SBS97(')) return signatureColors['SBS97'];
    if (signature.includes('Flat(')) return signatureColors['SBS_Flat'];
    if (signature.includes('Artefactes(')) return signatureColors['SBS_Artefactes'];
    
    // Default color if no match found
    return '#808080';
  };

  // Define signature annotations
  const signatureAnnotations = {
    'SBS1': 'SBS1(Deamination of 5meC)',
    'SBS2_13': 'SBS2/13(APOBEC)',
    'SBS2/13': 'SBS2/13(APOBEC)',
    'SBS4': 'SBS4(Tobacco smoking)',
    'SBS6': 'SBS6(Defective MMR)',
    'SBS7a': 'SBS7a(UV exposure)',
    'SBS7b': 'SBS7b(UV exposure)',
    'SBS8': 'SBS8(Unknown)',
    'SBS10a': 'SBS10a(POLE-exo*)',
    'SBS10b': 'SBS10b(POLE-exo*)',
    'SBS11': 'SBS11(TMZ treatment)',
    'SBS12': 'SBS12(Unknown)',
    'SBS14': 'SBS14(Defective MMR and POLE−exo*)',
    'SBS15': 'SBS15(Defective MMR)',
    'SBS19': 'SBS19(Unknown)',
    'SBS25': 'SBS25(chemotherapy)',
    'SBS30': 'SBS30(Defective BER)',
    'SBS32': 'SBS32(AZA treatment)',
    'SBS84': 'SBS84(AID)',
    'SBS87': 'SBS87(TP treatment)',
    'SBS89': 'SBS89(Unknown)',
    'SBS94': 'SBS94(Unknown)',
    'SBS97': 'SBS97(Unknown)',
    'SBS_Flat': 'Flat(SBS3/5/40a/40b)',
    'SBS_Artefactes': 'Artefactes(SBS50/51/57)'
  };

  // Expect data in format: 
  // { 
  //   tmbData: [{ CancerType, SBS, Count, Proportion, N, Presence, TMB_all, Label, CancerType_num }],
  //   dotData: [{ CancerType, SBS, Count, Proportion, N, Presence, TMB_all, Label, CancerType_num }] 
  // }
  // OR legacy format:
  // { 
  //   tmbData: [{ cancer, signature, tmb, sampleCount }],
  //   dotData: [{ cancer, signature, presence, sampleCount }] 
  // }
  const { tmbData = [], dotData = [] } = data;

  if (tmbData.length === 0) {
    return { traces: [], layout: {}, config: {} };
  }

  // Handle both new and legacy data formats
  const isNewFormat = tmbData.length > 0 && tmbData[0].hasOwnProperty('CancerType');
  
  let cancerOrder, tmbSignatures, dotSignatures;
  
  // Function to sort signatures in desired order: SBS1, SBS2, SBS3... then Flat, Artefactes
  const sortSignatures = (signatures) => {
    return signatures.sort((a, b) => {
      // Extract signature info for comparison
      const getSignatureInfo = (sig) => {
        // Handle SBS signatures (extract number for sorting)
        const sbsMatch = sig.match(/^SBS(\d+)/);
        if (sbsMatch) {
          return { type: 'SBS', number: parseInt(sbsMatch[1]) };
        }
        
        // Handle DBS signatures
        const dbsMatch = sig.match(/^DBS(\d+)/);
        if (dbsMatch) {
          return { type: 'DBS', number: parseInt(dbsMatch[1]) };
        }
        
        // Handle ID signatures  
        const idMatch = sig.match(/^ID(\d+)/);
        if (idMatch) {
          return { type: 'ID', number: parseInt(idMatch[1]) };
        }
        
        // Handle special signatures
        if (sig.includes('Flat')) return { type: 'SPECIAL', order: 1, name: 'Flat' };
        if (sig.includes('Artefactes')) return { type: 'SPECIAL', order: 2, name: 'Artefactes' };
        
        // Default for unknown signatures
        return { type: 'OTHER', name: sig };
      };
      
      const aInfo = getSignatureInfo(a);
      const bInfo = getSignatureInfo(b);
      
      // Sort by type first: SBS < DBS < ID < SPECIAL < OTHER
      const typeOrder = { 'SBS': 1, 'DBS': 2, 'ID': 3, 'SPECIAL': 4, 'OTHER': 5 };
      if (typeOrder[aInfo.type] !== typeOrder[bInfo.type]) {
        return typeOrder[aInfo.type] - typeOrder[bInfo.type];
      }
      
      // Within same type, sort by number for SBS/DBS/ID
      if (aInfo.type === bInfo.type && ['SBS', 'DBS', 'ID'].includes(aInfo.type)) {
        return aInfo.number - bInfo.number;
      }
      
      // For SPECIAL types, sort by predefined order
      if (aInfo.type === 'SPECIAL' && bInfo.type === 'SPECIAL') {
        return aInfo.order - bInfo.order;
      }
      
      // For OTHER types, sort alphabetically
      return a.localeCompare(b);
    });
  };

  if (isNewFormat) {
    // New format with CancerType, SBS, etc.
    
    // Sort cancer types by total count (descending)  
    const cancerTotalCounts = {};
    tmbData.forEach(item => {
      if (!cancerTotalCounts[item.CancerType]) {
        cancerTotalCounts[item.CancerType] = 0;
      }
      cancerTotalCounts[item.CancerType] += item.Count || 0;
    });

    cancerOrder = Object.entries(cancerTotalCounts)
      .sort((a, b) => b[1] - a[1])
      .map(([cancer, _]) => cancer);

    // Get all signatures and sort them properly
    // For TMB bar chart, reverse the order so higher signatures stack on top
    tmbSignatures = sortSignatures([...new Set(tmbData.map(item => item.SBS))]).reverse();
    // For dot plot, keep normal order (SBS1 at bottom, higher numbers above)
    dotSignatures = sortSignatures([...new Set(dotData.map(item => item.SBS))]);
    
    console.log('🎯 Sorted TMB signatures (reversed for stacking):', tmbSignatures);
    console.log('🎯 Sorted dot signatures:', dotSignatures);
  } else {
    // Legacy format
    
    // Sort cancer types by total TMB (descending)  
    const cancerTMBTotals = {};
    tmbData.forEach(item => {
      if (!cancerTMBTotals[item.cancer]) {
        cancerTMBTotals[item.cancer] = 0;
      }
      cancerTMBTotals[item.cancer] += item.tmb || 0;
    });

    cancerOrder = Object.entries(cancerTMBTotals)
      .sort((a, b) => b[1] - a[1])
      .map(([cancer, _]) => cancer);

    // Get all signatures and sort them properly
    // For TMB bar chart, reverse the order so higher signatures stack on top
    tmbSignatures = sortSignatures([...new Set(tmbData.map(item => item.signature))]).reverse();
    // For dot plot, keep normal order (SBS1 at bottom, higher numbers above)
    dotSignatures = sortSignatures([...new Set(dotData.map(item => item.signature))]);
    
    console.log('🎯 Sorted TMB signatures (legacy, reversed for stacking):', tmbSignatures);
    console.log('🎯 Sorted dot signatures (legacy):', dotSignatures);
  }

  const traces = [];

  // Create TMB bar chart traces (for top subplot)
  tmbSignatures.forEach(signature => {
    let signatureData, yValues, customData;
    
    if (isNewFormat) {
      signatureData = tmbData.filter(item => item.SBS === signature);
      yValues = cancerOrder.map(cancer => {
        const item = signatureData.find(d => d.CancerType === cancer);
        return item ? item.Count || 0 : 0;
      });
      customData = cancerOrder.map(cancer => {
        const item = signatureData.find(d => d.CancerType === cancer);
        return { 
          cancer: cancer,
          count: item ? item.Count || 0 : 0,
          proportion: item ? item.Proportion || 0 : 0,
          tmbAll: item ? item.TMB_all || 0 : 0
        };
      });
    } else {
      signatureData = tmbData.filter(item => item.signature === signature);
      yValues = cancerOrder.map(cancer => {
        const item = signatureData.find(d => d.cancer === cancer);
        return item ? item.tmb || 0 : 0;
      });
      customData = cancerOrder.map(cancer => ({ cancer }));
    }

    traces.push({
      type: 'bar',
      name: signatureAnnotations[signature] || signature,
      x: cancerOrder.map((_, i) => i + 1),
      y: yValues,
      xaxis: 'x',
      yaxis: 'y',
      marker: {
        color: getSignatureColor(signature),
        line: {
          width: 0.5,
          color: 'rgba(0,0,0,0.3)'
        }
      },
      showlegend: true,
      hovertemplate: isNewFormat ?
        '<b>Cancer:</b> %{customdata.cancer}<br>' +
        '<b>Signature:</b> ' + (signatureAnnotations[signature] || signature) + '<br>' +
        '<b>Count:</b> %{customdata.count}<br>' +
        '<b>Proportion:</b> %{customdata.proportion:.3f}<br>' +
        '<b>TMB All:</b> %{customdata.tmbAll:.3f}<br>' +
        '<extra></extra>' :
        '<b>Cancer:</b> %{customdata.cancer}<br>' +
        '<b>Signature:</b> ' + (signatureAnnotations[signature] || signature) + '<br>' +
        '<b>TMB:</b> %{y:.3f} mutations/Mb<br>' +
        '<extra></extra>',
      customdata: customData
    });
  });

  // Create dot plot traces (for bottom subplot) 
  dotSignatures.forEach((signature, sigIndex) => {
    let signatureData;
    const xValues = [];
    const yValues = [];
    const sizes = [];
    const customData = [];

    if (isNewFormat) {
      signatureData = dotData.filter(item => item.SBS === signature);
      cancerOrder.forEach((cancer, cancerIndex) => {
        const item = signatureData.find(d => d.CancerType === cancer);
        if (item && item.Presence > 0) {
          xValues.push(cancerIndex + 1);
          yValues.push(sigIndex + 1);
          sizes.push(Math.max(item.Presence * 20, 4));
          customData.push({
            cancer: cancer,
            signature: signatureAnnotations[signature] || signature,
            presence: item.Presence,
            sampleCount: item.N || 0
          });
        }
      });
    } else {
      signatureData = dotData.filter(item => item.signature === signature);
      cancerOrder.forEach((cancer, cancerIndex) => {
        const item = signatureData.find(d => d.cancer === cancer);
        if (item && item.presence > 0) {
          xValues.push(cancerIndex + 1);
          yValues.push(sigIndex + 1);
          sizes.push(Math.max(item.presence * 20, 4));
          customData.push({
            cancer: cancer,
            signature: signatureAnnotations[signature] || signature,
            presence: item.presence,
            sampleCount: item.sampleCount || 0
          });
        }
      });
    }

    if (xValues.length > 0) {
      traces.push({
        type: 'scatter',
        mode: 'markers',
        name: signatureAnnotations[signature] || signature,
        x: xValues,
        y: yValues,
        xaxis: 'x2',
        yaxis: 'y2',
        marker: {
          size: sizes,
          color: getSignatureColor(signature),
          line: {
            width: 1,
            color: 'rgba(0,0,0,0.3)'
          },
          opacity: 0.8
        },
        showlegend: false,
        hovertemplate: 
          '<b>Cancer:</b> %{customdata.cancer}<br>' +
          '<b>Signature:</b> %{customdata.signature}<br>' +
          '<b>Presence:</b> %{customdata.presence:.1%}<br>' +
          '<b>Sample Count:</b> %{customdata.sampleCount}<br>' +
          '<extra></extra>',
        customdata: customData
      });
    }
  });

  // Create sample count annotations for display between bar chart and dot plot
  const sampleCountAnnotations = cancerOrder.map((cancer, index) => {
    // Get sample count for this cancer type
    let sampleCount = 0;
    
    if (isNewFormat) {
      // Find any item with this cancer type to get the sample count (N field)
      const cancerItem = tmbData.find(item => item.CancerType === cancer);
      sampleCount = cancerItem ? cancerItem.N || 0 : 0;
    } else {
      // Legacy format - get from sampleCount field
      const cancerItem = tmbData.find(item => item.cancer === cancer);
      sampleCount = cancerItem ? cancerItem.sampleCount || 0 : 0;
    }

    return {
      x: index + 1,
      y: 0.77, // Position between bar chart (0.8-1.0) and dot plot (0-0.75)
      xref: 'x',
      yref: 'paper',
      text: `${sampleCount}`,
      showarrow: false,
      font: {
        size: 10,
        color: '#666666'
      },
      xanchor: 'center',
      yanchor: 'middle'
    };
  });

  const layout = {
    title: {
      text: '<b>Mutational Signature Analysis Across Cancer Types</b>',
      font: {
        family: 'Times New Roman',
        size: 18
      },
      x: 0.5
    },
    
    // TMB bar chart (top)
    xaxis: {
      domain: [0, 1],
      anchor: 'y',
      showticklabels: false,
      showgrid: false,
      range: [0.5, cancerOrder.length + 0.5],
      tickmode: 'array',
      tickvals: cancerOrder.map((_, i) => i + 1)
    },
    yaxis: {
      domain: [0.8, 1],
      anchor: 'x',
      title: {
        text: '<b>Total Mutation Burden (TMB)</b>',
        font: {
          family: 'Times New Roman',
          size: 12
        }
      },
      showgrid: true,
      gridcolor: 'rgba(128,128,128,0.2)'
    },
    
    // Dot plot (bottom)
    xaxis2: {
      domain: [0, 1],
      anchor: 'y2',
      tickmode: 'array',
      tickvals: cancerOrder.map((_, i) => i + 1),
      ticktext: cancerOrder,
      tickangle: -45, // Diagonal to the left
      tickfont: {
        size: 10
      },
      range: [0.5, cancerOrder.length + 0.5],
      showgrid: true,
      gridcolor: 'rgba(128,128,128,0.2)'
    },
    yaxis2: {
      domain: [0, 0.75],
      anchor: 'x2',
      title: {
        text: '<b>Mutational Signature</b>',
        font: {
          family: 'Times New Roman',
          size: 14
        }
      },
      tickmode: 'array',
      tickvals: dotSignatures.map((_, i) => i + 1),
      ticktext: dotSignatures.map(sig => signatureAnnotations[sig] || sig),
      tickfont: {
        size: 9
      },
      showgrid: true,
      gridcolor: 'rgba(128,128,128,0.2)',
      autorange: 'reversed'
    },
    
    autosize: true,
    height: 900,
    barmode: 'stack',
    margin: {
      l: 350,
      r: 50,
      t: 80,
      b: 200
    },
    plot_bgcolor: 'white',
    paper_bgcolor: 'white',
    annotations: sampleCountAnnotations
  };

  // Add size legend traces for "Proportion of tumors with the signature"
  // Create invisible traces that will show in legend to explain dot sizes
  const sizeLegendSizes = [10, 20, 30, 40]; // Different sizes
  const sizeLegendLabels = ['0.3', '0.5', '0.7', '0.9']; // Corresponding decimal proportions
  
  sizeLegendSizes.forEach((size, index) => {
    traces.push({
      type: 'scatter',
      mode: 'markers',
      x: [null], // No actual data points
      y: [null],
      marker: {
        size: size,
        color: 'rgba(128,128,128,0.7)', // Gray color for legend
        line: {
          width: 1,
          color: 'rgba(0,0,0,0.3)'
        }
      },
      name: sizeLegendLabels[index],
      legendgroup: 'proportion',
      legendgrouptitle: {
        text: 'Proportion of tumors with the signature'
      },
      showlegend: true,
      hoverinfo: 'skip'
    });
  });

  const config = {
    displayModeBar: true,
    responsive: true,
    displaylogo: false,
    toImageButtonOptions: {
      format: 'svg',
      filename: 'SATS_Combined_Plot',
      width: 1600,
      height: 900
    },
    modeBarButtonsToRemove: ['pan2d', 'lasso2d', 'select2d']
  };

  return { 
    traces, 
    layout, 
    config 
  };
}
