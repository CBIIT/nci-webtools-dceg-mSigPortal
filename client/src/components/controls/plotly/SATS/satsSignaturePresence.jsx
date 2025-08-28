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
    // SBS colors
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
    'SBS17': '#EF4C7D', // Added SBS17 with distinct color
    'SBS21': '#808080',
    'SBS31': '#808080',
    'SBS37': '#808080',
    'SBS39': '#808080',
    'SBS92': '#808080',
    // Add patterns for complex signature names
    'Flat': '#AEC7E8',
    'Artefactes': '#37474F',
    // DBS colors - using the exact order from R code DBS_col_12
    // Each DBS signature gets its corresponding color by position
    'DBS1': '#E64B35FF',  // 1st color
    'DBS2': '#4DBBD5FF',  // 2nd color
    'DBS3': '#00A087FF',  // 3rd color
    'DBS4': '#3C5488FF',  // 4th color
    'DBS5': '#F39B7FFF',  // 5th color
    'DBS6': '#8491B4FF',  // 6th color
    'DBS7': '#91D1C2FF',  // 7th color
    'DBS8': '#DC0000FF',  // 8th color
    'DBS9': '#7E6148FF',  // 9th color
    'DBS10': '#B09C85FF', // 10th color
    'DBS11': '#F6E8C3',   // 11th color
    'DBS12': '#4DBBD5FF',   // 12th color
    'DBS13': '#E64B35FF', // Cycle back to 1st color
    'DBS16': '#7E6148FF', // 6th color
    'DBS17': '#B09C85FF', // 7th color
    'DBS18': '#F6E8C3',   // 8th color
    'DBS19': '#808080', // 9th color
    // Additional DBS colors continue the original cycle
    'DBS20': '#B09C85FF', // 10th color
    'DBS21': '#F6E8C3'    // 11th color
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
    if (signature.includes('SBS17(')) return signatureColors['SBS17'];
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
    
    // DBS color mapping
    if (signature.includes('DBS1(') || signature.startsWith('DBS1')) return signatureColors['DBS1'];
    if (signature.includes('DBS2(') || signature.startsWith('DBS2')) return signatureColors['DBS2'];
    if (signature.includes('DBS3(') || signature.startsWith('DBS3')) return signatureColors['DBS3'];
    if (signature.includes('DBS4(') || signature.startsWith('DBS4')) return signatureColors['DBS4'];
    if (signature.includes('DBS5(') || signature.startsWith('DBS5')) return signatureColors['DBS5'];
    if (signature.includes('DBS6(') || signature.startsWith('DBS6')) return signatureColors['DBS6'];
    if (signature.includes('DBS7(') || signature.startsWith('DBS7')) return signatureColors['DBS7'];
    if (signature.includes('DBS8(') || signature.startsWith('DBS8')) return signatureColors['DBS8'];
    if (signature.includes('DBS9(') || signature.startsWith('DBS9')) return signatureColors['DBS9'];
    if (signature.includes('DBS10(') || signature.startsWith('DBS10')) return signatureColors['DBS10'];
    if (signature.includes('DBS11(') || signature.startsWith('DBS11')) return signatureColors['DBS11'];
    if (signature.includes('DBS12(') || signature.startsWith('DBS12')) return signatureColors['DBS12'];
    if (signature.includes('DBS13(') || signature.startsWith('DBS13')) return signatureColors['DBS13'];
    if (signature.includes('DBS14(') || signature.startsWith('DBS14')) return signatureColors['DBS14'];
    if (signature.includes('DBS15(') || signature.startsWith('DBS15')) return signatureColors['DBS15'];
    if (signature.includes('DBS16(') || signature.startsWith('DBS16')) return signatureColors['DBS16'];
    if (signature.includes('DBS17(') || signature.startsWith('DBS17')) return signatureColors['DBS17'];
    if (signature.includes('DBS18(') || signature.startsWith('DBS18')) return signatureColors['DBS18'];
    if (signature.includes('DBS19(') || signature.startsWith('DBS19')) return signatureColors['DBS19'];
    if (signature.includes('DBS20(') || signature.startsWith('DBS20')) return signatureColors['DBS20'];
    
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
    'SBS44': 'SBS44(Defective MMR)',
    'SBS84': 'SBS84(AID/APOBEC)',
    'SBS87': 'SBS87(Thiopurine treatment)',
    'SBS89': 'SBS89(Unknown)',
    'SBS94': 'SBS94(Unknown)',
    'SBS97': 'SBS97(Unknown)',
    'SBS_Flat': 'Flat(SBS3/5/40a/40b)',
    'SBS_Artefactes': 'Artefactes(SBS50/51/57)',
    // DBS annotations - matching R code
    'DBS1': 'DBS1(UV exposure)',
    'DBS2': 'DBS2(Tobacco smoking or other mutagens)',
    'DBS3': 'DBS3(POLE-exo*)',
    'DBS4': 'DBS4(Unknown)',
    'DBS6': 'DBS6(Unknown)',
    'DBS9': 'DBS9(Unknown)',
    'DBS11': 'DBS11(Unknown)',
    'DBS12': 'DBS12(Unknown)',
    'DBS13': 'DBS13(HR deficiency)',
    'DBS16': 'DBS16(Unknown)',
    'DBS17': 'DBS17(Unknown)',
    'DBS18': 'DBS18(Unknown)',
    'DBS19': 'DBS19(Unknown)',
    // Additional DBS signatures
    'DBS5': 'DBS5(Unknown)',
    'DBS7': 'DBS7(Unknown)',
    'DBS8': 'DBS8(Unknown)',
    'DBS10': 'DBS10(Unknown)',
    'DBS14': 'DBS14(Unknown)',
    'DBS15': 'DBS15(Unknown)',
    'DBS20': 'DBS20(Unknown)'
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
        // Handle SBS signatures with sub-types (e.g., SBS7a, SBS7b, SBS10a, SBS10b, SBS10c)
        const sbsSubMatch = sig.match(/^SBS(\d+)([abc]?)/);
        if (sbsSubMatch) {
          const number = parseInt(sbsSubMatch[1]);
          const subType = sbsSubMatch[2] || ''; // a, b, c, or empty string
          return { type: 'SBS', number: number, subType: subType };
        }
        
        // Handle regular SBS signatures
        const sbsMatch = sig.match(/^SBS(\d+)/);
        if (sbsMatch) {
          return { type: 'SBS', number: parseInt(sbsMatch[1]), subType: '' };
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
        // First sort by number
        if (aInfo.number !== bInfo.number) {
          return aInfo.number - bInfo.number;
        }
        
        // If same number, sort by subType (for SBS signatures like SBS7a, SBS7b)
        if (aInfo.type === 'SBS' && aInfo.subType !== bInfo.subType) {
          // Empty string (no subtype) comes first, then 'a', 'b', 'c'
          const subTypeOrder = { '': 0, 'a': 1, 'b': 2, 'c': 3 };
          const aOrder = subTypeOrder[aInfo.subType] !== undefined ? subTypeOrder[aInfo.subType] : 999;
          const bOrder = subTypeOrder[bInfo.subType] !== undefined ? subTypeOrder[bInfo.subType] : 999;
          return aOrder - bOrder;
        }
        
        return 0; // They are the same
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
    
    // Check if data already has CancerType_num (pre-sorted order)
    const hasCancerTypeNum = tmbData.some(item => item.CancerType_num !== undefined);
    
    if (hasCancerTypeNum) {
      // Use pre-sorted order from CancerType_num
      const cancerTypeData = {};
      tmbData.forEach(item => {
        if (!cancerTypeData[item.CancerType]) {
          cancerTypeData[item.CancerType] = {
            num: item.CancerType_num,
            tmbAll: item.TMB_all || 0
          };
        }
      });
      
      cancerOrder = Object.entries(cancerTypeData)
        .sort((a, b) => a[1].num - b[1].num) // Sort by CancerType_num (ascending)
        .map(([cancer, _]) => cancer);
        
      console.log('🎯 Using pre-sorted cancer order (CancerType_num):', cancerOrder);
      console.log('🎯 Cancer data:', cancerTypeData);
    } else {
      // Fallback: Use TMB_all for sorting
      const cancerTMBTotals = {};
      tmbData.forEach(item => {
        if (!cancerTMBTotals[item.CancerType]) {
          cancerTMBTotals[item.CancerType] = item.TMB_all || 0;
        }
      });

      cancerOrder = Object.entries(cancerTMBTotals)
        .sort((a, b) => b[1] - a[1]) // Sort by TMB_all descending
        .map(([cancer, _]) => cancer);
        
      console.log('🎯 Cancer order (sorted by TMB_all descending):', cancerOrder);
      console.log('🎯 Cancer TMB totals:', cancerTMBTotals);
    }

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
