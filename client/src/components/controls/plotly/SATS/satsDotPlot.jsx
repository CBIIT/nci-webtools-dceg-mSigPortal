import { groupBy } from 'lodash';

export default function SATSCombinedPlot(data, options = {}) {
  // Define signature color mapping based on the R code
  const signatureColors = {
    'SBS1': '#1F77B4',
    'SBS_Flat': '#AEC7E8', 
    'SBS_Artefactes': '#37474F',
    'SBS2_13': '#FF7F0E',
    'SBS2/13': '#FF7F0E', // Alternative naming
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
    'SBS92': '#808080'
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

  // Process the data
  const processedData = data.map(item => ({
    ...item,
    signatureDisplay: signatureAnnotations[item.signature] || item.signature,
    color: signatureColors[item.signature] || '#808080'
  }));

  // Group by cancer type and sort by total TMB
  const cancerGroups = groupBy(processedData, 'cancer');
  const cancerOrder = Object.entries(cancerGroups)
    .map(([cancer, items]) => ({
      cancer,
      totalTMB: items.reduce((sum, item) => sum + (item.tmb || 0), 0),
      sampleCount: items[0]?.sampleCount || 0
    }))
    .sort((a, b) => b.totalTMB - a.totalTMB)
    .map(item => item.cancer);

  // Create traces for the dot plot
  const signatures = [...new Set(processedData.map(item => item.signature))];
  const traces = [];

  signatures.forEach(signature => {
    const signatureData = processedData.filter(item => item.signature === signature);
    const xValues = [];
    const yValues = [];
    const sizes = [];
    const customData = [];
    const textLabels = [];

    cancerOrder.forEach((cancer, cancerIndex) => {
      const item = signatureData.find(d => d.cancer === cancer);
      if (item && item.presence > 0) {
        xValues.push(cancerIndex + 1);
        yValues.push(signatures.indexOf(signature) + 1);
        sizes.push(Math.max(item.presence * 50, 8)); // Scale size for visibility
        customData.push({
          cancer: cancer,
          signature: item.signatureDisplay,
          presence: item.presence,
          tmb: item.tmb || 0,
          sampleCount: item.sampleCount || 0
        });
        textLabels.push(cancer);
      }
    });

    if (xValues.length > 0) {
      traces.push({
        type: 'scatter',
        mode: 'markers',
        name: signatureAnnotations[signature] || signature,
        x: xValues,
        y: yValues,
        marker: {
          size: sizes,
          color: signatureColors[signature] || '#808080',
          line: {
            width: 1,
            color: 'rgba(0,0,0,0.3)'
          },
          opacity: 0.8
        },
        showlegend: false,
        hovertemplate: 
          '<b>Cancer Type:</b> %{customdata.cancer}<br>' +
          '<b>Signature:</b> %{customdata.signature}<br>' +
          '<b>Presence:</b> %{customdata.presence:.1%}<br>' +
          '<b>Sample Count:</b> %{customdata.sampleCount}<br>' +
          '<extra></extra>',
        customdata: customData
      });
    }
  });

  // Create layout
  const layout = {
    title: {
      text: '<b>Mutational Signature Presence Across Cancer Types</b>',
      font: {
        family: 'Times New Roman',
        size: 18
      },
      x: 0.5
    },
    xaxis: {
      title: {
        text: '<b>Cancer Type</b>',
        font: {
          family: 'Times New Roman',
          size: 14
        }
      },
      tickmode: 'array',
      tickvals: cancerOrder.map((_, i) => i + 1),
      ticktext: cancerOrder,
      tickangle: 45,
      tickfont: {
        size: 10
      },
      showgrid: true,
      gridcolor: 'rgba(128,128,128,0.2)',
      range: [0.5, cancerOrder.length + 0.5]
    },
    yaxis: {
      title: {
        text: '<b>Mutational Signature</b>',
        font: {
          family: 'Times New Roman',
          size: 14
        }
      },
      tickmode: 'array',
      tickvals: signatures.map((_, i) => i + 1),
      ticktext: signatures.map(sig => signatureAnnotations[sig] || sig).reverse(),
      tickfont: {
        size: 10
      },
      showgrid: true,
      gridcolor: 'rgba(128,128,128,0.2)',
      autorange: 'reversed'
    },
    autosize: true,
    height: 700,
    margin: {
      l: 350,
      r: 50,
      t: 80,
      b: 150
    },
    plot_bgcolor: 'white',
    paper_bgcolor: 'white',
    annotations: [
      {
        text: 'Circle size represents proportion of tumors with the signature',
        showarrow: false,
        x: 0.5,
        y: -0.25,
        xref: 'paper',
        yref: 'paper',
        font: {
          size: 12,
          color: 'gray'
        }
      }
    ]
  };

  const config = {
    displayModeBar: true,
    responsive: true,
    displaylogo: false,
    toImageButtonOptions: {
      format: 'svg',
      filename: 'SATS_Signature_Presence_DotPlot',
      width: 1400,
      height: 700
    },
    modeBarButtonsToRemove: ['pan2d', 'lasso2d', 'select2d']
  };

  return { 
    traces, 
    layout, 
    config 
  };
}
