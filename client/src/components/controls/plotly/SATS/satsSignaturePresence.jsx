import { groupBy } from 'lodash';

export default function SATSSignaturePresence(data, options = {}) {
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

  // Expect data in format: 
  // { 
  //   tmbData: [{ cancer, signature, tmb, sampleCount }],
  //   dotData: [{ cancer, signature, presence, sampleCount }] 
  // }
  const { tmbData = [], dotData = [] } = data;

  // Sort cancer types by total TMB (descending)  
  const cancerTMBTotals = {};
  tmbData.forEach(item => {
    if (!cancerTMBTotals[item.cancer]) {
      cancerTMBTotals[item.cancer] = 0;
    }
    cancerTMBTotals[item.cancer] += item.tmb || 0;
  });

  const cancerOrder = Object.entries(cancerTMBTotals)
    .sort((a, b) => b[1] - a[1])
    .map(([cancer, _]) => cancer);

  // Get all signatures
  const tmbSignatures = [...new Set(tmbData.map(item => item.signature))];
  const dotSignatures = [...new Set(dotData.map(item => item.signature))];

  const traces = [];

  // Create TMB bar chart traces (for top subplot)
  tmbSignatures.forEach(signature => {
    const signatureData = tmbData.filter(item => item.signature === signature);
    const yValues = cancerOrder.map(cancer => {
      const item = signatureData.find(d => d.cancer === cancer);
      return item ? item.tmb || 0 : 0;
    });

    traces.push({
      type: 'bar',
      name: signatureAnnotations[signature] || signature,
      x: cancerOrder.map((_, i) => i + 1),
      y: yValues,
      xaxis: 'x',
      yaxis: 'y',
      marker: {
        color: signatureColors[signature] || '#808080'
      },
      showlegend: false,
      hovertemplate: 
        '<b>Cancer:</b> %{customdata.cancer}<br>' +
        '<b>Signature:</b> ' + (signatureAnnotations[signature] || signature) + '<br>' +
        '<b>TMB:</b> %{y:.3f} mutations/Mb<br>' +
        '<extra></extra>',
      customdata: cancerOrder.map(cancer => ({ cancer }))
    });
  });

  // Create dot plot traces (for bottom subplot) 
  dotSignatures.forEach((signature, sigIndex) => {
    const signatureData = dotData.filter(item => item.signature === signature);
    const xValues = [];
    const yValues = [];
    const sizes = [];
    const customData = [];

    cancerOrder.forEach((cancer, cancerIndex) => {
      const item = signatureData.find(d => d.cancer === cancer);
      if (item && item.presence > 0) {
        xValues.push(cancerIndex + 1);
        yValues.push(sigIndex + 1);
        sizes.push(Math.max(item.presence * 40, 6));
        customData.push({
          cancer: cancer,
          signature: signatureAnnotations[signature] || signature,
          presence: item.presence,
          sampleCount: item.sampleCount || 0
        });
      }
    });

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
          color: signatureColors[signature] || '#808080',
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
      range: [0.5, cancerOrder.length + 0.5]
    },
    yaxis: {
      domain: [0.8, 1],
      anchor: 'x',
      title: {
        text: '<b>Mutations per Mb</b>',
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
    margin: {
      l: 350,
      r: 50,
      t: 80,
      b: 200
    },
    plot_bgcolor: 'white',
    paper_bgcolor: 'white',
    annotations: [
      {
        text: 'Circle size represents proportion of tumors with the signature',
        showarrow: false,
        x: 0.5,
        y: -0.15,
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
