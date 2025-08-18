import { groupBy } from 'lodash';

export default function SATSTMBBarChart(data, options = {}) {
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

  // Group by cancer type and sort by total TMB
  const cancerGroups = groupBy(data, 'cancer');
  const cancerOrder = Object.entries(cancerGroups)
    .map(([cancer, items]) => ({
      cancer,
      totalTMB: items.reduce((sum, item) => sum + (item.tmb || 0), 0)
    }))
    .sort((a, b) => b.totalTMB - a.totalTMB)
    .map(item => item.cancer);

  // Get all unique signatures
  const signatures = [...new Set(data.map(item => item.signature))];
  
  // Create traces for stacked bar chart
  const traces = signatures.map(signature => {
    const signatureData = data.filter(item => item.signature === signature);
    const yValues = cancerOrder.map(cancer => {
      const item = signatureData.find(d => d.cancer === cancer);
      return item ? item.tmb || 0 : 0;
    });

    return {
      type: 'bar',
      name: signatureAnnotations[signature] || signature,
      x: cancerOrder,
      y: yValues,
      marker: {
        color: signatureColors[signature] || '#808080'
      },
      hovertemplate: 
        '<b>Cancer:</b> %{x}<br>' +
        '<b>Signature:</b> ' + (signatureAnnotations[signature] || signature) + '<br>' +
        '<b>TMB:</b> %{y:.3f} mutations/Mb<br>' +
        '<extra></extra>',
      showlegend: true
    };
  });

  const layout = {
    title: {
      text: '<b>Mutations per Megabase by Cancer Type</b>',
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
      tickangle: 45,
      tickfont: {
        size: 10
      },
      categoryorder: 'array',
      categoryarray: cancerOrder
    },
    yaxis: {
      title: {
        text: '<b>Mutations per Mb</b>',
        font: {
          family: 'Times New Roman',
          size: 14
        }
      },
      showgrid: true,
      gridcolor: 'rgba(128,128,128,0.2)'
    },
    barmode: 'stack',
    autosize: true,
    height: 600,
    margin: {
      l: 80,
      r: 50,
      t: 80,
      b: 150
    },
    plot_bgcolor: 'white',
    paper_bgcolor: 'white',
    legend: {
      orientation: 'v',
      x: 1.02,
      y: 1,
      font: {
        size: 10
      }
    }
  };

  const config = {
    displayModeBar: true,
    responsive: true,
    displaylogo: false,
    toImageButtonOptions: {
      format: 'svg',
      filename: 'SATS_TMB_BarChart',
      width: 1200,
      height: 600
    },
    modeBarButtonsToRemove: ['pan2d', 'lasso2d', 'select2d']
  };

  return { 
    traces, 
    layout, 
    config 
  };
}
