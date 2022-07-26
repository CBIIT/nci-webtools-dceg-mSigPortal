export default function MSPrevalence(groupBySignature, groupBySample) {
  groupBySignature.sort((a, b) => b.samples.length - a.samples.length);
  console.log(groupBySignature);
  console.log(groupBySample);
  const traces = groupBySignature.map((group, groupIndex, array) => ({
    group: group,
    array: array,
    name: group.signatureName,
    type: 'bar',
    x: [group.signatureName],
    y: [group.samples.length / group.totalSamples],
    text: [
      Math.round((group.samples.length / group.totalSamples) * 1000) / 10 + '%',
    ],
    textposition: 'outside',
    hoverinfo: 'x+y',
    showlegend: false,
  }));
  console.log(traces);

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    bargap: 0.3,
    height: 450,
    // width: 1080,
    autosize: true,

    xaxis: {
      showline: true,
      tickangle: -90,
      tickfont: {
        family: 'Arial, monospace',
      },
      tickmode: 'array',

      linecolor: '#E0E0E0',
      linewidth: 1,
      mirror: 'all',
    },
    yaxis: {
      title: {
        text: '<b>Frequency (%)</b>',
        font: {
          family: 'Times New Roman',
        },
      },

      range: [0, 1.1],
      ticks: 'inside',
      tickcolor: '#D3D3D3',
      linecolor: '#D3D3D3',
      linewidth: 1,
      mirror: 'all',
      tickformat: ',.0%',
      showgrid: true,
      gridcolor: '#F5F5F5',
    },
  };

  return { traces: [...traces], layout: layout };
}
