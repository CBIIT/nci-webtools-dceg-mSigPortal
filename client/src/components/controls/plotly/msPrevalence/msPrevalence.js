export default function MSPrevalence(
  groupBySignature,
  groupBySample,
  mutation
) {
  groupBySignature.sort((a, b) => b.samples.length - a.samples.length);
  console.log(groupBySignature);
  console.log(groupBySample);
  console.log(mutation);

  let minumumNumber;
  mutation === 'null'
    ? (minumumNumber = 100)
    : (minumumNumber = parseInt(mutation));

  const tracesPie = {
    type: 'pie',
    labels: groupBySignature.map((group) => group.signatureName),
    values: groupBySignature.map((group) => group.samples.length),
    textposition: 'inside',
    textinfo: 'label+percent',
    showlegend: false,
  };
  console.log(tracesPie);

  const tracesBar = groupBySignature.map((group, groupIndex, array) => ({
    group: group,
    array: array,
    name: group.signatureName,
    type: 'bar',
    minumumNumber: minumumNumber,
    lenght: group.samples.filter((e) => e.exposure >= minumumNumber),

    x: [group.signatureName],
    y: [
      group.samples.filter((e) => e.exposure >= minumumNumber).length /
        group.totalSamples,
    ],
    text: [
      Math.round(
        (group.samples.filter((e) => e.exposure >= minumumNumber).length /
          group.totalSamples) *
          1000
      ) /
        10 +
        '%',
    ],
    textposition: 'outside',
    xaxis: 'x2',
    yaxis: 'y2',
    hoverinfo: 'x2+y2',
    showlegend: false,
    domain: {
      row: 0,
      column: 1,
    },
  }));
  console.log(tracesBar);

  const traces = [tracesPie, ...tracesBar];

  console.log(traces);

  const titleAnnotations = [
    {
      xref: 'paper',
      yref: 'paper',
      showarrow: false,
      x: 0.15,
      y: 1.1,
      xanchor: 'top',
      text: '<b>Prevalence by mutations</b>',
      font: {
        size: 18,
        family: 'Arial',
      },
    },
    {
      xref: 'x2 domain',
      yref: 'paper',
      showarrow: false,
      x: 0.5,
      y: 1.1,
      xanchor: 'top',
      text: '<b>Prevalence by samples</b>',
      font: {
        size: 18,
        family: 'Arial',
      },
    },
  ];

  const layout = {
    grid: { rows: 1, columns: 2 },
    hoverlabel: { bgcolor: '#FFF' },
    height: 450,
    //width: 500,
    autosize: true,
    title: '<b>Prevalence of Mutational Signatures</b>',
    xaxis2: {
      showline: true,
      tickangle: -90,
      tickfont: {
        family: 'Arial, monospace',
      },
      tickmode: 'array',

      linecolor: 'black',
      linewidth: 1,

      categoryorder: 'total descending',
    },
    yaxis2: {
      title: {
        text: '<b>Frequency (%)</b>',
        font: {
          family: 'Times New Roman',
        },
      },

      range: [0, 1.1],
      ticks: 'inside',
      tickcolor: '#D3D3D3',
      linecolor: 'black',
      linewidth: 1,

      tickformat: ',.0%',
      showgrid: true,
      gridcolor: '#F5F5F5',
    },
    annotations: titleAnnotations,
  };

  return { traces: traces, layout: layout };
}
