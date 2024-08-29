import { groupBy } from 'lodash';

// genome info is an object containing the start and end positions for each chromosome within the genome
export default function Rainfall(inputData, genomeInfo) {
  const subclassInfo = {
    ClassIA: { name: 'DBS', color: 'red' },
    ClassIB: { name: 'MBS', color: 'black' },
    ClassIC: { name: 'OMIKLI', color: 'green' },
    ClassII: { name: 'KATAEGIS', color: 'orange' },
    ClassIII: { name: 'OTHER', color: 'blue' },
    'Non-clust': { name: 'Non-clust', color: 'grey', text: 'white' },
    Simulation: { name: 'Simulation', color: 'grey', text: 'white' },
    Clust: { name: 'Clust', color: 'orange' },
  };
  const chrOrder = Object.keys(genomeInfo);

  // group data by classes and sort chromosomes
  const groupBySubclass = groupBy(inputData, (e) => e.subclass);
  const data = Object.entries(groupBySubclass).map(([subclass, d]) => ({
    name: subclassInfo[subclass].name,
    data: d.sort((a, b) => chrOrder.indexOf(a.chr) - chrOrder.indexOf(b.chr)),
  }));

  // extra margin for visualizing data and window size for bin
  const binSize = 10000000;

  // count the number of mutations along the length of the genome per bin size
  const sortedData = data.reduce((a, e) => [...a, ...e.data], []);
  const densityPlot = {
    name: 'Density',
    x: sortedData.map((e) => e.start + genomeInfo[e.chr].start),
    y: sortedData.map((e) => e['IMD']),
    histfunc: 'count',
    type: 'histogram',
    xbins: { size: binSize },
    color: 'blue',
    showlegend: false,
    xaxis: 'x',
    yaxis: 'y2',
  };

  const rainfallTraces = data.map((d) => {
    return {
      name: d.name,
      x: d.data.map((e) => e.start + genomeInfo[e.chr].start),
      y: d.data.map((e) => e['IMD']),
      text: d.data.map((e) =>
        [
          `Chromosome: ${e.chr}`,
          `IMD: ${e['IMD']}`,
          `Sample: ${e.sample || 'N/A'}`,
          `Subclass: ${e.subclass} (${d.name})`,
        ].join('<br>')
      ),
      hovertemplate: '%{text}<extra></extra>',
      marker: { color: subclassInfo[d.data[0].subclass].color },
      hoverlabel: {
        font: { color: subclassInfo[d.data[0].subclass].text },
      },
      mode: 'markers',
      type: 'scatter',
      showlegend: true,
    };
  });

  const layout = {
    autosize: true,
    height: 800,
    xaxis: {
      title: { text: 'Chromosome', font: { size: 18 }, standoff: 20 },
      mirror: true,
      showline: true,
      zeroline: false,
      showgrid: false,
      range: [-binSize, genomeInfo['Y'].end + binSize],
      tickmode: 'array',
      tickvals: Object.values(genomeInfo).map((e) => e.end - e.len / 2),
      ticktext: chrOrder,
      tickangle: 0,
    },
    yaxis: {
      title: {
        text: 'Distance between mutations in a single event (log<sub>10</sub>)',
        font: { size: 18 },
      },
      domain: [0, 0.78],
      mirror: true,
      showline: true,
      zeroline: false,
      ticks: 'outside',
      type: 'log',
      exponentformat: 'power',
    },
    yaxis2: {
      title: { text: 'Density', font: { size: 18 } },
      domain: [0.8, 1],
    },
    shapes: [
      // chromosome dividers
      ...Object.values(genomeInfo)
        .slice(0, -1)
        .map((e) => ({
          type: 'line',
          xref: 'x',
          yref: 'y',
          x0: -binSize,
          x1: genomeInfo['Y'].end + binSize,
          y0: Math.pow(10, 3),
          y1: Math.pow(10, 3),
          line: { width: 0.5, color: 'red' },
        })),
    ],
    legend: {
      orientation: 'h',
      xanchor: 'center',
      x: 0.5,
    },
  };

  const config = {
    responsive: true,
  };

  // data table - filter out subclass: non-clust
  const table = inputData
    .filter((e) => e.subclass != 'Non-clust')
    .map((e) => ({
      ...e,
      subclass: `${e.subclass} (${subclassInfo[e.subclass].name})`,
    }));
  return { traces: [densityPlot, ...rainfallTraces], layout, config, table };
}
