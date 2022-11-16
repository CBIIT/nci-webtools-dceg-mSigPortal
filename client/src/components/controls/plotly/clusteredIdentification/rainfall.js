import { groupBy } from 'lodash';

// genome info is an object containing the start and end positions for each chromosome within the genome
export default function Rainfall(inputData, genomeInfo) {
  function getRange(array) {
    let min = array[0];
    let max = array[0];
    for (let item of array) {
      if (item < min) min = item;
      if (item > max) max = item;
    }
    return [min, max];
  }

  const subclassInfo = {
    ClassIA: { name: 'DBS', color: 'red' },
    ClassIB: { name: 'MBS', color: 'black' },
    ClassIC: { name: 'OMIKLI', color: 'green' },
    ClassII: { name: 'KATAEGIS', color: 'orange' },
    ClassIII: { name: 'OTHER', color: 'blue' },
    'Non-clust': { name: 'Non-clust', color: 'grey' },
    Simulation: { name: 'Simulation', color: 'grey' },
    Clust: { name: 'Clust', color: 'orange' },
  };
  const subclassPriority = ['Non-clust', 'Simulation'];
  const chrOrder = Object.keys(genomeInfo);

  // group data by classes and sort chromosomes
  const groupBySubclass = groupBy(inputData, (e) => e.subclass);
  const data = Object.entries(groupBySubclass)
    .map(([subclass, d]) => ({
      name: subclassInfo[subclass].name,
      data: d.sort((a, b) => chrOrder.indexOf(a.chr) - chrOrder.indexOf(b.chr)),
    }))
    .sort((a, b) =>
      subclassPriority.includes(a.name)
        ? -1
        : subclassPriority.includes(b.name)
        ? 1
        : a.name.localeCompare(b.name)
    );

  // determine y coordinates for each bin
  const yCoordinates = inputData.map((e) => e['IMD']);
  const [_, yMax] = getRange(yCoordinates);

  const traces = data.map((d) => {
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
      marker: { color: d.data.map((e) => subclassInfo[e.subclass].color) },
      mode: 'markers',
      type: 'scatter',
      showlegend: true,
    };
  });

  const bufferMargin = 1000000;
  const layout = {
    autosize: true,
    height: 500,
    xaxis: {
      showgrid: false,
      showline: true,
      tickmode: 'array',
      tickvals: Object.values(genomeInfo).map((e) => e.end - e.len / 2),
      ticktext: chrOrder,
      tickangle: 0,
    },
    yaxis: {
      title: 'Distance between mutations in a single event (log<sub>10</sub>)',
      zeroline: false,
      ticks: 'outside',
      // range: [-bufferMargin, yMax + bufferMargin],
      type: 'log',
      // exponentformat: 'power',
      // tickformat: '',
    },
    shapes: [
      // chromosome dividers
      // ...Object.values(genomeInfo).map((e) => ({
      //   type: 'line',
      //   x0: e.end,
      //   x1: e.end,
      //   y0: 0,
      //   y1: yMax,
      //   line: { width: 1 },
      // })),
    ],
  };

  const config = {
    responsive: true,
  };
  return { traces, layout, config };
}
