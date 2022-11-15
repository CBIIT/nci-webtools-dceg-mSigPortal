import { groupBy } from 'lodash';
import chromosomePositions from './chrPositions.json';

export default function Rainfall(inputData) {
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

  const chrPositions = chromosomePositions[inputData[0].genome];
  const chrOrder = Object.keys(chrPositions);
  const sortChr = (a, b) => chrOrder.indexOf(a.chr) - chrOrder.indexOf(b.chr);

  const groupBySubclass = groupBy(inputData, (e) => e.subclass);
  const data = Object.entries(groupBySubclass)
    .map(([subclass, d]) => ({
      name: subclassInfo[subclass].name,
      data: d.sort(sortChr),
    }))
    .sort((a, b) =>
      subclassPriority.includes(a.name)
        ? -1
        : subclassPriority.includes(b.name)
        ? 1
        : a.name.localeCompare(b.name)
    );

  const traces = data.map((d) => {
    return {
      name: d.name,
      x: d.data.map((e) => e.start + chrPositions[e.chr].start),
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

  const layout = {
    autosize: true,
    height: 500,
    // legend: {
    //   title: {
    //     text: '<b>Classes</b>',
    //   },
    //   x: 1.01,
    //   y: 0.5,
    // },
    xaxis: {
      mirror: true,
      showgrid: false,
      showline: true,
      tickmode: 'array',
      tickvals: Object.values(chrPositions).map((e) => e.center),
      ticktext: chrOrder,
      tickangle: 0,
    },
    yaxis: {
      title: 'Distance between mutations in a single event (log<sub>10</sub>)',
      ticks: 'outside',
      zeroline: false,
      showline: true,
      mirror: true,
      // exponentformat: 'power',
      // tickformat: '',
    },
  };

  const config = {
    responsive: true,
  };
  return { traces, layout, config };
}
