import { groupByCustom } from '../../utils/utils';
export default function MsAssociation(data, arg) {
  console.log(data);
  console.log(arg);
  const [signatureName1, signatureName2] = arg.signatureName.split(';');
  console.log(signatureName1);
  console.log(signatureName2);
  const groupBySignatureName = groupByCustom(data, (e) => e.signatureName);
  const signatureName1data = groupBySignatureName.get(signatureName1);
  const signatureName2data = groupBySignatureName.get(signatureName2);

  console.log(groupBySignatureName);
  console.log(signatureName1data);
  console.log(signatureName2data);
  function median(numbers) {
    const sorted = Array.from(numbers).sort((a, b) => a - b);
    const middle = Math.floor(sorted.length / 2);

    if (sorted.length % 2 === 0) {
      return (sorted[middle - 1] + sorted[middle]) / 2;
    }

    return sorted[middle];
  }

  function getAvg(grades) {
    const total = grades.reduce((acc, c) => acc + c, 0);
    return total / grades.length;
  }
  const minX = Math.min(
    ...signatureName1data.map((e) => Math.log(e['exposure'] + 1))
  );
  console.log(minX);

  const maxX = Math.max(
    ...signatureName1data.map((e) => Math.log(e['exposure'] + 1))
  );
  console.log(maxX);

  const minY = Math.min(
    ...signatureName2data.map((e) => Math.log(e['exposure'] + 1))
  );
  console.log(minY);

  const maxY = Math.max(
    ...signatureName2data.map((e) => Math.log(e['exposure'] + 1))
  );
  console.log(maxY);

  const medianY = median([minY, maxY]);
  console.log(medianY);

  const medianY2 = median([
    ...signatureName2data.map((e) => Math.log(e['exposure'] + 1)),
  ]);
  console.log(medianY2);

  const avgY0 = getAvg([
    ...signatureName2data.map((e) => Math.log(e['exposure'] + 1)),
  ]);
  console.log(avgY0);

  const avgY1 = getAvg([minY, maxY]);
  console.log(avgY1);
  const traceSig1 = {
    x: signatureName1data.map((e) => Math.log(e['exposure'] + 1)),
    name: signatureName1,
    type: 'histogram',
    histnorm: 'density',
    nbinsx: signatureName1data.length - 1,
    yaxis: 'y2',
    marker: { color: '#019E72', line: { color: 'black', width: 1 } },
  };

  const traceSig2 = {
    y: signatureName2data.map((e) => Math.log(e['exposure'] + 1)),
    name: signatureName2,
    type: 'histogram',
    histnorm: 'density',
    nbinsy: signatureName2data.length - 1,
    xaxis: 'x2',
    marker: { color: '#D55E00', line: { color: 'black', width: 1 } },
  };

  const traceMain = {
    x: signatureName1data.map((e) => Math.log(e['exposure'] + 1)),
    y: signatureName2data.map((e) => Math.log(e['exposure'] + 1)),
    mode: 'markers',
    type: 'scatter',
    marker: {
      color: '#A3A3A3',
      size: 10,
    },
    opacity: 0.9,
  };

  const traceLine = {
    x: [minX, maxX],
    y: [avgY1, avgY0],
    mode: 'lines',
  };

  const traces = [traceMain, traceLine, traceSig1, traceSig2];
  const layout = {
    showlegend: false,
    hoverlabel: { bgcolor: '#FFF' },
    height: 900,
    bargap: 0,
    autosize: true,
    title: {
      text: '<b>Mutational Signature Association</b>',
    },
    xaxis: {
      domain: [0.0, 0.83],
      showgrid: true,
      title: {
        text: '<b>Number of mutations in ' + signatureName1 + ' (log10)</b>',
      },
    },
    yaxis: {
      anchor: 'x',
      domain: [0.0, 0.83],
      title: {
        text: '<b>Number of mutations in ' + signatureName2 + ' (log10)</b>',
      },
      showgrid: true,
    },

    xaxis2: { domain: [0.85, 1], zerolinecolor: '#EBEBEB' },
    yaxis2: { anchor: 'x', domain: [0.85, 1], zerolinecolor: '#EBEBEB' },
  };
  return { traces: traces, layout: layout };
}
