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

  const traceSig1 = {
    x: signatureName1data.map((e) => Math.log(e['exposure'] + 1)),
    name: signatureName1,
    type: 'histogram',
    histnorm: 'density',
    nbinsx: signatureName1data.length,
    yaxis: 'y2',
    marker: { color: '#019E72' },
  };
  // const traceSig1 = {
  //   type: 'violin',
  //   spanmode: 'soft',
  //   hoveron: 'points+kde',
  //   // span: [0, 5],
  //   side: 'positive', //positive side means right for vertical violin plots
  //   //y: 4,
  //   x: signatureName1data.map((e) => e['exposure']),
  //   points: 'none',
  //   box: {
  //     visible: false,
  //   },
  //   boxpoints: false,
  //   line: {
  //     color: 'black',
  //   },
  //   fillcolor: '#019E72',
  //   scalemode: 'count',
  //   marker: {
  //     line: {
  //       width: 0,
  //     },
  //     symbol: 'line-ns',
  //   },
  //   opacity: 0.6,
  //   meanline: {
  //     visible: true,
  //   },
  //   y0: signatureName1,
  //   xaxis: 'x',
  //   yaxis: 'y3',
  // };

  const traceSig2 = {
    y: signatureName2data.map((e) => Math.log(e['exposure'] + 1)),
    name: signatureName2,
    type: 'histogram',
    histnorm: 'density',
    nbinsy: signatureName2data.length,
    xaxis: 'x2',
    marker: { color: '#D55E00' },
  };

  // const traceSig2 = {
  //   type: 'violin',
  //   spanmode: 'soft',
  //   hoveron: 'points+kde',
  //   // span: [0, 5],
  //   side: 'positive', //positive side means right for vertical violin plots
  //   //y: 4,
  //   y: signatureName2data.map((e) => e['exposure']),
  //   points: 'none',
  //   box: {
  //     visible: false,
  //   },
  //   boxpoints: false,
  //   line: {
  //     color: 'black',
  //   },
  //   fillcolor: '#D55E00',
  //   scalemode: 'count',
  //   marker: {
  //     line: {
  //       width: 0,
  //     },
  //     symbol: 'line-ns',
  //   },
  //   opacity: 0.6,
  //   meanline: {
  //     visible: true,
  //   },
  //   x0: signatureName2,
  //   xaxis: 'x2',
  //   yaxis: 'y',
  // };

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
  console.log(traceSig1);
  const traces = [traceMain, traceSig1, traceSig2];
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
      domain: [0.0, 0.75],
      showgrid: true,
      title: {
        text: '<b>Number of mutations in ' + signatureName1 + ' (log10)</b>',
      },
    },
    yaxis: {
      anchor: 'x',
      domain: [0.0, 0.75],
      title: {
        text: '<b>Number of mutations in ' + signatureName2 + ' (log10)</b>',
      },
      showgrid: true,
    },

    xaxis2: { domain: [0.85, 1] },
    yaxis2: { anchor: 'x', domain: [0.85, 1] },
  };
  return { traces: traces, layout: layout };
}
