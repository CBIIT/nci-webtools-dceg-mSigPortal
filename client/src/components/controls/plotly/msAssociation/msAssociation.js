import { groupBy } from 'lodash';
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
    x: signatureName1data.map((e) => e['exposure']),
    type: 'histogram',
    histnorm: 'density',
    nbinsx: signatureName1data.length,
    xaxis: 'x3',
    yaxis: 'y3',
  };

  const traceSig2 = {
    y: signatureName2data.map((e) => e['exposure']),
    type: 'histogram',
    histnorm: 'density',
    nbinsx: signatureName2data.length,
    xaxis: 'x2',
    yaxis: 'y2',
  };

  const traceMain = {
    x: signatureName1data.map((e) => e['exposure']),
    y: signatureName2data.map((e) => e['exposure']),
    mode: 'markers',
    type: 'scatter',
  };
  console.log(traceSig1);
  const traces = [traceMain, traceSig1, traceSig2];
  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: 350,

    autosize: true,
    title: {
      text: '<b>Mutational Signature Association</b>',
    },

    xaxis2: { domain: [0.8, 1] },
    yaxis2: { anchor: 'x2', domain: [0.0, 0.75] },

    xaxis3: { domain: [0.0, 0.75] },
    yaxis3: { anchor: 'x3', domain: [0.75, 1] },
  };
  return { traces: traces, layout: layout };
}
