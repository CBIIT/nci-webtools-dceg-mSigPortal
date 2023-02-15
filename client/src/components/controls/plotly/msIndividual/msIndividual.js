import { groupBy } from 'lodash';
export default function MsIndividual(data, arg) {
  console.log('MS Individual Plot');
  console.log(data);
  console.log(arg);
  const exposureData = data[0].data;
  console.log(exposureData);
  const signatureData = data[1].data;
  console.log(signatureData);
  const segmatrixData = data[2].data;
  console.log(segmatrixData);

  const exposure_groupBySignature = groupBy(
    exposureData.filter((o) => o['exposure'] > 0),
    'signatureName'
  );

  const signature_groupBySignature = groupBy(signatureData, 'signatureName');
  console.log(signature_groupBySignature);
  console.log(exposure_groupBySignature);

  const seqmatrix_groupByMutationType = groupBy(segmatrixData, 'mutationType');
  console.log(seqmatrix_groupByMutationType);
  var trace1 = {};

  var trace2 = {};

  var traces = [trace1, trace2];
  const layout = {
    showlegend: true,
    hoverlabel: { bgcolor: '#FFF' },
    height: 900,
  };
  return { traces: traces, layout: layout };
}
