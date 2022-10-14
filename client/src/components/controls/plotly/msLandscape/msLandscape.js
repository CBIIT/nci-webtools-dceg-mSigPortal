import { groupBy } from 'lodash';

export default function MsLandscape(data, arg) {
  console.log(data);
  console.log(arg);
  const rawDataActivity = data[0]['data']; //exposure data
  const rawDataSignature = data[1]['data']; //signature data
  const rawDataSpectrum = data[2]['data']; //seqmatrix data
  console.log('activity/exposure data');
  console.log(rawDataActivity);
  console.log('signature data');
  console.log(rawDataSignature);
  console.log('spectrum/seqmatrix data');
  console.log(rawDataSpectrum);

  const groupBySample_activity = groupBy(rawDataActivity, 'sample');
  console.log(groupBySample_activity);
  const groupByMutationType_signature = groupBy(
    rawDataSignature,
    'mutationType'
  );
  console.log(groupByMutationType_signature);

  const groupByMutationType_spectrum = groupBy(rawDataSpectrum, 'mutationType');
  console.log(groupByMutationType_spectrum);

  const traces = [];

  const shapes = [];

  const lines = [];

  let annotations = [];

  const layout = {
    // width: totalCancer > 1 ? null : 350,
    autosize: true,
    height: 500,
    showlegend: false,
    xaxis: {
      showticklabels: false,
      tickfont: {
        size: 10,
      },
      //autorange: false,
      //range: [0, totalCancer],
      // linecolor: 'black',
      linewidth: 2,
      mirror: true,
      tickmode: 'array',
      showgrid: false,
      //tickvals: flatSorted.map((_, i) => i),
      //ticktext: flatSorted.map((_, i) => i),
    },
    yaxis: {
      title: 'Number of Mutations per Megabase<br>(log10)',
      zeroline: false,
      //showline: true,
      linecolor: 'black',
      linewidth: 2,
      mirror: true,
      automargin: true,
      autorange: true,
      //range: [-Math.floor(yMax), Math.floor(yMax)],
    },

    shapes: [...shapes, ...lines],
    annotations: annotations,
  };

  //console.log('layout:');
  //console.log(layout);

  var config = {
    //responsive: true,
  };

  return { traces: [...traces], layout: layout, config };
  //return { traces, layout };
}
