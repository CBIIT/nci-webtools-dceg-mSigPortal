export default function MsLandscape(data, tmbTabName, signatureName) {
  console.log(data);
  console.log(signatureName);
  function average(arr) {
    const sum = arr.reduce((a, b) => a + b, 0);
    return sum / arr.length || 0;
  }

  const totalCancer = data.length;

  const absYValue = data
    .map((o) => o.samples.map((e) => Math.abs(e.burden)))
    .flat();
  const yMax = Math.max(...absYValue);

  const traces = [];

  //console.log('traces:--');
  //console.log(traces);

  const topLabel = [];
  console.log('top label:--');
  console.log(topLabel);

  const bottoLabel1 = [];

  const bottoLabelline = [];

  const bottoLabel2 = [];

  const shapes = [];
  //console.log('shapes:--');
  //console.log(shapes);

  const lines = [];

  const signatureNameAnnotation = [];
  let annotations = [];
  signatureName != null
    ? (annotations = [
        ...topLabel,
        ...bottoLabel1,
        ...bottoLabel2,
        signatureNameAnnotation,
      ])
    : (annotations = [...topLabel, ...bottoLabel1, ...bottoLabel2]);

  const layout = {
    width: totalCancer > 1 ? null : 350,
    autosize: true,
    height: 500,
    showlegend: false,
    xaxis: {
      showticklabels: false,
      tickfont: {
        size: 10,
      },
      autorange: false,
      range: [0, totalCancer],
      linecolor: 'black',
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
      range: [-Math.floor(yMax), Math.floor(yMax)],
    },

    shapes: [...shapes, ...lines, ...bottoLabelline],
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
