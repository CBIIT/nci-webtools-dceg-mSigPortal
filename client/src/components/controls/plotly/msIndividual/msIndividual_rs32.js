export default function RS32(rawData, arg) {
  console.log(rawData);
  console.log(arg);
  const annotation = {
    xref: 'paper',
    yref: 'paper',

    x: 0.5,
    y: 0.5,
    text:
      'Signature SetName: <b>' +
      arg.params_activity.signatureSetName +
      '</b> is not supported in MS Individual',
    font: {
      size: 15,
      color: 'red',
    },
    showarrow: false,
    align: 'center',
  };
  const traces = {};
  const layout = {
    height: 100,
    autosize: true,
    xaxis: {
      showticklabels: false,
      showline: false,
      zeroline: false,
    },
    yaxis: { showticklabels: false, zeroline: false },
    annotations: [annotation],
  };
  return { traces, layout };
}
