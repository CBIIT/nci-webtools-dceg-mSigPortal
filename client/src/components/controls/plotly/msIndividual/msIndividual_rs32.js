export default function RS32(rawData, arg) {
  console.log(rawData);
  console.log(arg);
  const annotation = {};
  const traces = {};
  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: 500,
    autosize: true,
    xaxis: {
      showticklabels: true,
      showline: true,
      tickangle: -90,
      tickfont: { size: 11 },
      tickmode: 'array',

      linecolor: 'black',
      linewidth: 1,
      mirror: 'all',
    },
    yaxis: {
      title: {
        text: '<b>Percentage(%)</b>',
        font: {
          family: 'Times New Roman',
          size: 18,
        },
      },
      autorange: false,
      tickformat: ',.1%',
      linecolor: 'black',
      linewidth: 1,
      mirror: true,
    },

    shapes: [],
    annotations: [],
  };
  return { traces, layout };
}
