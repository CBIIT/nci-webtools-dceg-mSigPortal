export default function mutationalPatternBar(data) {
  console.log(data);

  const traces = {
    marker: {
      colorscale: [
        [0.0, 'rgb(19,43,67)'],
        [0.1, 'rgb(19,43,67)'],
        [0.1, 'rgb(21,47,72)'],
        [0.2, 'rgb(21,47,72)'],
        [0.2, 'rgb(26,57,86)'],
        [0.3, 'rgb(26,57,86)'],
        [0.3, 'rgb(27,59,88)'],
        [0.4, 'rgb(27,59,88)'],
        [0.4, 'rgb(42,88,128)'],
        [0.5, 'rgb(42,88,128)'],
        [0.5, 'rgb(52,108,155)'],
        [0.6, 'rgb(52,108,155)'],
        [0.6, 'rgb(73,151,213)'],
        [0.7, 'rgb(73,151,213)'],
        [0.7, 'rgb(82,169,237)'],
        [0.8, 'rgb(82,169,237)'],
        [0.8, 'rgb(85,175,245)'],
        [0.9, 'rgb(85,175,245)'],
        [0.9, 'rgb(86,177,247)'],
        [1.0, 'rgb(86,177,247)'],
      ],
      color: data.map((patterndata, index, array) => patterndata.data.length),
    },

    type: 'bar',
    hoverinfo: 'x+y',
    x: data.map((patterndata) => patterndata.pattern),
    y: data.map((patterndata, index, array) => patterndata.data.length),
  };

  console.log(traces);
  var layout = {
    autosize: true,
    height: 500,
    xaxis: {
      showline: true,
      tickangle: -90,
      tickfont: {
        family: 'Arial',
        color: 'black',
      },
      tickmode: 'array',
      categoryorder: 'total descending',
    },
    yaxis: {
      title: {
        text: '<b>Frequency</b>',
        font: {
          family: 'Times New Roman',
        },
      },
    },
  };
  var config = {
    responsive: true,
  };

  return { traces: [traces], layout: layout, config };
}
