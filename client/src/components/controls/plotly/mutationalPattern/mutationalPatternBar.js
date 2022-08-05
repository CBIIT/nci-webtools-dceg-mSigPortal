export default function mutationalPatternBar(data) {
  console.log(data);

  const traces = data.map((patterndata, index, array) => ({
    data: patterndata,
    name: patterndata.pattern,
    colorscale: 'Greens',
    type: 'bar',
    hoverinfo: 'x+y',
    x: [patterndata.pattern],
    y: [patterndata.data.length],
    showlegend: false,
  }));

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

  return { traces: traces, layout: layout, config };
}
