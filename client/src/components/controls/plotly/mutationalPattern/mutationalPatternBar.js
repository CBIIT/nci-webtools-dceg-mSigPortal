import { groupBy } from 'lodash';

export default function mutationalPatternBar(inputData) {
  const groupByPattern = groupBy(inputData, 'pattern');
  const data = Object.entries(groupByPattern).map(([pattern, data]) => {
    return { pattern, data };
  });

  const maxVal = Math.max(...data.map((o) => o.data.length));
  const traces = {
    marker: {
      colorscale: [
        [0, 'rgb(19, 43, 67)'],
        [0.025, 'rgb(19, 43, 67)'],
        [0.05, 'rgb(21, 46, 72)'],
        [0.075, 'rgb(22, 50, 76)'],
        [0.1, 'rgb(24, 53, 81)'],
        [0.125, 'rgb(26, 57, 85)'],
        [0.15, 'rgb(28, 60, 90)'],
        [0.175, 'rgb(29, 64, 95)'],
        [0.2, 'rgb(31, 67, 99)'],
        [0.225, 'rgb(33, 70, 104)'],
        [0.25, 'rgb(35, 74, 109)'],
        [0.275, 'rgb(36, 77, 113)'],
        [0.3, 'rgb(38, 81, 118)'],
        [0.325, 'rgb(40, 84, 122)'],
        [0.35, 'rgb(42, 88, 127)'],
        [0.375, 'rgb(43, 91, 132)'],
        [0.4, 'rgb(45, 95, 136)'],
        [0.425, 'rgb(47, 98, 141)'],
        [0.45, 'rgb(49, 101, 145)'],
        [0.475, 'rgb(50, 105, 150)'],
        [0.5, 'rgb(52, 108, 155)'],
        [0.525, 'rgb(54, 112, 159)'],
        [0.55, 'rgb(56, 115, 164)'],
        [0.575, 'rgb(57, 119, 169)'],
        [0.6, 'rgb(59, 122, 173)'],
        [0.625, 'rgb(61, 125, 178)'],
        [0.65, 'rgb(63, 129, 182)'],
        [0.675, 'rgb(64, 132, 187)'],
        [0.7, 'rgb(66, 136, 192)'],
        [0.725, 'rgb(68, 139, 196)'],
        [0.75, 'rgb(70, 143, 201)'],
        [0.775, 'rgb(71, 146, 205)'],
        [0.8, 'rgb(73, 150, 210)'],
        [0.825, 'rgb(75, 153, 215)'],
        [0.85, 'rgb(77, 156, 219)'],
        [0.875, 'rgb(78, 160, 224)'],
        [0.9, 'rgb(80, 163, 229)'],
        [0.925, 'rgb(82, 167, 233)'],
        [0.95, 'rgb(84, 170, 238)'],
        [0.975, 'rgb(85, 174, 242)'],
        [1.0, 'rgb(87, 177, 247)'],
      ],

      color: data.map((patterndata, index, array) => patterndata.data.length),
    },

    type: 'bar',
    hoverinfo: 'x+y',
    x: data.map((patterndata) => patterndata.pattern),
    y: data.map((patterndata, index, array) => patterndata.data.length),
  };

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
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
      autorange: false,
      range: [0, maxVal + maxVal * 0.2],
      linecolor: 'black',
    },
    title: {
      text: '<b>Frequency of Mutational Pattern</b>',
      font: {
        family: 'Times New Roman',
        size: 18,
      },
    },
  };
  const config = {
    displayModeBar: true,
    responsive: true,
    displaylogo: false,
    toImageButtonOptions: {
      format: 'svg',
      filename: 'Frequency of Mutational Pattern',
    },
  };

  return { traces: [traces], layout, config };
}
