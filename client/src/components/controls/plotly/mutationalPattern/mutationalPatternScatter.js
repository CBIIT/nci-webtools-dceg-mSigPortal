import { groupBy } from 'lodash';
import { mutationalPatternColors } from '../../utils/colors';

export default function mutationalPatternScatter(inputData, arg) {
  const { pattern } = arg;
  const type =
    pattern.substring(1, 2) + pattern.substring(3, 4) + pattern.substring(5, 6);
  const subtype1 = pattern.split('>')[0][0];
  const subtype2 = pattern.split('>')[0].slice(-1);
  // const subtype1 = pattern.substring(0, 1);
  // const subtype2 = pattern.substring(2, 3);
  const pattern1 = type + ' context';
  const pattern2 = pattern + ' other context';
  const tmpdata0 = Object.values(
    groupBy(inputData, (e) => `${e.study}_${e.sample}`)
  ).map((samples) => {
    return {
      study: `${samples[0].study}`,
      sample: `${samples[0].sample}`,
      type: type,
      cancer: samples[0].cancer,
      total: samples.reduce((acc, e) => acc + e.mutations, 0),
    };
  });

  const mutationTypeFilter = inputData.filter(
    (e) => e.mutationType.substring(2, 5) === type
  );

  const groupByStudySampleType = groupBy(
    mutationTypeFilter,
    (e) => `${e.study}_${e.sample}`
  );

  const tmpdata1 = Object.values(groupByStudySampleType).map((samples) => {
    return {
      study: `${samples[0].study}`,
      sample: `${samples[0].sample}`,
      n0: samples.reduce((acc, e) => acc + e.mutations, 0),
    };
  });
  function iupac(base) {
    let result = [];
    if (base === 'T' || base === 'U') {
      result = ['T'];
    } else if (base === 'M') {
      result = ['A', 'C'];
    } else if (base === 'R') {
      result = ['A', 'G'];
    } else if (base === 'W') {
      result = ['A', 'T'];
    } else if (base === 'S') {
      result = ['C', 'G'];
    } else if (base === 'Y') {
      result = ['C', 'T'];
    } else if (base === 'K') {
      result = ['G', 'T'];
    } else if (base === 'V') {
      result = ['A', 'C', 'G'];
    } else if (base === 'H') {
      result = ['A', 'C', 'T'];
    } else if (base === 'D') {
      result = ['A', 'G', 'T'];
    } else if (base === 'B') {
      result = ['C', 'G', 'T'];
    } else if (base === 'N') {
      result = ['G', 'A', 'T', 'C'];
    } else {
      result = [base];
    }
    return result;
  }
  const mutationTypeSubTypesFilter = inputData.filter(
    (e) =>
      e.mutationType.substring(2, 5) === type &&
      iupac(subtype1).includes(e.mutationType.substring(0, 1)) &&
      iupac(subtype2).includes(e.mutationType.substring(6, 7))
  );

  const groupByStudySampleTypeFilter = groupBy(
    mutationTypeSubTypesFilter,
    (e) => `${e.study}_${e.sample}`
  );

  const tmpdata2 = Object.values(groupByStudySampleTypeFilter).map(
    (samples) => {
      return {
        study: `${samples[0].study}`,
        sample: `${samples[0].sample}`,
        cancer: `${samples[0].cancer}`,
        type: type,
        n1: samples.reduce((acc, e) => acc + e.mutations, 0),
      };
    }
  );
  const merge = (a1, a2) =>
    a1.map((itm) => ({
      ...a2.find(
        (item) => item.sample === itm.sample && item.study === itm.study && item
      ),
      ...itm,
    }));

  const result0 = merge(tmpdata1, tmpdata0);
  const result1 = merge(tmpdata2, result0);

  const result = result1
    .map((e) => {
      let n2 = e.n0 - e.n1;

      return {
        study: e.study,
        sample: e.sample,
        cancer: e.cancer,
        total: e.total,
        type: pattern,
        n0: e.n0,
        n1: e.n1 / e.total,
        n2: n2 / e.total,
      };
    })
    .filter((e) => e.total > 200);

  const groupByCancer = groupBy(result, (e) => `${e.cancer}`);

  const maxTotal = Math.max(...result.map((o) => o.total));

  function format_output(output) {
    let n = Math.log(output) / Math.LN10;
    let x = Math.round(n);
    if (x < 0) x = 0;
    let out = '1';
    for (let i = 0; i <= x; i++) out += '0';
    return out / 10;
  }

  const maxMutationFilter = format_output(maxTotal);

  const data1 = result.filter((o) => o.total < maxMutationFilter / 100);
  const data2 = result.filter(
    (o) => o.total < maxMutationFilter / 10 && o.total > maxMutationFilter / 100
  );
  const data3 = result.filter((o) => o.total > maxMutationFilter / 10);

  const additionalData = [];

  if (maxMutationFilter >= 100000) {
    additionalData.push(300);
    additionalData.push(1000);
    additionalData.push(3000);
    additionalData.push(10000);
    additionalData.push(30000);
    additionalData.push(100000);
  } else if (maxMutationFilter >= 10000) {
    additionalData.push(300);
    additionalData.push(1000);
    additionalData.push(3000);
    additionalData.push(10000);
  } else if (maxMutationFilter >= 3000) {
    additionalData.push(300);
    additionalData.push(1000);
    additionalData.push(3000);
  } else {
    additionalData.push(300);
    additionalData.push(1000);
  }

  const tracesSize = [];

  let trace1 = {
    name: (maxMutationFilter / 100).toLocaleString(undefined),
    x: data1.map((e) => e.n1),
    y: data1.map((e) => e.n2),
    mode: 'markers',
    type: 'scatter',
    opacity: 1,
    marker: {
      color: 'white',
      line: {
        color: 'black',
        width: 1,
      },
      size: 5,
    },
    showlegend: true,
    legendrank: 2,
    legendgroup: 'b',
    legendgrouptitle: {
      text: 'Number of mutations',
    },
  };
  let trace2 = {
    name: (maxMutationFilter / 10).toLocaleString(undefined),
    x: data2.map((e) => e.n1),
    y: data2.map((e) => e.n2),
    mode: 'markers',
    type: 'scatter',
    opacity: 1,
    marker: {
      color: 'white',
      line: {
        color: 'black',
        width: 1,
      },
      size: 10,
    },
    showlegend: true,
    legendrank: 2,
    legendgroup: 'b',
    legendgrouptitle: {
      text: 'Number of mutations',
    },
  };
  let trace3 = {
    name: maxMutationFilter.toLocaleString(undefined),
    x: data3.map((e) => e.n1),
    y: data3.map((e) => e.n2),
    mode: 'markers',
    type: 'scatter',
    opacity: 1,
    marker: {
      color: 'white',
      line: {
        color: 'black',
        width: 1,
      },
      size: 15,
    },
    showlegend: true,
    legendrank: 2,
    legendgroup: 'b',
    legendgrouptitle: {
      text: 'Number of mutations',
    },
  };

  additionalData.forEach((additionalFilter, index) => {
    const additionalDataFiltered = result.filter(
      (o) => o.total > additionalFilter
    );
    const markerSizes = [];

    if (additionalData.length === 6) {
      markerSizes.push(5, 7, 9, 11, 13, 15);
    } else if (additionalData.length === 4) {
      markerSizes.push(5, 8, 12, 15);
    } else if (additionalData.length === 3) {
      markerSizes.push(5, 10, 15);
    } else {
      markerSizes.push(10);
    }

    const additionalTrace = {
      name: additionalFilter.toLocaleString(undefined),
      x: additionalDataFiltered.map((e) => e.n1),
      y: additionalDataFiltered.map((e) => e.n2),
      customdata: additionalDataFiltered.map((e) => ({
        sample: e.sample,
        total: e.total,
      })),
      mode: 'markers',
      type: 'scatter',
      opacity: 1,
      marker: {
        color: 'white',
        line: {
          color: 'black',
          width: 1,
        },
        size: markerSizes[index],
      },

      hoverinfo: 'skip',
      showlegend: true,
      legendgroup: 'size',
      legendgrouptitle: {
        text: 'Number of mutations',
      },
    };
    tracesSize.push(additionalTrace);
  });

  // Create an array to store all the traces
  let scatterTraces = [];
  let scatterLegdend = [];
  // Create an array of keys from the originalColorPallet object

  // Loop through each group in the groupByCancer result
  for (let group in groupByCancer) {
    let result = groupByCancer[group];

    let trace = {
      name:
        result[0]?.study !== undefined
          ? result[0].study + '@' + result[0].cancer
          : 'Input',
      x: result.map((e) => e.n1),
      y: result.map((e) => e.n2),
      customdata: result.map((e) => ({
        sample: e.sample,
        total: e.total,
        study:
          e?.study !== undefined
            ? result[0].study + '@' + result[0].cancer
            : 'Input',
      })),
      mode: 'markers',
      type: 'scatter',
      opacity: 1,
      marker: {
        color:
          scatterTraces.length < mutationalPatternColors.length
            ? mutationalPatternColors[scatterTraces.length]
            : mutationalPatternColors[
                scatterTraces.length % mutationalPatternColors.length
              ], // Use colors from mutationalPatternColors array in a circular manner
        line: {
          color: 'black',
          width: 1,
        },
        // size: result.map((e, i, a) =>
        //   e.total < maxMutationFilter / 100
        //     ? 5
        //     : e.total < maxMutationFilter / 10
        //     ? 10
        //     : 15
        // ),
        size: result.map((e, i, a) =>
          getMarkerSize(additionalData.length, e.total)
        ),
      },
      hovertemplate:
        '<b>Sample: ' +
        '</b>' +
        '%{customdata.sample} <br>' +
        '<b>Study: </b>' +
        '%{customdata.study} <br>' +
        '<b>' +
        pattern2 +
        ':</b>' +
        ' %{x} <br>' +
        '<b>' +
        pattern1 +
        ':</b>' +
        ' %{y}<br>' +
        '<b>Total:</b> %{customdata.total}<extra></extra>',
      showlegend: false,
      legendrank: 1,
      legendgroup: 'a',
      legendgrouptitle: {
        text: 'Study',
      },
    };

    let traceLedgend = {
      name:
        result[0]?.study !== undefined
          ? result[0].study + '@' + result[0].cancer
          : 'Input',
      x: [null],
      y: [null],
      mode: 'markers',
      type: 'scatter',
      opacity: 1,
      marker: {
        color: scatterTraces.length === 0 ? 'green' : getRandomColor(), // Generate a random color for each trace
        line: {
          color: 'black',
          width: 1,
        },
        size: 12,
      },

      showlegend: true,
      legendrank: 1,
      legendgroup: 'a',
      legendgrouptitle: {
        text: 'Study',
      },
    };

    scatterTraces.push(trace);
    scatterLegdend.push(traceLedgend);
  }
  // Function to generate a random color
  function getRandomColor() {
    return '#' + Math.floor(Math.random() * 16777215).toString(16);
  }

  // Add additional colors to mutationalPatternColors if needed
  while (mutationalPatternColors.length < scatterTraces.length) {
    mutationalPatternColors.push(getRandomColor());
  }

  function getMarkerSize(additionalDataLength, total) {
    if (additionalDataLength === 6) {
      if (total < 300) {
        return 5;
      } else if (total < 1000) {
        return 7;
      } else if (total < 3000) {
        return 9;
      } else if (total < 10000) {
        return 11;
      } else if (total < 30000) {
        return 13;
      } else {
        return 15;
      }
    } else if (additionalDataLength === 4) {
      if (total < 300) {
        return 5;
      } else if (total < 1000) {
        return 8;
      } else if (total < 3000) {
        return 11;
      } else {
        return 15;
      }
    } else if (additionalDataLength === 3) {
      if (total < 300) {
        return 5;
      } else if (total < 1000) {
        return 10;
      } else {
        return 15;
      }
    } else {
      return 10;
    }
  }

  const traces = [
    // trace1,
    // trace2,
    // trace3,
    ...tracesSize,
    ...scatterTraces,
    ...scatterLegdend,
  ];
  let layout = {
    height: 1000,
    width: 1000,
    hoverlabel: { bgcolor: '#FFF' },
    legend: {
      x: 1,
      y: 0.5,
    },
    xaxis: {
      title: pattern2,
      range: [-0.1, 1.1],
      dtick: 0.25,
      zeroline: false,
    },
    yaxis: {
      title: pattern1,
      range: [-0.1, 1.1],
      dtick: 0.25,
      zeroline: false,
    },
    title: {
      text: '<b>Proportion of Mutational Pattern Context Compared to Other Contexts with the same SBS Mutation</b>',
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
      filename: 'Proportion of Mutational Pattern',
    },
  };

  return { traces: traces, layout: layout, config };
}
