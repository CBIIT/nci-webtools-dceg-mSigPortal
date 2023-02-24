import { groupBy } from 'lodash';

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
    legendgroup: 'size',
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
    legendgroup: 'size',
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
    legendgroup: 'size',
    legendgrouptitle: {
      text: 'Number of mutations',
    },
  };

  let trace4 = {
    name:
      result[0]?.study != 'undefined'
        ? result[0].study + '@' + result[0].cancer
        : 'Input',
    x: result.map((e) => e.n1),
    y: result.map((e) => e.n2),
    customdata: result.map((e) => ({ sample: e.sample, total: e.total })),
    mode: 'markers',
    type: 'scatter',
    opacity: 1,
    marker: {
      color: 'green',
      line: {
        color: 'black',
        width: 1,
      },
      size: result.map((e, i, a) =>
        e.total < maxMutationFilter / 100
          ? 5
          : e.total < maxMutationFilter / 10
          ? 10
          : 15
      ),
    },
    hovertemplate:
      '<b>Sample: ' +
      '</b>' +
      '%{customdata.sample} <br>' +
      '<b>' +
      pattern2 +
      ':</b>' +
      ' %{x} <br>' +
      '<b>' +
      pattern1 +
      ':</b>' +
      ' %{y}<br>' +
      '<b>Total:</b> %{customdata.total}<extra></extra>',
    showlegend: true,
    legendgrouptitle: {
      text: 'Study',
    },
  };

  const traces = [trace1, trace2, trace3, trace4];
  let layout = {
    height: 700,
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
