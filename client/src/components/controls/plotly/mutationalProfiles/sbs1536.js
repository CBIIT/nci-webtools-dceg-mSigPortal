import { sbsColor } from '../../utils/colors.js';
import { createSampleAnnotation } from './utils.js';
export default function SBS1536(data, title = '') {
  const colors = sbsColor;

  const heatmapColorscale = [
    [0, 'rgb(56,56,156'],
    [0.2, 'rgb(56,56,156'],
    [0.2, 'rgb(106,106,128'],
    [0.4, 'rgb(106,106,128'],
    [0.4, 'rgb(155,146,98'],
    [0.6, 'rgb(155,146,98'],
    [0.6, 'rgb(205,186,69'],
    [0.8, 'rgb(205,186,69'],
    [0.8, 'rgb(255,255,39)'],
    [1, 'rgb(255,255,39)'],
  ];

  const totalMutations = data.reduce((a, e) => a + parseInt(e.mutations), 0);
  const chunks = (a, size) =>
    Array.from(new Array(Math.ceil(a.length / size)), (_, i) =>
      a.slice(i * size, i * size + size)
    );

  // const maxValMutation = Math.max(...data.map((o) => o.mutations));
  // console.log("maxValMutation:---");
  // console.log(maxValMutation);

  const groupByMutationInner = data.reduce((groups, e, i) => {
    const mutation = e.mutationType.substring(1, 8);
    const signature = {
      mutationType: e.mutationType,
      contribution: e.mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  const groupByMutationOuter = data.reduce((groups, e, i) => {
    const mutation =
      e.mutationType[0] + e.mutationType[e.mutationType.length - 1];
    const signature = {
      mutationType: e.mutationType,
      contribution: e.mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  /////---------------- Bar Chart ---------------------------------////

  const totalMutationsGroup = Object.entries(groupByMutationInner).map(
    ([mutation, signatures], groupIndex, array) => ({
      mutationType: mutation,
      signatures: signatures,
      total: signatures.reduce((a, e) => a + parseInt(e.contribution), 0),
    })
  );

  const groupByTotal = totalMutationsGroup.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.mutationType.match(mutationRegex)[1];
    const signature = {
      mutationType: e.mutationType,
      contribution: e.total,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  const flatSorted = Object.values(groupByTotal).flat();
  //const mutationTitle = Object.keys(groupByTotal).flat();

  const maxVal = Math.max(...flatSorted.map((o) => o.contribution));

  const tracesBar = Object.entries(groupByTotal).map(
    ([mutation, signatures], groupIndex, array) => ({
      name: mutation,
      type: 'bar',
      marker: { color: colors[mutation] },
      //   x: signatures.map((e) => e.mutationType),
      //x: signatures.map((e, i) => groupIndex * signatures.length + i),
      x: signatures.map(
        (e, i) =>
          array
            .slice(0, groupIndex)
            .reduce((x0, [_, sigs]) => x0 + sigs.length, 0) + i
      ),
      y: signatures.map((e) => e.contribution),
      hoverinfo: 'x+y',
      showlegend: false,
    })
  );

  ////// ------- Heat Map 1 -----////
  const heatmapY2 = [];
  const heatmapZ2 = [];
  const heatmapX2 = [];

  const groupByMutationFront = data.reduce((groups, e, i) => {
    const mutation = e.mutationType.substring(0, e.mutationType.length - 1);
    const signature = {
      mutationType: e.mutationType,
      contribution: e.mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  const mutationSumFront = Object.entries(groupByMutationFront).map(
    ([key, value]) => ({
      mutationType: key,
      contribution: value.reduce((a, e) => a + parseInt(e.contribution), 0),
    })
  );

  const arrayMutationSumFront = Object.values(mutationSumFront).flat();
  const groupByMutation2 = arrayMutationSumFront.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.mutationType.match(mutationRegex)[1];
    const signature = {
      mutationType: e.mutationType,
      contribution: e.contribution,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  Object.entries(groupByMutation2).forEach(
    ([key, value], groupIndex, array) => {
      heatmapY2.push(Object.entries(value).map(([k, v]) => v.mutationType));
      heatmapZ2.push(
        Object.entries(value).map(([k, v]) => v.contribution / totalMutations)
      );
      heatmapX2.push(
        value.map(
          (e, i) =>
            array
              .slice(0, groupIndex)
              .reduce((x0, [_, sigs]) => x0 + sigs.length, 0) + i
        )
      );
    }
  );

  let heatmapY2_c = [
    heatmapY2[0][0].charAt(0) + '--N',
    heatmapY2[0][16].charAt(0) + '--N',
    heatmapY2[0][32].charAt(0) + '--N',
    heatmapY2[0][48].charAt(0) + '--N',
  ];

  let heatMapZ2_0 = chunks(heatmapZ2[0], 16);
  let heatMapZ2_1 = chunks(heatmapZ2[1], 16);
  let heatMapZ2_2 = chunks(heatmapZ2[2], 16);
  let heatMapZ2_3 = chunks(heatmapZ2[3], 16);
  let heatMapZ2_4 = chunks(heatmapZ2[4], 16);
  let heatMapZ2_5 = chunks(heatmapZ2[5], 16);

  const heatMapZFinal2 = [
    heatMapZ2_0,
    heatMapZ2_1,
    heatMapZ2_2,
    heatMapZ2_3,
    heatMapZ2_4,
    heatMapZ2_5,
  ];
  const maxZ2 = Math.max(...heatMapZFinal2.flat(Infinity));
  const traceHeatMap2 = heatMapZFinal2.map((num, index, array) => ({
    colorbar: { len: 0.2, y: 0.625 },
    colorscale: heatmapColorscale,
    zmin: 0,
    zmax: maxZ2 + maxZ2 * 0.1,
    z: num,
    y: heatmapY2_c,
    type: 'heatmap',
    hoverongaps: false,
    xaxis: 'x',
    yaxis: 'y2',
    x: num.map(
      (e, i) =>
        array.slice(0, index).reduce((x0, [_, sigs]) => x0 + sigs.length, 0) + i
    ),
    xgap: 0.1,
    ygap: 0.1,
    hovertemplate: 'x: %{x}<br>y: %{y}<br>Value: %{z}<extra></extra>',
  }));

  ////// ------------------- Heat Map 2 --------------------------------////

  const heatmapY3 = [];
  const heatmapZ3 = [];
  const heatmapX3 = [];
  const groupByMutationBack = data.reduce((groups, e, i) => {
    const mutation = e.mutationType.substring(1, e.mutationType.length);
    const signature = {
      mutationType: e.mutationType,
      contribution: e.mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  // console.log("groupByMutationBack:---");
  // console.log(groupByMutationBack);

  const mutationSumBack = Object.entries(groupByMutationBack).map(
    ([key, value]) => ({
      mutationType: key,
      contribution: value.reduce((a, e) => a + parseInt(e.contribution), 0),
    })
  );

  // console.log("mutationsumBack:---");
  // console.log(mutationsumBack);

  // sort by the last letter
  mutationSumBack.sort((a, b) =>
    a.mutationType.substring(a.mutationType.length - 1, a.mutationType.length) <
    b.mutationType.substring(b.mutationType.length - 1, b.mutationType.length)
      ? -1
      : b.mutationType.substring(
          b.mutationType.length - 1,
          b.mutationType.length
        ) <
        a.mutationType.substring(
          a.mutationType.length - 1,
          a.mutationType.length
        )
      ? 1
      : 0
  );

  const arrayMutationSumBack = Object.values(mutationSumBack).flat();

  const groupByMutation3 = arrayMutationSumBack.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.mutationType.match(mutationRegex)[1];
    const signature = {
      mutationType: e.mutationType,
      contribution: e.contribution,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  Object.entries(groupByMutation3).forEach(
    ([key, value], groupIndex, array) => {
      heatmapY3.push(Object.entries(value).map(([k, v]) => v.mutationType));
      heatmapZ3.push(
        Object.entries(value).map(([k, v]) => v.contribution / totalMutations)
      );
      heatmapX3.push(
        value.map(
          (e, i) =>
            array
              .slice(0, groupIndex)
              .reduce((x0, [_, sigs]) => x0 + sigs.length, 0) + i
        )
      );
    }
  );

  let heatmapY3_c = [
    'N--' + heatmapY3[0][0].charAt(heatmapY3[0][0].length - 1),
    'N--' + heatmapY3[0][16].charAt(heatmapY3[0][16].length - 1),
    'N--' + heatmapY3[0][32].charAt(heatmapY3[0][32].length - 1),
    'N--' + heatmapY3[0][48].charAt(heatmapY3[0][48].length - 1),
  ];

  let heatMapZ3_0 = chunks(heatmapZ3[0], 16);
  let heatMapZ3_1 = chunks(heatmapZ3[1], 16);
  let heatMapZ3_2 = chunks(heatmapZ3[2], 16);
  let heatMapZ3_3 = chunks(heatmapZ3[3], 16);
  let heatMapZ3_4 = chunks(heatmapZ3[4], 16);
  let heatMapZ3_5 = chunks(heatmapZ3[5], 16);

  const heatMapZFinal3 = [
    heatMapZ3_0,
    heatMapZ3_1,
    heatMapZ3_2,
    heatMapZ3_3,
    heatMapZ3_4,
    heatMapZ3_5,
  ];

  const maxZ3 = Math.max(...heatMapZFinal2.flat(Infinity));
  const traceHeatMap3 = heatMapZFinal3.map((num, index, array) => ({
    colorbar: { len: 0.2, y: 0.44 },
    colorscale: heatmapColorscale,
    zmin: 0,
    zmax: maxZ3 + maxZ3 * 0.1,
    z: num,
    y: heatmapY3_c,
    type: 'heatmap',
    hoverongaps: false,
    xaxis: 'x',
    yaxis: 'y3',
    x: num.map(
      (e, i) =>
        array.slice(0, index).reduce((x0, [_, sigs]) => x0 + sigs.length, 0) + i
    ),
    xgap: 0.1,
    ygap: 0.1,
    hovertemplate: 'x: %{x}<br>y: %{y}<br>Value: %{z}<extra></extra>',
  }));
  ////--------------------- Heat Map Total --------------------------//
  const heatmapY = [];
  const heatmapZ = [];
  const heatmapX = [];

  let heatMapZ0 = [];
  let heatMapZ1 = [];
  let heatMapZ2 = [];
  let heatMapZ3 = [];
  let heatMapZ4 = [];
  let heatMapZ5 = [];

  Object.entries(groupByMutationOuter).forEach(
    ([key, value], groupIndex, array) => {
      value.sort((a, b) =>
        a.mutationType.substring(3, 6) < b.mutationType.substring(3, 6)
          ? -1
          : b.mutationType.substring(3, 6) < a.mutationType.substring(3, 6)
          ? 1
          : 0
      );
      //console.log(value);
      heatmapY.push(key.charAt(0) + '--' + key.charAt(key.length - 1));

      //console.log(totalMutations);
      //console.log(value);
      //console.log(Object.entries(value).map(([k, v]) => v.contribution));
      heatmapZ.push(
        Object.entries(value).map(([k, v]) => v.contribution / totalMutations)
      );
      heatmapX.push(
        value.map(
          (e, i) =>
            array
              .slice(0, groupIndex)
              .reduce((x0, [_, sigs]) => x0 + sigs.length, 0) + i
        )
      );
    }
  );

  heatmapZ.forEach((item, index) => {
    heatMapZ0.push(item.slice().splice(0, 16));
    heatMapZ1.push(item.slice().splice(16, 16));
    heatMapZ2.push(item.slice().splice(32, 16));
    heatMapZ3.push(item.slice().splice(48, 16));
    heatMapZ4.push(item.slice().splice(64, 16));
    heatMapZ5.push(item.slice().splice(80, 16));
  });

  const heatMapZFinal = [
    heatMapZ0,
    heatMapZ1,
    heatMapZ2,
    heatMapZ3,
    heatMapZ4,
    heatMapZ5,
  ];

  const maxZ = Math.max(...heatMapZFinal.flat(Infinity));
  const traceHeatMap = heatMapZFinal.map((num, index, array) => ({
    colorbar: { len: 0.38, y: 0.17 },
    colorscale: heatmapColorscale,
    zmin: 0,
    zmax: maxZ + maxZ * 0.1,
    z: num,
    y: heatmapY,

    type: 'heatmap',
    hoverongaps: false,
    xaxis: 'x',
    yaxis: 'y4',
    x: num.map(
      (e, i) =>
        array.slice(0, index).reduce((x0, [_, sigs]) => x0 + sigs.length, 0) + i
    ),
    test: heatmapY.map((a) => a.replace('--', `%{x}`)),
    array: array,
    num: num,
    xgap: 0.1,
    ygap: 0.1,
    hovertemplate: 'x: %{x}<br>y: %{y}<br>Value: %{z}<extra></extra>',
  }));

  const traces = [
    ...tracesBar,
    ...traceHeatMap,
    ...traceHeatMap2,
    ...traceHeatMap3,
  ];

  //console.log('traces:');
  //console.log(traces);
  const annotations = Object.entries(groupByTotal).map(
    ([mutation, signatures], groupIndex, array) => ({
      xref: 'x',
      yref: 'paper',
      xanchor: 'bottom',
      yanchor: 'bottom',
      x:
        array
          .slice(0, groupIndex)
          .reduce((x0, [_, sigs]) => x0 + sigs.length, 0) +
        (signatures.length - 1) * 0.5,
      y: 1.04,
      text: `<b>${mutation}</b>`,
      showarrow: false,
      font: {
        size: 18,
      },
      align: 'center',
    })
  );

  const sampleAnnotation = createSampleAnnotation(data, false, 0.95);

  const shapes = Object.entries(groupByTotal).map(
    ([mutation, _], groupIndex, array) => ({
      type: 'rect',
      xref: 'x',
      yref: 'paper',
      x0: array
        .slice(0, groupIndex)
        .reduce((x0, [_, sigs]) => x0 + sigs.length, -0.4),
      x1: array
        .slice(0, groupIndex + 1)
        .reduce((x0, [_, sigs]) => x0 + sigs.length, -0.6),
      y0: 1.04,
      y1: 1.01,
      fillcolor: colors[mutation],
      line: {
        width: 0,
      },
    })
  );

  const xannotations = flatSorted.map((num, index) => ({
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: index,
    y: -0.04,
    text: num.mutationType.replace(/\[(.*)\]/, '-'),
    showarrow: false,
    font: {
      size: 7.5,
      family: 'Courier New, monospace',
    },
    align: 'center',
    num: num,
    index: index,
    textangle: -90,
  }));

  const yLabelAnnotation = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'top',
    yanchor: 'top',
    x: -0.045,
    y: 1.02,
    text: '<b>Number of Single Base Substitutions</b>',
    showarrow: false,
    font: {
      size: 10,
      family: 'Times New Roman',
    },
    align: 'center',
    textangle: -90,
  };

  const layout = {
    title: `<b>${title}</b>`,
    hoverlabel: { bgcolor: '#FFF' },
    height: 800,
    width: 1080,
    grid: {
      rows: 4,
      columns: 1,
    },
    xaxis: {
      showticklabels: false,
      showline: true,
      tickangle: -90,
      tickfont: {
        family: 'Courier New, monospace',
        size: 8,
      },
      tickmode: 'array',
      tickvals: flatSorted.map((_, i) => i),
      ticktext: flatSorted.map((e) => e.mutationType),
      linecolor: '#D3D3D3',
      linewidth: 1,
      mirror: 'all',
      tickformat: maxVal > 1000 ? '~s' : '',
      ticks: '',
    },
    yaxis: {
      //title: 'Number of Single Base Substitutions',
      autorange: false,
      range: [0, maxVal + maxVal * 0.2],
      linecolor: '#D3D3D3',
      linewidth: 1,
      mirror: 'all',
      domain: [0.72, 1],
      tickformat: maxVal > 1000 ? '~s' : '',
      tickfont: {
        family: 'Courier New, monospace',
        size: 8,
      },
      showgrid: true,
      gridcolor: '#F5F5F5',
    },
    yaxis2: {
      autorange: true,
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      anchor: 'x',
      domain: [0.54, 0.715],
      tickfont: {
        family: 'Courier New, monospace',
        size: 8,
      },
    },
    yaxis3: {
      autorange: true,
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      anchor: 'x',
      tickfont: {
        family: 'Courier New, monospace',
        size: 8,
      },
      domain: [0.36, 0.535],
    },
    yaxis4: {
      autorange: true,
      linecolor: '#D3D3D3',
      linewidth: 1,
      ticks: '',
      mirror: 'all',
      anchor: 'x',
      //dtick: 1,
      tickfont: {
        family: 'Courier New, monospace',
        size: 8,
      },
      domain: [0, 0.35],
    },

    // autosize: false,
    shapes: shapes,
    annotations: [
      ...annotations,
      sampleAnnotation,
      ...xannotations,
      yLabelAnnotation,
    ],
  };
  // console.log("layout");
  //console.log(layout);

  return { traces, layout };
}
