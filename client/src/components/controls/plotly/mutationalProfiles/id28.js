import { createSampleAnnotation } from './utils';

export default function ID28(data, title = '') {
  const colors = {
    '1:Del:C': '#FBBD6F',
    '1:Del:T': '#FE8002',
    '1:Ins:C': '#AEDD8A',
    '1:Ins:T': '#35A12E',
    'o:': '#1764AA',
  };
  const annotationColors = {
    '1:Del:C': 'black',
    '1:Del:T': 'white',
    '1:Ins:C': 'black',
    '1:Ins:T': 'white',
    'o:': '#1764AA',
  };
  const arrayIDAnnXTop = ['1bp Deletion', '1bp Insertion', '>1bp'],
    arrayIDAnnXBot = ['Homopolymer Length', 'Homopolymer Length', 'Type'],
    arrayIDAnnXLabel = [5, 17, 25],
    arrayIDAnnotationTop = [],
    arrayIDAnnotationBot = [];

  const totalMutations = data.reduce((a, e) => a + parseInt(e.mutations), 0);

  const maxVal = Math.max(...data.map((o) => o.mutations));

  const data1 = data.slice(0, data.length - 4);
  const data2 = data.slice(-4);
  //data2.push(data2.shift());

  // group data by dominant mutation
  const groupByMutation = data1.reduce((groups, e, i) => {
    const mutationRegex = /^.{0,7}/;
    const mutation = e.mutationType.match(mutationRegex)[0];

    const signature = {
      mutationType: e.mutationType,
      contribution: e.mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  //console.log(groupByMutation);

  const groupByFirstGroup = Object.fromEntries(
    Object.entries(groupByMutation).slice(0, 4)
  );

  const arrayID1 = Object.keys(groupByFirstGroup).map(function (key) {
    return groupByFirstGroup[key];
  });

  const arrayID2_Mod = data2.map((element) => ({
    mutationType: 'o:' + element.mutationType,
    contribution: element.mutations,
  }));

  const groupO = arrayID2_Mod.reduce((r, a) => {
    let m;
    m = a.mutationType.match(a.mutationType.substr(0, 1));
    const s = {
      mutationType: a.mutationType,
      contribution: a.contribution,
    };
    r[m] = r[m] ? [...r[m], a] : [s];
    return r;
  }, {});

  const arrayID2 = Object.keys(groupO).map(function (key) {
    return groupO[key];
  });
  const arrayID = [...arrayID1, ...arrayID2];
  const flatSorted = Object.values(arrayID).flat();

  Object.values(arrayID).forEach((group) => {
    if (group.length > 1) {
      arrayIDAnnotationTop.push(
        group[Math.floor(group.length / 2)].mutationType
      );
    } else {
      arrayIDAnnotationTop.push(group[0].mutationType);
    }
    group.forEach((e) => {
      //arrayIDAnnotationBot.push(e.mutationType);
      let lastNum = e.mutationType.substring(
        e.mutationType.length - 1,
        e.mutationType.length
      );
      let newNum;
      if (e.mutationType.substring(0, 1) === 'o') {
        newNum = e.mutationType;
      } else {
        if (e.mutationType.substring(2, 5) === 'Del') {
          lastNum = +lastNum + 1;
        }

        if ((e.mutationType.substring(2, 5) === 'Del') & (+lastNum > 5)) {
          newNum = lastNum + '+';
        } else if (
          (e.mutationType.substring(2, 5) !== 'Del') &
          (+lastNum > 4)
        ) {
          newNum = lastNum + '+';
        } else {
          newNum = lastNum + '';
        }
      }

      arrayIDAnnotationBot.push(newNum);
    });
  });

  const traces = Object.entries(arrayID).map(
    ([mutation, signatures], groupIndex, array) => ({
      name: mutation,
      signature: signatures,
      type: 'bar',
      marker: {
        color:
          signatures[0].mutationType.substring(0, 1) === 'o'
            ? '#1764AA'
            : colors[
                signatures[0].mutationType.substring(
                  0,
                  signatures[0].mutationType.length - 2
                )
              ],
      },
      x: signatures.map(
        (e, i) =>
          array
            .slice(0, groupIndex)
            .reduce((x0, [_, sigs]) => x0 + sigs.length, 0) + i
      ),
      y: signatures.map((e) => e.contribution),
      customdata: signatures.map((e) => ({
        mutationType:
          e.mutationType.substring(2, 5) === 'Del' ? 'Deletion' : 'Insertion',
        xval:
          e.mutationType.substring(2, 5) === 'Del'
            ? +e.mutationType.slice(-1) + 1
            : e.mutationType.slice(-1),
      })),
      hovertemplate: '%{y} indels<extra></extra>',
      //hoverinfo: 'x+y',
      showlegend: false,
    })
  );

  const annotations1 = Object.entries(arrayID).map(
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
      y: 1.01,
      text:
        groupIndex < 4
          ? `<b>${signatures[0].mutationType.substring(
              signatures[0].mutationType.length - 3,
              signatures[0].mutationType.length - 2
            )}</b>`
          : `<b>${signatures[0].mutationType.substring(0, 1)}</b>`,
      showarrow: false,
      font: {
        size: 14,
        color:
          signatures[0].mutationType.substring(0, 1) === 'o'
            ? '#1764AA'
            : annotationColors[
                signatures[0].mutationType.substring(
                  0,
                  signatures[0].mutationType.length - 2
                )
              ],
      },
      align: 'center',
      signatures: signatures,
      mutation: mutation,
      groupIndex: groupIndex,
    })
  );
  const annotations2 = arrayIDAnnotationBot.map((num, index) => ({
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: index,
    y: num.substring(0, 1) === 'o' ? -0.12 : -0.1,
    xnum: parseInt(num.substring(num.length - 1, num.length)),
    text:
      num.substring(0, 1) === 'o'
        ? num === 'o:complex'
          ? '<b>comp</b>'
          : num === 'o:MH'
          ? '<b>MH</b>'
          : '<b>' + num.substring(num.length - 3, num.length) + '</b>'
        : '<b>' + num + '</b>',
    showarrow: false,
    font: {
      size: num.substring(0, 1) === 'o' ? 11 : 14,
    },
    textangle: num.substring(0, 1) === 'o' ? -90 : 0,
    align: 'center',
    num: num,
    index: index,
  }));

  const annotationsIDTopLabel = arrayIDAnnXLabel.map((num, index) => ({
    xref: 'x',
    yref: 'paper',
    x: num,
    xanchor: 'bottom',
    y: 1.07,
    yanchor: 'bottom',
    text: '<b>' + arrayIDAnnXTop[index] + '</b>',
    showarrow: false,
    font: {
      size: 14,
    },
    align: 'center',
  }));

  const annotationsIDBotLabel = arrayIDAnnXLabel.map((num, index) => ({
    xref: 'x',
    yref: 'paper',
    x: num,
    xanchor: 'bottom',
    y: -0.16,
    yanchor: 'bottom',
    text: '<b>' + arrayIDAnnXBot[index] + '</b>',
    showarrow: false,
    font: {
      size: 16,
      family: 'Times New Roman',
    },
    align: 'center',
  }));

  const sampleAnnotation = createSampleAnnotation(data);

  const shapes1 = Object.entries(arrayID).map(
    ([mutation, signatures], groupIndex, array) => ({
      type: 'rect',
      xref: 'x',
      yref: 'paper',
      x0: array
        .slice(0, groupIndex)
        .reduce((x0, [_, sigs]) => x0 + sigs.length, -0.4),
      x1: array
        .slice(0, groupIndex + 1)
        .reduce((x0, [_, sigs]) => x0 + sigs.length, -0.6),
      y0: 1.07,
      y1: 1.01,
      fillcolor:
        signatures[0].mutationType.substring(0, 1) === 'o'
          ? '#1764AA'
          : colors[
              signatures[0].mutationType.substring(
                0,
                signatures[0].mutationType.length - 2
              )
            ],
      line: {
        width: 0,
      },
      mutation: mutation,
      signature: signatures[0].mutationType.substring(
        0,
        signatures[0].mutationType.length - 2
      ),
    })
  );

  const shapes2 = Object.entries(arrayID).map(
    ([mutation, signatures], groupIndex, array) => ({
      type: 'rect',
      xref: 'x',
      yref: 'paper',
      x0: array
        .slice(0, groupIndex)
        .reduce((x0, [_, sigs]) => x0 + sigs.length, -0.4),
      x1: array
        .slice(0, groupIndex + 1)
        .reduce((x0, [_, sigs]) => x0 + sigs.length, -0.6),
      y0: -0.01,
      y1: -0.05,
      fillcolor:
        signatures[0].mutationType.substring(0, 1) === 'o'
          ? '#1764AA'
          : colors[
              signatures[0].mutationType.substring(
                0,
                signatures[0].mutationType.length - 2
              )
            ],
      line: {
        width: 0,
      },
      mutation: mutation,
      signature: signatures[0].mutationType.substring(
        0,
        signatures[0].mutationType.length - 2
      ),
    })
  );

  const layout = {
    title: `<b>${title}</b>`,
    hoverlabel: { bgcolor: '#FFF' },
    height: 600,
    width: 750,
    xaxis: {
      showticklabels: false,
      showline: true,
      tickangle: -90,
      tickfont: {
        size: 10,
      },
      tickmode: 'array',
      tickvals: flatSorted.map((_, i) => i),
      ticktext: flatSorted.map((e) => e.mutationType),
      linecolor: 'black',
      linewidth: 1,
      mirror: true,
    },
    yaxis: {
      title: {
        text: '<b>Number of Indels</b>',
        font: {
          family: 'Times New Roman',
          size: 20,
        },
      },
      autorange: false,
      range: [0, maxVal + maxVal * 0.25],
      linecolor: 'black',
      linewidth: 1,
      mirror: true,
    },

    shapes: [...shapes1, ...shapes2],
    annotations: [
      ...annotations1,
      ...annotations2,
      ...annotationsIDTopLabel,
      ...annotationsIDBotLabel,
      sampleAnnotation,
    ],
  };

  return { traces, layout };
}
