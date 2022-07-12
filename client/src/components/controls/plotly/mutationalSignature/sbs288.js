export default function SBS288(data, sample) {
  const colors = {
    'C>A': '#03BCEE',
    'C>G': 'black',
    'C>T': '#E32926',
    'T>A': '#CAC9C9',
    'T>C': '#A1CE63',
    'T>G': '#EBC6C4',
  };
  //   console.log("data--:");
  //   console.log(data);

  const arrayDataT = [];
  const arrayDataU = [];
  const arrayDataN = [];

  Object.values(data).forEach((group) => {
    if (group.mutationType.substring(0, 1) === 'T') {
      arrayDataT.push(group);
    } else if (group.mutationType.substring(0, 1) === 'U') {
      arrayDataU.push(group);
    } else {
      arrayDataN.push(group);
    }
  });

  //   console.log("---- arrayDataT ---- ");
  //   console.log(arrayDataT);
  //   console.log("---- arrayDataU ---- ");
  //   console.log(arrayDataU);
  //   console.log("---- arrayDataN ---- ");
  //   console.log(arrayDataN);

  const numberWithCommas = (x) =>
    x.toString().replace(/\B(?<!\.\d*)(?=(\d{3})+(?!\d))/g, ',');

  const totalMutations = data.reduce((a, e) => a + parseInt(e.mutations), 0);
  const maxVal = Math.max(...data.map((o) => o.mutations));
  //   console.log("maxVal--:");
  //   console.log(maxVal);

  //// ------ bar char left  --------- //

  const groupByMutationWoFirstLetter = data.reduce((groups, e, i) => {
    const mutation = e.mutationType.substring(2, e.mutationType.length);
    const signature = {
      mutationType: e.mutationType,
      contribution: e.mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});
  //   console.log("groupByMutationWoFirstLetter");
  //   console.log(groupByMutationWoFirstLetter);

  const totalMutationsGroup = Object.entries(groupByMutationWoFirstLetter).map(
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
  //   console.log("groupByTotal");
  //   console.log(groupByTotal);

  const flatSortedTotal = Object.values(groupByTotal).flat();
  const maxValTotal = Math.max(...flatSortedTotal.map((o) => o.contribution));

  const tracesBarTotal = Object.entries(groupByTotal).map(
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

  //////------------- bar chart right ---------------//////

  const groupByMutationT = arrayDataT.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.mutationType.match(mutationRegex)[1];
    const signature = {
      mutationType: e.mutationType.substring(2, e.mutationType.length),
      contribution: e.mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];

    return groups;
  }, {});

  //   console.log("groupByMutationT:");
  //   console.log(groupByMutationT);

  const totalMutationsGroupT = Object.entries(groupByMutationT).map(
    ([mutation, signatures], groupIndex, array) => ({
      mutationType: mutation,
      signatures: signatures,
      total: signatures.reduce((a, e) => a + parseInt(e.contribution), 0),
    })
  );

  //   console.log("totalMutationsGroupT:");
  //   console.log(totalMutationsGroupT);

  const groupByMutationU = arrayDataU.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.mutationType.match(mutationRegex)[1];
    const signature = {
      mutationType: e.mutationType.substring(2, e.mutationType.length),
      contribution: e.mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  //   console.log("groupByMutationU:");
  //   console.log(groupByMutationU);

  const totalMutationsGroupU = Object.entries(groupByMutationU).map(
    ([mutation, signatures], groupIndex, array) => ({
      mutationType: mutation,
      signatures: signatures,
      total: signatures.reduce((a, e) => a + parseInt(e.contribution), 0),
    })
  );

  //   console.log("totalMutationsGroupU:");
  //   console.log(totalMutationsGroupU);

  const groupByMutationN = arrayDataN.reduce((groups, e, i) => {
    const mutationRegex = /\[(.*)\]/;
    const mutation = e.mutationType.match(mutationRegex)[1];
    const signature = {
      mutationType: e.mutationType.substring(2, e.mutationType.length),
      contribution: e.mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});
  //   console.log("groupByMutationN:");
  //   console.log(groupByMutationN);

  const totalMutationsGroupN = Object.entries(groupByMutationN).map(
    ([mutation, signatures], groupIndex, array) => ({
      mutationType: mutation,
      signatures: signatures,
      total: signatures.reduce((a, e) => a + parseInt(e.contribution), 0),
    })
  );

  //   console.log("totalMutationsGroupN:");
  //   console.log(totalMutationsGroupN);

  const groupByFirstLetter = data.reduce((groups, e, i) => {
    const mutation = e.mutationType.substring(0, 1);
    const signature = {
      mutationType: e.mutationType,
      contribution: e.mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];

    return groups;
  }, {});

  //   console.log("groupByFirstLetter:");
  //   console.log(groupByFirstLetter);

  const totalGroupByFirstLetter = Object.entries(groupByFirstLetter).map(
    ([mutation, signatures], groupIndex, array) => ({
      mutationType: mutation,
      signatures: signatures,
      total: signatures.reduce((a, e) => a + parseInt(e.contribution), 0),
    })
  );

  //   console.log("totalGroupByFirstLetter:");
  //   console.log(totalGroupByFirstLetter);

  const flatSortedT = Object.values(totalMutationsGroupT).flat();
  const flatSortedU = Object.values(totalMutationsGroupU).flat();
  const flatSortedN = Object.values(totalMutationsGroupN).flat();
  const flatSortedFirstLetter = Object.values(totalGroupByFirstLetter).flat();

  //   console.log(flatSortedT);
  //   console.log(flatSortedU);
  //   console.log(flatSortedN);
  //   console.log(flatSortedFirstLetter);

  Object.entries(flatSortedFirstLetter).forEach(
    ([key, value], groupIndex, array) => {
      if (value.mutationType === 'T') {
        value.mutationType = 'All';
        flatSortedT.push(value);
      } else if (value.mutationType === 'U') {
        value.mutationType = 'All';
        flatSortedU.push(value);
      } else {
        value.mutationType = 'All';
        flatSortedN.push(value);
      }
    }
  );

  // sort by the first letter
  //   flatSortedT.sort((a, b) =>
  //     a.mutationType.substring(0, 1) < b.mutationType.substring(0, 1)
  //       ? -1
  //       : b.mutationType.substring(0, 1) < a.mutationType.substring(0, 1)
  //       ? 1
  //       : 0
  //   );

  //   console.log("New---");
  //   console.log(flatSortedT);
  //   console.log(flatSortedU);
  //   console.log(flatSortedN);
  const tracesT = {
    name: 'Transcrribed',
    type: 'bar',
    marker: { color: '#004765' },

    x: flatSortedT.map((element, index, array) => element.total),
    y: flatSortedT.map(
      (element, index, array) => `<b>${element.mutationType}<b>`
    ),
    xaxis: 'x2',
    yaxis: 'y2',
    hoverinfo: 'x2+y2',
    showlegend: true,
    orientation: 'h',
  };

  const tracesU = {
    name: 'Untranscribed',
    type: 'bar',
    marker: { color: '#E32925' },
    x: flatSortedU.map((element, index, array) => element.total),
    y: flatSortedU.map(
      (element, index, array) => `<b>${element.mutationType}<b>`
    ),
    xaxis: 'x2',
    yaxis: 'y2',
    hoverinfo: 'x2+y2',
    showlegend: true,
    orientation: 'h',
  };

  const tracesN = {
    name: 'Nontranscribed',
    type: 'bar',
    marker: { color: '#008001' },
    x: flatSortedN.map((element, index, array) => element.total),
    y: flatSortedN.map(
      (element, index, array) => `<b>${element.mutationType}<b>`
    ),
    xaxis: 'x2',
    yaxis: 'y2',
    hoverinfo: 'x2+y2',
    showlegend: true,
    orientation: 'h',
  };

  const traces = [...tracesBarTotal, tracesN, tracesU, tracesT];
  //   console.log("traces:");
  //   console.log(traces);
  const xannotations = flatSortedTotal.map((num, index) => ({
    xref: 'x',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: index,
    y: -0.1,
    text: num.mutationType.replace(
      /\[(.*)\]/,
      num.mutationType.substring(2, 3)
    ),
    showarrow: false,
    font: {
      size: 10,
      color: colors[num.mutationType.substring(2, 5)],
    },
    align: 'center',
    num: num,
    index: index,
    textangle: -90,
  }));

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

  const sampleAnnotation = {
    xref: 'paper',
    yref: 'paper',
    xanchor: 'bottom',
    yanchor: 'bottom',
    x: 0,
    y: 0.92,
    text:
      '<b>' + sample + ': ' + numberWithCommas(totalMutations) + ' subs </b>',
    showarrow: false,
    font: {
      size: 18,
    },
    align: 'center',
  };

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
      y0: 1.05,
      y1: 1.01,
      fillcolor: colors[mutation],
      line: {
        width: 0,
      },
      mutation: mutation,
    })
  );
  // console.log("shapes");
  // console.log(shapes);

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    bargap: 0.3,
    legend: {
      x: 1,
      xanchor: 'right',
      y: 0,
      traceorder: 'reversed',
    },
    grid: {
      rows: 1,
      columns: 2,
      pattern: 'independent',
    },
    xaxis: {
      showticklabels: false,
      showline: true,
      tickangle: -90,
      tickfont: {
        size: 10,
      },
      tickmode: 'array',
      tickvals: flatSortedTotal.map((_, i) => i),
      // ticktext: flatSorted.map((e) =>
      //   e.mutationType.replace(/\[(.*)\]/, e.mutationType.substring(2, 3))
      // ),
      ticktext: flatSortedTotal.map((e) => e.mutationType),
      linecolor: 'black',
      linewidth: 1,
      mirror: true,
      domain: [0, 0.75],
    },
    yaxis: {
      title: 'Number of Single Base Substitutions',
      autorange: false,
      range: [0, maxValTotal + maxValTotal * 0.15],
      linecolor: 'black',
      linewidth: 1,
      mirror: true,
    },

    xaxis2: {
      showticklabels: true,
      showline: true,
      tickfont: {
        size: 12,
      },
      domain: [0.8, 1],
    },
    yaxis2: {
      showline: true,
      tickangle: 0,
      tickfont: {
        size: 12,
      },
      anchor: 'x2',
      categoryorder: 'category descending',
    },

    shapes: shapes,
    annotations: [...annotations, sampleAnnotation, ...xannotations],
  };
  //   console.log("layout");
  //   console.log(layout);

  return { traces, layout };
}
