export default function SBS24(data, sample) {
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
  const numberWithCommas = (x) =>
    x.toString().replace(/\B(?<!\.\d*)(?=(\d{3})+(?!\d))/g, ',');

  const groupByMutation = data.reduce((groups, e, i) => {
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
  const flatSorted = Object.values(groupByMutation).flat();

  //group data by 1st letter
  const dataT = [];
  const dataU = [];
  Object.entries(flatSorted).forEach(([key, value], groupIndex, array) => {
    if (value.mutationType.substring(0, 1) === 'T') {
      dataT.push(value);
    } else if (value.mutationType.substring(0, 1) === 'U') {
      dataU.push(value);
    }
  });

  const totalMutations = [...dataT, ...dataU].reduce(
    (a, e) => a + parseInt(e.contribution),
    0
  );

  const dataTU = [...dataT, ...dataU];

  const maxVal = Math.max(...dataTU.map((o) => o.contribution));
  //   console.log("dataT");
  //   console.log(dataT);
  //   console.log("dataU");
  //   console.log(dataU);

  const tracesT = {
    name: 'Transcrribed',
    type: 'bar',
    marker: { color: '#004765' },

    x: dataT.map((element, index, array) => element.contribution),
    y: dataU.map(
      (element, index, array) =>
        element.mutationType.substring(2, element.mutationType.length) + ' '
    ),
    hoverinfo: 'y+x',
    orientation: 'h',
  };

  //console.log(tracesT);
  const tracesU = {
    name: 'Untranscribed',
    type: 'bar',
    marker: { color: '#E32925' },
    x: dataU.map((element, index, array) => element.contribution),
    y: dataU.map(
      (element, index, array) =>
        element.mutationType.substring(2, element.mutationType.length) + ' '
    ),
    hoverinfo: 'y+x',
    orientation: 'h',
  };
  //console.log(tracesU);

  const traces = [tracesU, tracesT];

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    bargap: 0.3,
    legend: {
      x: 1,
      xanchor: 'right',
      y: 1,
      traceorder: 'reversed',
    },
    title: {
      text:
        '<b>' + sample + ': ' + numberWithCommas(totalMutations) + ' subs </b>',
      font: {
        size: 24,
      },
      xref: 'paper',
      x: 0.05,
    },
    xaxis: {
      title: {
        text: '<b>Number of Single Base Substitution</b>',
        font: {
          size: 18,
        },
      },
      tickfont: {
        size: 16,
      },
      autorange: false,
      range: [0, maxVal + maxVal * 0.15],
    },
    yaxis: {
      tickfont: {
        size: 16,
      },
      categoryorder: 'category descending',
    },
  };
  // console.log("layout");
  // console.log(layout);

  return { traces, layout };
}
