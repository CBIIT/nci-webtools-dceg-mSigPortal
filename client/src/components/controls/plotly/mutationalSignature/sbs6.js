export default function SBS6(data, sample) {
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

  const totalMutations = data.reduce((a, e) => a + parseInt(e.mutations), 0);

  // group data by dominant mutation
  const groupByMutation = data.reduce((groups, e, i) => {
    const mutation = e.mutationType.substring(0, e.mutationType.length);
    const signature = {
      mutationType: e.mutationType,
      contribution: e.mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  const traces = Object.entries(groupByMutation).map(
    ([mutation, signatures], groupIndex, array) => ({
      name: mutation,
      type: 'bar',
      marker: { color: colors[mutation] },
      y: signatures.map((e) => e.mutationType + ' '),
      x: signatures.map((e) => e.contribution),
      hoverinfo: 'x+y',
      orientation: 'h',
    })
  );

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    showlegend: false,
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
    },
    yaxis: {
      tickfont: {
        size: 16,
      },
      categoryorder: 'category descending',
    },
  };
  return { traces, layout };
}
