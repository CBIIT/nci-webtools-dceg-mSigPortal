export default function SBS6(data, sample) {
  const colors = {
    "C>A": "#03BCEE",
    "C>G": "black",
    "C>T": "#E32926",
    "T>A": "#CAC9C9",
    "T>C": "#A1CE63",
    "T>G": "#EBC6C4",
  };
  console.log("data--:");
  console.log(data);
  const numberWithCommas = (x) =>
    x.toString().replace(/\B(?<!\.\d*)(?=(\d{3})+(?!\d))/g, ",");

  const totalMutations = data.reduce((a, e) => a + parseInt(e.Mutations), 0);

  // group data by dominant mutation
  const groupByMutation = data.reduce((groups, e, i) => {
    const mutation = e.MutationType.substring(0, e.MutationType.length);
    const signature = {
      mutationType: e.MutationType,
      contribution: e.Mutations,
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  const traces = Object.entries(groupByMutation).map(
    ([mutation, signatures], groupIndex, array) => ({
      name: mutation,
      type: "bar",
      marker: { color: colors[mutation] },
      y: signatures.map((e) => e.mutationType + " "),
      x: signatures.map((e) => e.contribution),
      hoverinfo: "x+y",
      orientation: "h",
    })
  );
  const sampleAnnotation = {
    xref: "paper",
    yref: "paper",
    xanchor: "bottom",
    yanchor: "bottom",
    x: 0,
    y: 0.9,
    text:
      "<b>" + sample + ": " + numberWithCommas(totalMutations) + " subs </b>",
    showarrow: false,
    font: {
      size: 24,
    },
    align: "center",
  };

  const layout = {
    hoverlabel: { bgcolor: "#FFF" },
    showlegend: false,
    title: {
      text:
        "<b>" + sample + ": " + numberWithCommas(totalMutations) + " subs </b>",
      font: {
        size: 24,
      },
      xref: "paper",
      x: 0.05,
    },
    xaxis: {
      title: {
        text: "<b>Number of Single Base Substitution</b>",
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
      categoryorder: "category descending",
    },

    annotations: sampleAnnotation,
  };

  console.log("layout");
  console.log(layout);
  return { traces, layout };
}
