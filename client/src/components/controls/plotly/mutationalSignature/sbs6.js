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
  const maxVal = Math.max(...data.map((o) => o.Mutations));
  // console.log("maxVal--:");
  // console.log(maxVal);

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
  const flatSorted = Object.values(groupByMutation).flat();

  console.log("groupByMutation:");
  console.log(groupByMutation);
  console.log("FlatSorted");
  console.log(flatSorted);

  const traces = Object.entries(groupByMutation).map(
    ([mutation, signatures], groupIndex, array) => ({
      name: mutation,
      type: "bar",
      marker: { color: colors[mutation] },
      y: signatures.map((e) => `<b>${e.mutationType}<b>`),
      x: signatures.map((e) => e.contribution),
      hoverinfo: "x+y",
      showlegend: false,
      orientation: "h",
    })
  );
  //const traces = [trace1, trace2];
  console.log("traces:");
  console.log(traces);

  const layout = {
    hoverlabel: { bgcolor: "#FFF" },
    title: {
      title: {
        text:
          "<b>" +
          sample +
          ": " +
          numberWithCommas(totalMutations) +
          " subs </b>",
        font: { size: 24 },
      },
    },
    xaxis: {
      title: "<b>Number of Single Base Substitutions</b>",
      showticklabels: true,
      showline: true,
      tickfont: {
        size: 10,
      },
    },
    yaxis: {
      categoryorder: "category descending",
    },
  };
  console.log("layout");
  console.log(layout);

  return { traces, layout };
}
