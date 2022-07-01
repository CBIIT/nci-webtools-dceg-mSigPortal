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

  const xval = [];
  const yval = [];
  data.forEach((array) => {
    yval.push(array.MutationType);
    xval.push(array.Mutations);
  });
  //   console.log(xval);
  //   console.log(yval);
  //   const traces = data.map((element, index) => ({
  //     element: element,
  //     idex: index,
  //     name: element.MutationType,
  //     type: "bar",
  //     //marker: { color: colors[element.MutationType] },
  //     x: xval.push(element.MutationType),
  //     y: yval.push(element.Mutations),
  //     hoverinfo: "x+y",
  //   }));
  //   console.log("traces:");
  //   console.log(traces);
  const trace1 = {};
  const trace2 = {
    type: "bar",
    x: xval,
    y: yval,
    hoverinfo: "x+y",
    orientation: "h",
  };
  const traces = [trace1, trace2];
  console.log("traces:");
  console.log(traces);

  //   const annotations = Object.entries(groupByMutation).map(
  //     ([mutation, signatures], groupIndex, array) => ({
  //       xref: "x",
  //       yref: "paper",
  //       xanchor: "bottom",
  //       yanchor: "bottom",
  //       x:
  //         array
  //           .slice(0, groupIndex)
  //           .reduce((x0, [_, sigs]) => x0 + sigs.length, 0) +
  //         (signatures.length - 1) * 0.5,
  //       y: 1.04,
  //       text: `<b>${mutation}</b>`,
  //       showarrow: false,
  //       font: {
  //         size: 18,
  //       },
  //       align: "center",
  //     })
  //   );

  //   const xannotations = flatSorted.map((num, index) => ({
  //     xref: "x",
  //     yref: "paper",
  //     xanchor: "bottom",
  //     yanchor: "bottom",
  //     x: index,
  //     y: -0.1,
  //     text: num.mutationType.replace(
  //       /\[(.*)\]/,
  //       num.mutationType.substring(2, 3)
  //     ),
  //     showarrow: false,
  //     font: {
  //       size: 10,
  //       color: colors[num.mutationType.substring(2, 5)],
  //     },
  //     align: "center",
  //     num: num,
  //     index: index,
  //     textangle: -90,
  //   }));

  //   const sampleAnnotation = {
  //     xref: "paper",
  //     yref: "paper",
  //     xanchor: "bottom",
  //     yanchor: "bottom",
  //     x: 0,
  //     y: 0.9,
  //     text:
  //       "<b>" + sample + ": " + numberWithCommas(totalMutations) + " subs </b>",
  //     showarrow: false,
  //     font: {
  //       size: 18,
  //     },
  //     align: "center",
  //   };

  //   const shapes = Object.entries(groupByMutation).map(
  //     ([mutation, _], groupIndex, array) => ({
  //       type: "rect",
  //       xref: "x",
  //       yref: "paper",
  //       x0: array
  //         .slice(0, groupIndex)
  //         .reduce((x0, [_, sigs]) => x0 + sigs.length, -0.4),
  //       x1: array
  //         .slice(0, groupIndex + 1)
  //         .reduce((x0, [_, sigs]) => x0 + sigs.length, -0.6),
  //       y0: 1.05,
  //       y1: 1.01,
  //       fillcolor: colors[mutation],
  //       line: {
  //         width: 0,
  //       },
  //       mutation: mutation,
  //     })
  //   );
  //   // console.log("shapes");
  //   // console.log(shapes);

  const layout = {
    hoverlabel: { bgcolor: "#FFF" },
    xaxis: {
      showticklabels: true,
      showline: true,
      tickfont: {
        size: 10,
      },
    },
  };
  console.log("layout");
  console.log(layout);

  return { traces, layout };
}
