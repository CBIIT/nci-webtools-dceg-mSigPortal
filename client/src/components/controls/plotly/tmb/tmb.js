import { groupBy } from "lodash";

export default function TMB(data, study = "PCAWG") {
  // Calculate the number of mutations per megabase for each study
  const genome = { PCAWG: "GRCh37", TCGA: "GRCh37" };
  const genomeSize = { GRCh37: 3101976562 / Math.pow(10, 6) };
  const burden = (exposure) => Math.log10(exposure / genomeSize[genome[study]]);

  function average(arr) {
    const sum = arr.reduce((a, b) => a + b, 0);
    return sum / arr.length || 0;
  }
  const groupByCancer = groupBy(data, "cancer");

  // calculate burden per cancer/sample
  const cancerBurden = Object.entries(groupByCancer)
    .map(([cancer, values]) => {
      const groupBySample = groupBy(values, 'sample');
      // sum exposure values per sample group and calculate burden
      const tmbs = Object.entries(groupBySample).map(([_, e]) =>
        burden(e.reduce((sum, e) => e.exposure + sum, 0))
      );
      return { cancer, tmbs };
    })
    .sort((a, b) => (a.cancer > b.cancer ? 1 : b.cancer > a.cancer ? -1 : 0));

  console.log(cancerBurden);

  const traces = cancerBurden.map((element, index, array) => ({
    element: element,
    index: index,
    array: array,
    name: `${element.cancer}`,
    type: "scatter",
    marker: { symbol: "circle-open", size: 4, color: "black" },
    mode: "markers",
    y: element.tmbs.map((e) => e),
    average: average(element.tmbs.map((e) => e)),
    hovertemplate: "Burden: %{y}<br>",
    // x: element.tmbs.map(
    //   (e, i) =>
    //     array
    //       .slice(0, index)
    //       .reduce((x0, curr) => x0 + curr.tmbs.length, 0) +
    //     i +
    //     0.5
    // ),
    x:
      array.length > 1
        ? element.tmbs.map(
            (e, i) => index + 0.7 + (0.3 / element.tmbs.length) * i
          )
        : element.tmbs.map(
            (e, i) => index + 0.07 + (0.8 / element.tmbs.length) * i
          ),
    showlegend: true,
  }));

  // var trace1 = {
  //   x: [1, 2, 3, 4, 5],
  //   y: [1, 6, 3, 6, 1],
  //   mode: "markers",
  //   type: "scatter",
  //   name: "Team A",
  //   text: ["A-1", "A-2", "A-3", "A-4", "A-5"],
  //   marker: { size: 12 },
  // };

  // var trace2 = {
  //   x: [1.5, 2.5, 3.5, 4.5, 5.5],
  //   y: [4, 1, 7, 1, 4],
  //   mode: "markers",
  //   type: "scatter",
  //   name: "Team B",
  //   text: ["B-a", "B-b", "B-c", "B-d", "B-e"],
  //   marker: { size: 12 },
  // };
  // const traces = [trace1, trace2];
  console.log("traces:--");
  console.log(traces);

  const topLabel = cancerBurden.map((element, index, array) => ({
    xref: "x",
    yref: "paper",
    xanchor: "bottom",
    yanchor: "bottom",
    // x:
    //   array
    //     .slice(0, index)
    //     .reduce((x0, curr) => x0 + curr.tmbs.length, 0) +
    //   (element.tmbs.length - 1) * 0.5,
    //x: array.length * 0.5,
    x: array.length > 1 ? index : (index + index + 1) * 0.5,
    y: 1.0,
    text: `${element.cancer}`,
    showarrow: false,
    // font: {
    //   size: 12,
    // },
    align: "right",
    textangle: 45,
  }));
  console.log("top label:--");
  console.log(topLabel);

  const bottoLabel1 = cancerBurden.map((element, index, array) => ({
    xref: "x",
    yref: "paper",
    xanchor: "bottom",
    yanchor: "bottom",
    x: (index + index + 1) * 0.5,
    y: -0.16,
    text: `${element.tmbs.length}<br> - <br>`,
    showarrow: false,
    font: {
      size: 12,
      color: "blue",
    },
    align: "center",
  }));

  const bottoLabel2 = cancerBurden.map((element, index, array) => ({
    xref: "x",
    yref: "paper",
    xanchor: "bottom",
    yanchor: "bottom",
    x: (index + index + 1) * 0.5,
    y: -0.21,
    text: `${element.tmbs.length}`,
    showarrow: false,
    font: {
      size: 12,
      color: "green",
    },
    align: "center",
  }));

  // const shapes = cancerBurden.map((element, index, array) => ({
  //   type: "rect",
  //   xref: "x",
  //   yref: "paper",

  //   // x0: array
  //   //   .slice(0, index)
  //   //   .reduce((x0, curr) => x0 + curr.tmbs.length, 0),

  //   // x1: array
  //   //   .slice(0, index + 1)
  //   //   .reduce((x0, curr) => x0 + curr.tmbs.length, 0),
  //   x0: index,
  //   x1: index + 1,
  //   y0: 0,
  //   y1: 1,
  //   fillcolor: index % 2 === 0 ? "gray" : "#F8F8F8",
  //   line: {
  //     width: 0,
  //   },
  //   opacity: 0.2,
  // }));
  // console.log("shapes:--");
  // console.log(shapes);

  const lines = cancerBurden.map((element, index, array) => ({
    type: "line",
    xref: "x",
    yref: "y",

    // x0: array
    //   .slice(0, index)
    //   .reduce((x0, curr) => x0 + curr.tmbs.length, 0),

    // x1: array
    //   .slice(0, index + 1)
    //   .reduce((x0, curr) => x0 + curr.tmbs.length, 0),
    x0: index + 0.1,
    x1: index + 0.9,

    y0: average(element.tmbs.map((e) => e)),
    y1: average(element.tmbs.map((e) => e)),
    line: {
      width: 1,
      color: "red",
    },
  }));

  const layout = {
    // title: {
    //   text: "Tumor Mutational Burden Separated by Signatures",
    //   yanchor: "top",
    // },
    xaxis: {
      showticklabels: true,
      tickfont: {
        size: 10,
      },
      linecolor: "black",
      linewidth: 2,
      mirror: true,
      tickmode: "array",
      showgrid: false,
      //tickvals: flatSorted.map((_, i) => i),
      //ticktext: flatSorted.map((_, i) => i),
    },
    yaxis: {
      title: "Number of Mutations per Megabase<br>(log10)",
      autorange: true,
      zeroline: false,
      //showline: true,
      linecolor: "black",
      linewidth: 2,
      mirror: true,
      automargin: true,
    },

    //shapes: [...shapes, ...lines],
    //annotations: [...topLabel, ...bottoLabel1, ...bottoLabel2],
  };

  console.log("layout:");
  console.log(layout);

  var config = {
    //responsive: true,
  };

  return { traces: [...traces], layout: layout, config };
  //return { traces, layout };
}
