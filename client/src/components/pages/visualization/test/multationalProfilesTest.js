import React, { useEffect, useState } from "react";
import barData from "./data.json";
import Plot from "react-plotly.js";
import { chordTranspose } from "d3";

export default function MultationalProfilesTest() {
  console.log(barData);
  const regex = /(\[.*\])/;
  const regex2 = /[\[\]']+/g;
  function average(arr) {
    const sum = arr.reduce((a, b) => a + b, 0);
    return sum / arr.length || 0;
  }

  const groupByMutation = barData.reduce((groups, v) => {
    const mutation = v.mutationType.match(regex)[0];
    const signature = {
      mutationType: v.mutationType,
      contribution: v.contributions[0],
    };
    groups[mutation] = groups[mutation]
      ? [...groups[mutation], signature]
      : [signature];
    return groups;
  }, {});

  console.log(groupByMutation);

  const colors = {
    "[C>A]": "#03BCEE",
    "[C>G]": "black",
    "[C>T]": "#E32926",
    "[T>A]": "#CAC9C9",
    "[T>C]": "#A1CE63",
    "[T>G]": "#EBC6C4",
  };

  const data1 = Object.entries(groupByMutation).map(
    ([mutation, signatures], groupIndex, array) => ({
      name: mutation.replace(regex2, ""),
      type: "scatter",
      marker: { color: colors[mutation], symbol: "circle-open-dot" },
      mode: "markers",
      x: signatures.map(
        (e, i) =>
          array
            .slice(0, groupIndex)
            .reduce((x0, [_, sigs]) => x0 + sigs.length, 0) +
          i +
          0.5
      ),
      y: signatures.map((e) => e.contribution),
      test0: signatures.map((e, i) => array.slice(0, groupIndex)),
      testx0: signatures.map((e, i) =>
        array.slice(0, groupIndex).reduce((x0, curr) => x0, 0)
      ),
      testcurr: signatures.map((e, i) =>
        array.slice(0, groupIndex).reduce((x0, curr) => curr, 0)
      ),
      test1: signatures.map((e, i) =>
        array.slice(0, groupIndex).reduce((x0, curr) => curr.length, 0)
      ),
      test2: signatures.map((e, i) =>
        array.slice(0, groupIndex).reduce((x0, [_, curr]) => curr.length, 0)
      ),
      hoverinfo: "x+y",
      showlegend: false,
      array: array,
    })
  );

  console.log(data1);

  // const data2 = Object.entries(groupByMutation).map(
  //   ([mutation, signatures]) => ({
  //     name: mutation,
  //     type: "bar",
  //     marker: { color: colors[mutation] },
  //     x: signatures.map((e) => e.mutationType),
  //     y: [0.01],
  //     xaxis: "x2",
  //     yaxis: "y2",
  //     text: mutation,
  //     showlegend: false,
  //     hoverinfo: "none",
  //   })
  // );

  const shapes1 = Object.entries(groupByMutation).map(
    ([mutation, signatures], groupIndex, array) => ({
      type: "rect",
      xref: "x",
      yref: "paper",
      //x0: signatures.map((e) => e.mutationType)[0],
      //x1: signatures.map((e) => e.mutationType)[signatures.length - 1],
      x0: array
        .slice(0, groupIndex)
        .reduce((x0, [_, sigs]) => x0 + sigs.length, 0),
      x1: array
        .slice(0, groupIndex + 1)
        .reduce((x0, [_, sigs]) => x0 + sigs.length, 0),
      y0: 1.03,

      y1: 1,
      fillcolor: colors[mutation],
      line: {
        width: 0,
      },
    })
  );

  console.log(shapes1);

  const shapes2 = Object.entries(groupByMutation).map(
    ([mutation, signatures], groupIndex, array) => ({
      type: "rect",
      xref: "x",
      yref: "paper",
      //x0: signatures.map((e) => e.mutationType)[0],
      y0: 1,
      //x1: signatures.map((e) => e.mutationType)[signatures.length - 1],
      x0: array
        .slice(0, groupIndex)
        .reduce((x0, [_, sigs]) => x0 + sigs.length, 0),
      x1: array
        .slice(0, groupIndex + 1)
        .reduce((x0, [_, sigs]) => x0 + sigs.length, 0),
      y1: 0,
      fillcolor: colors[mutation],
      line: {
        width: 0,
      },
      opacity: 0.2,
    })
  );

  console.log(shapes2);

  const shapes = [...shapes1, ...shapes2];
  const annotations = Object.entries(groupByMutation).map(
    ([mutation, signatures], groupIndex, array) => ({
      xref: "x",
      yref: "paper",
      x:
        array
          .slice(0, groupIndex)
          .reduce((x0, [_, sigs]) => x0 + sigs.length, 0) +
        (signatures.length - 1) * 0.5,
      xanchor: "bottom",
      y: 1.05,
      yanchor: "bottom",
      text: "<b>" + mutation.replace(regex2, "") + "</b>",
      showarrow: false,
      font: {
        color: colors[mutation],
        size: 18,
      },
      align: "center",
    })
  );

  const data = [...data1];
  console.log(data);
  console.log(shapes);
  console.log(annotations);

  var layout = {
    xaxis: {
      title: "Substitution",
      showline: true,
      showticklabels: true,
      tickangle: -90,
      tickfont: {
        size: 10,
      },
    },
    yaxis: {
      title: "Mutation probability",
      //tickformat: ".1%",
      autorange: true,
    },

    shapes: shapes,
    annotations: annotations,
  };
  console.log(layout);

  return (
    <Plot
      data={data}
      layout={layout}
      useResizeHandler={true}
      style={{ width: "100%", height: "100%" }}
    />
  );
}
