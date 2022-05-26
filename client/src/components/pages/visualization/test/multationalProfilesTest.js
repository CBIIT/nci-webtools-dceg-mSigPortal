import React, { useEffect, useState } from "react";
import barData from "./data.json";
import Plot from "react-plotly.js";

export default function BarChart() {
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
    "[C>A]": "blue",
    "[C>G]": "black",
    "[C>T]": "red",
    "[T>A]": "grey",
    "[T>C]": "green",
    "[T>G]": "pink",
  };

  const data1 = Object.entries(groupByMutation).map(
    ([mutation, signatures]) => ({
      name: mutation.replace(regex2, ""),
      type: "bar",
      marker: { color: colors[mutation] },
      x: signatures.map((e) => e.mutationType),
      y: signatures.map((e) => e.contribution),
      hoverinfo: "x+y",
      showlegend: false,
    })
  );

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

  const shapes = Object.entries(groupByMutation).map(
    ([mutation, signatures]) => ({
      type: "rect",
      xref: "x",
      yref: "paper",
      x0: signatures.map((e) => e.mutationType)[0],
      y0: 1.03,
      x1: signatures.map((e) => e.mutationType)[signatures.length - 1],
      y1: 1,
      fillcolor: colors[mutation],
      line: {
        width: 0,
      },
    })
  );

  const annotations = Object.entries(groupByMutation).map(
    ([mutation, signatures]) => ({
      xref: "x",
      yref: "paper",
      x: signatures.map((e) => e.mutationType)[
        Math.round(signatures.length / 2)
      ],
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

  return (
    <Plot
      data={data}
      layout={layout}
      useResizeHandler={true}
      style={{ width: "100%", height: "100%" }}
    />
  );
}
