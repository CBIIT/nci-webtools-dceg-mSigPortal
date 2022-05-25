import React, { useEffect, useState } from "react";
import barData from "./data.json";
import Plot from "react-plotly.js";

export default function BarChart() {
  console.log(barData);
  const traces = [];
  const valuesX = [];
  const valuesY = [];
  const valueTopX = [];

  const regex = /(\[.*\])/;
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
      name: mutation,
      type: "bar",
      marker: { color: colors[mutation] },
      x: signatures.map((e) => e.mutationType),
      y: signatures.map((e) => e.contribution),
      //showlegend: false,
    })
  );

  const data2 = Object.entries(groupByMutation).map(
    ([mutation, signatures]) => ({
      name: mutation,
      type: "bar",
      marker: { color: colors[mutation] },
      x: signatures.map((e) => e.mutationType),
      y: [0.01],
      xaxis: "x2",
      yaxis: "y2",
      text: mutation,
      showlegend: false,
      hoverinfo: "none",
    })
  );

  const shapes = Object.entries(groupByMutation).map(
    ([mutation, signatures]) => ({
      type: "rect",
      xref: "x",
      yref: "paper",
      x0: signatures.map((e) => e.mutationType)[0],
      y0: 1.03,
      x1: signatures.map((e) => e.mutationType)[15],
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
      x: signatures.map((e) => e.mutationType)[8],
      xanchor: "bottom",
      y: 1.05,
      yanchor: "bottom",
      text: mutation,
      showarrow: false,
      font: {
        color: colors[mutation],
        text: "bold",
      },
    })
  );

  const data = [...data1];
  console.log(data);
  console.log(data2);
  console.log(shapes);

  // const groups = Object.entries(groupByMutation);
  // for (let i = 0; i < groups.length; i++) {
  //   valueTopX.push(groups[i][0]);
  //   traces[i] = [];
  //   valuesX[i] = [];
  //   valuesY[i] = [];
  //   for (let k = 0; k < groups[i][1].length; k++) {
  //     valuesX[i].push(groups[i][1][k].mutationType);
  //     valuesY[i].push(groups[i][1][k].contribution);
  //     traces[i].push({
  //       x: valuesX[i],
  //       y: valuesY[i],
  //       type: "bar",
  //       name: groups[i][0],
  //     });
  //   }
  // }

  //console.log(groups);
  //console.log(traces);

  var layout = {
    grid: {
      //rows: 2,
      columns: 1,
      roworder: "bottom to top",
      //subplots: [["xy"], ["xy2"]],
      ygap: 0.1,
    },
    xaxis: {
      title: "Substitution",
      showline: true,
      tickangle: -90,
      automargin: true,
    },
    yaxis: {
      title: "Mutation probability",
      //tickformat: ".1%",
      autorange: true,
    },
    yaxis2: { visible: false, scaleanchor: "y", automargin: true },
    xaxis2: { visible: false },
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
