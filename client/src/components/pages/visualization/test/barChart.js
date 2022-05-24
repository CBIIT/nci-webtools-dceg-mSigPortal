import React, { useEffect, useState } from "react";
import barData from "./data.json";
//import Plot from "../../../controls/plot/plot";
import Plot from "react-plotly.js";

export default function BarChart() {
  console.log(barData);
  const traces = [];
  const valuesX = [];
  const valuesY = [];

  for (let i = 0; i < barData[0].contributions.length; i++) {
    valuesX[i] = [];
    valuesY[i] = [];
    for (let j = 0; j < barData.length; j++) {
      valuesX[i].push(barData[j].mutationType);
      valuesY[i].push(barData[j].contributions[i]);
    }
  }

  console.log(valuesX);
  console.log(valuesY);

  for (let i = 0; i < 6; i++) {
    traces.push({
      x: valuesX[i],
      y: valuesY[i],
      type: "bar",
      name: valuesX[i][i],
      text: valuesY[i].map(String),
    });
  }
  var xValue = ["A"];
  var xValue2 = ["F", "G", "H", "K", "L"];
  var xValue3 = ["M", "N", "O", "P", "J"];

  var yValue = [20, 40, 10, 15, 15];
  var yValue2 = [15, 5, 25, 15, 80];

  var trace1 = {
    x: xValue,
    y: yValue,

    type: "bar",
    name: "Blue",
    text: yValue.map(String),
    textposition: "auto",
    hoverinfo: "none",
    marker: {
      color: "blue",
    },
  };

  var trace2 = {
    x: xValue2,
    y: yValue2,

    type: "bar",
    name: "Orange",
    text: yValue2.map(String),
    textposition: "auto",
    hoverinfo: "none",
    marker: {
      color: "orange",
    },
  };

  var trace3 = {
    x: xValue3,
    y: yValue2,

    type: "bar",
    name: "Green",
    text: yValue2.map(String),
    textposition: "auto",
    hoverinfo: "none",
    marker: {
      color: "green",
    },
  };

  var trace4 = {
    x: xValue,
    y: [10],
    text: "A>E",
    textposition: "auto",
    xaxis: "x2",
    yaxis: "y2",
    type: "bar",
    marker: {
      color: "blue",
    },
    showlegend: false,
  };
  var trace5 = {
    x: xValue2,
    y: [10],
    text: "F>L",
    textposition: "auto",
    xaxis: "x2",
    yaxis: "y2",
    type: "bar",
    marker: {
      color: "orange",
    },
    showlegend: false,
  };

  var trace6 = {
    x: xValue3,
    y: [10],
    text: "M>J",
    textposition: "auto",
    xaxis: "x2",
    yaxis: "y2",
    type: "bar",
    marker: {
      color: "green",
    },
    showlegend: false,
  };

  var data = [traces[0], traces[1]];

  var layout = {
    grid: {
      rows: 2,
      columns: 1,
      roworder: "bottom to top",
      ygap: 0.0001,
    },
    yaxis2: { visible: false, scaleanchor: "y" },
    xaxis2: { visible: false },
  };

  return (
    <Plot
      data={traces}
      layout={layout}
      useResizeHandler={true}
      style={{ width: "100%", height: "100%" }}
    />
  );
}
