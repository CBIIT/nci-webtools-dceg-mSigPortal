import React, { useEffect, useState } from 'react';
import barData from './data.json';
//import Plot from "../../../controls/plot/plot";
import Plot from 'react-plotly.js';

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
    'C>A': 'blue',
    'C>G': 'black',
    'C>T': 'red',
    'T>A': 'grey',
    'T>C': 'green',
    'T>G': 'pink',
  };

  const data = Object.entries(groupByMutation).map(
    ([mutation, signatures]) => ({
      name: mutation,
      type: 'bar',
      marker: { color: colors[mutation] },
      x: signatures.map((e) => e.mutationType),
      y: signatures.map((e) => e.contribution),
    })
  );

  // Object.entries(groupByMutation).forEach(([key, value]) => {
  //   valueTopX.push(key);
  //   console.log(value.length);

  //   value.forEach((element, index) => {
  //     valuesX.push(element.mutationType);
  //     valuesY.push(element.contribution);
  //     console.log(element);
  //   });
  // });

  const groups = Object.entries(groupByMutation);
  for (let i = 0; i < groups.length; i++) {
    valueTopX.push(groups[i][0]);
    traces[i] = [];
    valuesX[i] = [];
    valuesY[i] = [];
    for (let k = 0; k < groups[i][1].length; k++) {
      valuesX[i].push(groups[i][1][k].mutationType);
      valuesY[i].push(groups[i][1][k].contribution);
      traces[i].push({
        x: valuesX[i],
        y: valuesY[i],
        type: 'bar',
        name: groups[i][0],
      });
    }
  }

  console.log(groups);
  console.log(traces);

  var xValue = ['A'];
  var xValue2 = ['F', 'G', 'H', 'K', 'L'];
  var xValue3 = ['M', 'N', 'O', 'P', 'J'];

  var yValue = [20, 40, 10, 15, 15];
  var yValue2 = [15, 5, 25, 15, 80];

  var trace1 = {
    x: xValue,
    y: yValue,

    type: 'bar',
    name: 'Blue',
    text: yValue.map(String),
    textposition: 'auto',
    hoverinfo: 'none',
    marker: {
      color: 'blue',
    },
  };

  var trace2 = {
    x: xValue2,
    y: yValue2,

    type: 'bar',
    name: 'Orange',
    text: yValue2.map(String),
    textposition: 'auto',
    hoverinfo: 'none',
    marker: {
      color: 'orange',
    },
  };

  var trace3 = {
    x: xValue3,
    y: yValue2,

    type: 'bar',
    name: 'Green',
    text: yValue2.map(String),
    textposition: 'auto',
    hoverinfo: 'none',
    marker: {
      color: 'green',
    },
  };

  var trace4 = {
    x: xValue,
    y: [10],
    text: 'A>E',
    textposition: 'auto',
    xaxis: 'x2',
    yaxis: 'y2',
    type: 'bar',
    marker: {
      color: 'blue',
    },
    showlegend: false,
  };
  var trace5 = {
    x: xValue2,
    y: [10],
    text: 'F>L',
    textposition: 'auto',
    xaxis: 'x2',
    yaxis: 'y2',
    type: 'bar',
    marker: {
      color: 'orange',
    },
    showlegend: false,
  };

  var trace6 = {
    x: xValue3,
    y: [10],
    text: 'M>J',
    textposition: 'auto',
    xaxis: 'x2',
    yaxis: 'y2',
    type: 'bar',
    marker: {
      color: 'green',
    },
    showlegend: false,
  };

  // var data = [
  //   traces[0][0],
  //   traces[1][0],
  //   traces[2][0],
  //   traces[3][0],
  //   traces[4][0],
  //   traces[5][0],
  // ];

  var layout = {
    grid: {
      rows: 2,
      columns: 1,
      roworder: 'bottom to top',
      ygap: 0.0001,
    },
    xaxis: { showline: true, tickangle: -90 },
    yaxis: {
      //  title: 'Percentage of Single Base Substitution',
      tickformat: '.1%',
    },
    yaxis2: { visible: false, scaleanchor: 'y' },
    xaxis2: { visible: false },
  };

  return (
    <Plot
      data={data}
      layout={layout}
      useResizeHandler={true}
      style={{ width: '100%', height: '100%' }}
    />
  );
}
