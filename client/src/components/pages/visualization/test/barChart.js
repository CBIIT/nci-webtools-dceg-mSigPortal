import React, { useEffect, useState } from "react";
import barData from "./data.json";
import Plot from "../../../controls/plot/plot";
//import Plot from "react-plotly.js";

export default function BarChart() {
  return (
    <Plot
      data={[
        {
          x: [1, 2, 3],
          y: [2, 6, 3],
          type: "scatter",
          mode: "lines+markers",
          marker: { color: "red" },
        },
        { type: "bar", x: [1, 2, 3], y: [2, 5, 3] },
      ]}
      layout={{ width: 320, height: 240, title: "A Fancy Plot" }}
    />
  );
}
