export default function profilerSummary(data) {
  const traces = data.map((e, i, array) => {
    return {
      name: e.name,
      x: array[0].samples.map((s) => s.sample),
      y: e.samples.map((s) => s.logTotalMutations),
      mode: 'lines+markers',
      type: 'scatter',
    };
  });

  const bottomMargin = (labels) => {
    const longest = labels.reduce((a, e) => (a > e.length ? a : e.length), 0);
    if (longest < 10) return 60;
    else return longest * 6;
  };

  const layout = {
    autosize: true,
    margin: { t: 20, b: bottomMargin(traces[0].x) },
    legend: {
      title: {
        text: '<b>Profile</b>',
      },
    },
    xaxis: {
      showline: true,
      mirror: true,
      tickangle: 45,
      range: [-1, traces[0].x.length],
    },
    yaxis: {
      title: 'log<sub>10</sub>(Mutations)',
      ticks: 'outside',
      zeroline: false,
      showline: true,
      mirror: true,
    },
  };

  const config = {
    responsive: true,
  };
  return { traces, layout, config };
}
