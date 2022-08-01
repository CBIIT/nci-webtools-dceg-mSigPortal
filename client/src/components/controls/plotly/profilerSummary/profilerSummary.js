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

  const layout = {
    autosize: true,
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
