export default function SBS24(data, sample) {
  // filter transcribed and untranscribed data
  const transcribed = data.filter((e) => /^T:/.test(e.mutationType));
  const untranscribed = data.filter((e) => /^U:/.test(e.mutationType));

  const totalMutations =
    transcribed.reduce((total, e) => total + e.mutations, 0) +
    untranscribed.reduce((total, e) => total + e.mutations, 0);
  const maxMutation = Math.max(
    ...[
      ...transcribed.map((e) => e.mutations),
      ...untranscribed.map((e) => e.mutations),
    ]
  );

  const transcribedTraces = {
    name: 'Transcribed Strand',
    type: 'bar',
    marker: { color: '#004765' },
    x: transcribed.map((e) => e.mutations),
    y: transcribed.map((e) => e.mutationType.slice(-3)),
    hoverinfo: 'x+y',
    orientation: 'h',
  };

  const untranscribedTraces = {
    name: 'Untranscribed Strand',
    type: 'bar',
    marker: { color: '#E32925' },
    x: untranscribed.map((e) => e.mutations),
    y: untranscribed.map((e) => e.mutationType.slice(-3)),
    hoverinfo: 'x+y',
    orientation: 'h',
  };

  const traces = [transcribedTraces, untranscribedTraces];

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    bargap: 0.3,
    height: 600,
    width: 750,
    legend: {
      x: 1,
      xanchor: 'right',
      y: 1,
      traceorder: 'reversed',
      bgcolor: '#FFFFFF',
      bordercolor: '#E1E1E1',
      borderwidth: 1,
    },
    title: {
      text:
        '<b>' +
        sample +
        ': ' +
        totalMutations.toLocaleString(undefined) +
        ' transcribed subs </b>',
      font: {
        size: 24,
      },
      xref: 'paper',
      x: 0.05,
    },
    xaxis: {
      title: {
        text: '<b>Number of Single Base Substitutions</b>',
        font: {
          size: 18,
        },
      },
      tickfont: {
        size: 16,
      },
      ticks: 'outside',
      linecolor: '#E0E0E0',
      linewidth: 1,
      showgrid: false,
      autorange: false,
      range: [0, maxMutation * 1.2],
      tickformat: '~s',
    },
    yaxis: {
      tickfont: {
        size: 16,
      },
      linecolor: '#E0E0E0',
      linewidth: 1,
      tickformat: '~s',
      categoryorder: 'category descending',
    },
  };

  return { traces, layout };
}
