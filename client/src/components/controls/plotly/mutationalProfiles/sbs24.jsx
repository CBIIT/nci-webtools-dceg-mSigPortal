export default function SBS24(data, title = '') {
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
    //hoverinfo: 'x+y',
    hovertemplate: '<b>Transcribed Strand</b><br>%{y}, %{x} <extra></extra>',
    orientation: 'h',
  };

  const untranscribedTraces = {
    name: 'Untranscribed Strand',
    type: 'bar',
    marker: { color: '#E32925' },
    x: untranscribed.map((e) => e.mutations),
    y: untranscribed.map((e) => e.mutationType.slice(-3)),
    //hoverinfo: 'x+y',
    hovertemplate: '<b>Untranscribed Strand</b><br>%{y}, %{x} <extra></extra>',
    orientation: 'h',
  };

  const traces = [untranscribedTraces, transcribedTraces];

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
        '<b>' + data[0].sample ||
        data[0].signatureName +
          ': ' +
          totalMutations.toLocaleString(undefined) +
          ' transcribed subs</b>',
      font: {
        size: 26,
        family: 'Arial',
      },
      xref: 'paper',
      x: 0.01,
      y: 0.9,
    },
    xaxis: {
      title: {
        text: '<b>Number of Single Base Substitutions</b>',
        font: {
          size: 20,
          family: 'Times New Roman',
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
      tickformat: maxMutation > 1000 ? '~s' : '',
    },
    yaxis: {
      tickfont: {
        size: 16,
      },
      linecolor: '#E0E0E0',
      linewidth: 1,
      tickformat: maxMutation > 1000 ? '~s' : '',
      categoryorder: 'category descending',
    },
  };

  return { traces, layout };
}
