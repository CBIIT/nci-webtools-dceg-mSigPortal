import { sbsColor } from '../../utils/colors';
export default function SBS6(data) {
  const colors = sbsColor;

  const totalMutations = data.reduce((total, e) => total + e.mutations, 0);
  const maxMutation = Math.max(...data.map((e) => e.mutations));

  const traces = data.map((e) => ({
    name: e.mutationType,
    type: 'bar',
    marker: { color: colors[e.mutationType] },
    y: [e.mutationType],
    x: [e.mutations],
    //hoverinfo: 'y+x',
    hovertemplate: '%{y}, %{x} <extra></extra>',
    orientation: 'h',
  }));

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    showlegend: false,
    height: 600,
    width: 750,
    title: {
      text:
        '<b>' + data[0].sample ||
        data[0].signatureName +
          ': ' +
          totalMutations.toLocaleString(undefined) +
          ' subs</b>',
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
          size: 22,
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
      range: [0, maxMutation * 1.25],
      tickformat: maxMutation > 1000 ? '~s' : '',
    },
    yaxis: {
      tickfont: {
        size: 16,
      },
      linecolor: '#E0E0E0',
      linewidth: 1,
      categoryorder: 'category descending',
    },
    tickformat: maxMutation > 1000 ? '~s' : '',
  };
  return { traces, layout };
}
