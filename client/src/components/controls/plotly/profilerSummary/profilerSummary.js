import { groupBy } from 'lodash';

export default function profilerSummary(data) {
  const groupBySampleProfile = groupBy(data, (e) => `${e.sample}_${e.profile}`);
  const summaraizeMutations = Object.entries(groupBySampleProfile)
    .map(([sampleProfile, samples]) => {
      const [sample, profile] = sampleProfile.split('_');
      const matrix = [...new Set(samples.map((e) => e.matrix))].sort(
        (a, b) => a - b
      );
      const mutations = samples.reduce((acc, e) => acc + e.mutations, 0);
      return {
        sample,
        profile,
        matrix,
        mutations: Math.log10(mutations),
      };
    })
    .sort((a, b) => a.mutations - b.mutations);

  const groupByProfile = groupBy(summaraizeMutations, 'profile');

  const orderByMean = Object.entries(groupByProfile)
    .map(([profile, samples]) => ({
      name: `${profile}: ${samples[0].matrix.join('/')}`,
      samples,
      mean: samples.reduce((acc, e) => acc + e.mutations, 0) / samples.length,
    }))
    .sort((a, b) => b.mean - a.mean);

  const traces = orderByMean.map((e, i, array) => {
    return {
      name: e.name,
      x: array[0].samples.map((s) => s.sample),
      y: e.samples.map((s) => s.mutations),
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
      title: 'Samples',
      showline: true,
      mirror: true,
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
