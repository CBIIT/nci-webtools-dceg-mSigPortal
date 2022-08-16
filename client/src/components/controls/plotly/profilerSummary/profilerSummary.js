import { groupBy } from 'lodash';

export default function profilerSummary(rawData) {
  const groupByProfileMatrix = groupBy(
    rawData,
    (e) => `${e.profile}_${e.matrix}`
  );
  const data = Object.values(groupByProfileMatrix)
    .map((samples) => {
      return {
        name: `${samples[0].profile}: ${samples[0].matrix}`,
        samples: samples.sort(
          (a, b) => a.logTotalMutations - b.logTotalMutations
        ),
        mean:
          samples.reduce(
            (acc, e) => acc + parseFloat(e.meanTotalMutations),
            0
          ) / samples.length,
      };
    })
    .sort((a, b) => b.mean - a.mean);

  // sort samples of other profiles to match the sample order of the profile with the largest mean
  const sampleOrder = data[0].samples.map((s) => s.sample);
  const sortedSamples = data.map((e) => ({
    ...e,
    samples: e.samples.sort(
      (a, b) => sampleOrder.indexOf(a.sample) - sampleOrder.indexOf(b.sample)
    ),
  }));

  const traces = sortedSamples.map((e, i, array) => {
    return {
      name: e.name,
      x: e.samples.map((s) => s.sample),
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
    margin: { t: 40, b: bottomMargin(traces[0].x) },
    legend: {
      title: {
        text: '<b>Profile</b>',
      },
      x: 1.01,
      y: 0.5,
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
