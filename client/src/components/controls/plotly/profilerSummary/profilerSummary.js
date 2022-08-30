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

  // sort samples of other profiles to match the sample order of the profile with the largest mean mutation value
  const topSamples = data[0].samples.map((s) => s.sample);
  const allSamples = [...new Set(rawData.map((e) => e.sample))].sort(
    (a, b) => topSamples.indexOf(a) - topSamples.indexOf(b)
  );
  const sortedData = data.map((e) => ({
    ...e,
    samples: e.samples.sort(
      (a, b) => allSamples.indexOf(a.sample) - allSamples.indexOf(b.sample)
    ),
  }));

  const traces = sortedData.map((e, i, array) => {
    return {
      name: e.name,
      x: e.samples.map((s) => s.sample),
      y: e.samples.map((s) => s.logTotalMutations),
      mode: 'lines+markers',
      type: 'scatter',
    };
  });

  // find the longest label to calculate extra height margin
  const labels = traces[0].x;
  const longest = labels.reduce((a, e) => (a > e.length ? a : e.length), 0);
  const extraMargin = longest < 10 ? 60 : longest * 7;

  const layout = {
    autosize: true,
    height: 500 + extraMargin,
    margin: { t: 40, b: extraMargin },
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
      type: 'category',
      categoryorder: 'array',
      categoryarray: allSamples,
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
