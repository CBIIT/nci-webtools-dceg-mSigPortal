import { groupBy } from 'lodash';
export default function MsIndividual(data, arg) {
  console.log('MS Individual Plot');
  console.log(data);
  console.log(arg);
  const exposureData = data[0].data;
  console.log(exposureData);
  const signatureData = data[1].data;
  console.log(signatureData);
  const segmatrixData = data[2].data;
  console.log(segmatrixData);

  const exposure_groupBySignature = groupBy(
    exposureData.filter((o) => o['exposure'] > 0),
    'signatureName'
  );
  console.log(exposure_groupBySignature);

  const signatureNames = Object.keys(exposure_groupBySignature).map((e) => e);
  console.log(signatureNames);
  const exposureSum = Object.values(exposure_groupBySignature).reduce(
    (n, { exposure }) => n + exposure,
    0
  );
  console.log(exposureSum);
  const signature_groupBySignature = groupBy(signatureData, 'signatureName');
  console.log(signature_groupBySignature);

  const seqmatrix_groupByMutationType = groupBy(segmatrixData, 'mutationType');
  console.log(seqmatrix_groupByMutationType);

  const colors = {
    'C>A': '#03BCEE',
    'C>G': 'black',
    'C>T': '#E32926',
    'T>A': '#CAC9C9',
    'T>C': '#A1CE63',
    'T>G': '#EBC6C4',
  };

  const mutationRegex = /\[(.*)\]/;
  const mutationLabels = (e) => `<b>${e.mutation}</b>`;
  const formatTickLabels = (mutationGroups) =>
    mutationGroups
      .map(({ mutation, data }) =>
        data.map((e) => {
          const color = colors[mutation];
          const regex = /^(.)\[(.).{2}\](.)$/;
          const match = e.mutationType.match(regex);

          return `${match[1]}<span style="color:${color}"><b>${match[2]}</b></span>${match[3]}`;
        })
      )
      .flat();

  function groupDataByMutation(
    data,
    groupRegex,
    mutationGroupSort = false,
    mutationTypeSort = false
  ) {
    const groupByMutation = data.reduce((acc, e) => {
      const mutation = e.mutationType.match(groupRegex)[1];
      acc[mutation] = acc[mutation] ? [...acc[mutation], e] : [e];
      return acc;
    }, {});

    const groupedData = Object.entries(groupByMutation).map(
      ([mutation, data]) => ({
        mutation,
        data: mutationTypeSort ? data.sort(mutationTypeSort) : data,
      })
    );

    return mutationGroupSort
      ? groupedData.sort(mutationGroupSort)
      : groupedData;
  }

  function getTotalMutations(data) {
    return data.reduce(
      (total, e) => total + (e.mutations || e.contribution || 0),
      0
    );
  }

  function getMaxMutations(data) {
    return Math.max(...data.map((e) => e.mutations || e.contribution || 0));
  }

  function getRss(sampleDifferenceData) {
    const squareDiff = sampleDifferenceData.map((e) => Math.pow(e || 0, 2));
    return squareDiff.reduce((a, b, i) => a + b, 0).toExponential(3);
  }

  function getCosineSimilarity(data1, data2) {
    function dotp(x, y) {
      function dotp_sum(a, b) {
        return a + b;
      }
      function dotp_times(a, i) {
        return x[i] * y[i];
      }
      return x.map(dotp_times).reduce(dotp_sum, 0);
    }

    function cosineSimilarity(A, B) {
      var similarity =
        dotp(A, B) / (Math.sqrt(dotp(A, A)) * Math.sqrt(dotp(B, B)));
      return similarity;
    }
    return cosineSimilarity(
      data1.map((e) => e || 0),
      data2.map((e) => e || 0)
    ).toFixed(3);
  }

  var trace1 = {};

  var trace2 = {};

  var traces = [trace1, trace2];
  const layout = {
    showlegend: true,
    hoverlabel: { bgcolor: '#FFF' },
    height: 900,
  };
  return { traces: traces, layout: layout };
}
