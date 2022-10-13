import { groupBy } from 'lodash';

export default function MsDecomposition(data, arg) {
  console.log(data);
  console.log(arg);
  const groupByCancer = groupBy(data, 'cancer');
  console.log(groupByCancer);

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

  function calculate_simularities(
    original_genomes,
    signature,
    signature_activities
  ) {}
  const traces = {};

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    bargap: 0.3,
    height: 450,
    // width: 1080,
  };

  return { traces, layout };
}
