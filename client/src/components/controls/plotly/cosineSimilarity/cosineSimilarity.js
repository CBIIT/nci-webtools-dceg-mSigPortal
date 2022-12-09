export default function CosineSimilarity(inputData) {
  const { data, sampleOrder, signatureOrder } = inputData;

  const traces = [
    {
      x: signatureOrder,
      y: sampleOrder,
      z: sampleOrder.map((sample) =>
        signatureOrder
          .map(
            (signature) =>
              data.filter(
                (e) => e.signature == signature && e.sample == sample
              )[0].similarity
          )
          .flat()
      ),
      type: 'heatmap',
      colorscale: 'Viridis',
      colorbar: {
        title: 'Cosine Similarity',
        orientation: 'h',
        len: 0.5,
      },
    },
  ];

  const layout = {
    autosize: true,
    // width: 900,
    height: 1000,
    xaxis: {
      automargin: true,
    },
    yaxis: {
      automargin: true,
    },
  };

  const config = {
    responsive: true,
  };

  return { traces, layout, config };
}
