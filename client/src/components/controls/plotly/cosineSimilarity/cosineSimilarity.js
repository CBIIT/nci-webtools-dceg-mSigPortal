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
        y: 1.01,
        len: 450,
        lenmode: 'pixels',
      },
    },
  ];

  const layout = {
    autosize: true,
    width: signatureOrder.length < 10 ? 400 + 50 * signatureOrder.length : null,
    height: 300 + 50 * sampleOrder.length,
    xaxis: {
      automargin: true,
      showticklabels: signatureOrder.length > 100 ? false : true,
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
