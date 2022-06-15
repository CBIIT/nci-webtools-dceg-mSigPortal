export default function TMB(data, study) {
  const genome = { PCAWG: "GRCh37", TCGA: "GRCh37" };
  const genomeSize = { GRCh37: 3101976562 / Math.pow(10, 6) };

  function average(arr) {
    const sum = arr.reduce((a, b) => a + b, 0);
    return sum / arr.length || 0;
  }
  const groupByCancer = data.reduce((groups, e) => {
    const cancerName = e.Cancer_Type;
    groups[cancerName] = groups[cancerName] ? [...groups[cancerName], e] : [e];
    return groups;
  }, {});

  const groupBySamples = Object.entries(groupByCancer).map(
    ([cancerName, cancerArray]) => {
      const sampleGroups = cancerArray.reduce((sampleGroup, e) => {
        const sampleName = e.Sample;
        sampleGroup[sampleName] = sampleGroup[sampleName]
          ? [...sampleGroup[sampleName], e]
          : [e];
        return sampleGroup;
      }, {});

      return { cancer: cancerName, samples: sampleGroups };
    },
    {}
  );
  console.log("data:");
  console.log(data);
  console.log("group by cancer:");
  console.log(groupByCancer);

  const burdenQuotient = genomeSize[genome[study]];
  const cancerBurden = groupBySamples.map(({ cancer, samples }) => {
    const sampleBurden = Object.entries(samples).map(
      ([sampleName, sampleArray]) => ({
        sample: sampleName,
        burden: Math.log10(
          sampleArray.reduce((sum, e) => e.Exposure + sum, 0) / burdenQuotient
        ),
      })
    );
    sampleBurden.sort((a, b) => a.burden - b.burden);

    return { cancer, sampleBurden };
  });

  const flatSorted = Object.values(cancerBurden).flat();

  console.log("Cancer Budern:");
  console.log(cancerBurden);

  console.log("flatSorted:");
  console.log(flatSorted);

  const traces = cancerBurden.map((element, index, array) => ({
    element: element,
    index: index,
    array: array,
    name: `${element.cancer}`,
    type: "scatter",
    marker: { symbol: "circle-open", size: 4, color: "black" },
    mode: "markers",
    y: element.sampleBurden.map((e) => e.burden),
    average: average(element.sampleBurden.map((e) => e.burden)),
    hovertemplate: "Burden: %{y}<br>",
    // x: element.sampleBurden.map(
    //   (e, i) =>
    //     array
    //       .slice(0, index)
    //       .reduce((x0, curr) => x0 + curr.sampleBurden.length, 0) +
    //     i +
    //     0.5
    // ),
    x:
      array.lenght > 1
        ? element.sampleBurden.map(
            (e, i) => index + 0.7 + (0.3 / element.sampleBurden.length) * i
          )
        : element.sampleBurden.map(
            (e, i) => index + 0.07 + (0.8 / element.sampleBurden.length) * i
          ),
    showlegend: false,
  }));
  console.log("traces:--");
  console.log(traces);

  const topLabel = cancerBurden.map((element, index, array) => ({
    xref: "x",
    yref: "paper",
    xanchor: "bottom",
    yanchor: "bottom",
    // x:
    //   array
    //     .slice(0, index)
    //     .reduce((x0, curr) => x0 + curr.sampleBurden.length, 0) +
    //   (element.sampleBurden.length - 1) * 0.5,
    //x: array.length * 0.5,
    x: array.length > 1 ? index : (index + index + 1) * 0.5,
    y: 1.0,
    text: `${element.cancer}`,
    showarrow: false,
    // font: {
    //   size: 12,
    // },
    align: "right",
    textangle: 45,
  }));
  console.log("top label:--");
  console.log(topLabel);

  const bottoLabel1 = cancerBurden.map((element, index, array) => ({
    xref: "x",
    yref: "paper",
    xanchor: "bottom",
    yanchor: "bottom",
    x: (index + index + 1) * 0.5,
    y: -0.12,
    text: `${element.sampleBurden.length}<br> - <br>`,
    showarrow: false,
    font: {
      size: 12,
      color: "blue",
    },
    align: "center",
  }));

  const bottoLabel2 = cancerBurden.map((element, index, array) => ({
    xref: "x",
    yref: "paper",
    xanchor: "bottom",
    yanchor: "bottom",
    x: (index + index + 1) * 0.5,
    y: -0.17,
    text: `${element.sampleBurden.length}`,
    showarrow: false,
    font: {
      size: 12,
      color: "green",
    },
    align: "center",
  }));

  const shapes = cancerBurden.map((element, index, array) => ({
    type: "rect",
    xref: "x",
    yref: "paper",

    // x0: array
    //   .slice(0, index)
    //   .reduce((x0, curr) => x0 + curr.sampleBurden.length, 0),

    // x1: array
    //   .slice(0, index + 1)
    //   .reduce((x0, curr) => x0 + curr.sampleBurden.length, 0),
    x0: index,
    x1: index + 1,
    y0: 0,
    y1: 1,
    fillcolor: index % 2 === 0 ? "gray" : "#F8F8F8",
    line: {
      width: 0,
    },
    opacity: 0.2,
  }));

  const lines = cancerBurden.map((element, index, array) => ({
    type: "line",
    xref: "x",
    yref: "y",

    // x0: array
    //   .slice(0, index)
    //   .reduce((x0, curr) => x0 + curr.sampleBurden.length, 0),

    // x1: array
    //   .slice(0, index + 1)
    //   .reduce((x0, curr) => x0 + curr.sampleBurden.length, 0),
    x0: index + 0.1,
    x1: index + 0.9,

    y0: average(element.sampleBurden.map((e) => e.burden)),
    y1: average(element.sampleBurden.map((e) => e.burden)),
    line: {
      width: 1,
      color: "red",
    },
  }));

  console.log("shapes:--");
  console.log(shapes);

  const layout = {
    // title: {
    //   text: "Tumor Mutational Burden Separated by Signatures",
    //   yanchor: "top",
    // },
    xaxis: {
      showticklabels: false,
      tickfont: {
        size: 10,
      },
      linecolor: "black",
      linewidth: 2,
      mirror: true,
      tickmode: "array",
      showgrid: false,
      //tickvals: flatSorted.map((_, i) => i),
      //ticktext: flatSorted.map((_, i) => i),
    },
    yaxis: {
      title: "Number of Mutations per Megabase<br>(log10)",
      autorange: true,
      zeroline: false,
      //showline: true,
      linecolor: "black",
      linewidth: 2,
      mirror: true,
      automargin: true,
    },

    shapes: [...shapes, ...lines],
    annotations: [...topLabel, ...bottoLabel1, ...bottoLabel2],
  };
  console.log("layout:");
  console.log(layout);

  var config = {
    //responsive: true,
  };

  return { traces: [...traces], layout: layout, config };
  //return { traces, layout };
}
