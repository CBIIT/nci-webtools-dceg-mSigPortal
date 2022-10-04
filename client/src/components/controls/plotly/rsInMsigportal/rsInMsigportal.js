import { groupBy } from 'lodash';
export default function RsInMsigportal(rawData) {
  console.log(rawData);

  //   const grouped = groupBy(
  //     rawData,
  //     (item) => `"${item.species}_${item.profile}_${item.matrix}"`
  //   );
  //   console.log(grouped);

  const data = rawData.reduce((acc, e, i) => {
    const human = '(GRCh37/38)';
    const species = e.species.includes(human);

    acc[species] = acc[species] ? [...acc[species], e] : [e];
    return acc;
  }, {});

  const humanData = data.true;
  console.log(humanData);
  const grouped = groupBy(
    humanData,
    (item) => `"${item.profile}${item.matrix}"`
  );
  console.log(grouped);

  const tracesPie0 = {
    type: 'pie',
    // labels: grouped[0].map((group) => group.signatureName),
    // values: grouped[0].map((group) =>
    //   group.samples.reduce((a, b) => a + b.exposure, 0)
    // ),

    textposition: 'inside',
  };

  const tracePies = Object.entries(grouped).map(
    ([key, element], index, array) => ({
      type: 'pie',
      e: element,

      labels: element.map((e) => e.signatureSetName),
      values: element.map((e) => parseInt(e.count)),
      name: key,
    })
  );
  console.log(tracePies);
  const traces = [tracePies[0]];
  console.log(traces);
  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: 500,
    autosize: true,
    grid: { rows: 2, columns: 4 },
  };
  var config = {
    //responsive: true,
  };
  return { traces, layout, config };
}
