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
  const dataArray = [];
  Object.values(grouped).map((e) => dataArray.push(e));
  console.log(dataArray);

  const tracePies = Object.entries(grouped).map(
    ([key, element], index, array) => ({
      type: 'pie',
      e: element,
      textposition: 'inside',
      labels: element.map((e) => e.signatureSetName),
      values: element.map((e) => parseInt(e.count)),
      name: key,
      domain: {
        row: index <= 4 ? 0 : 1,
        column: index <= 4 ? index : index - 5,
      },
    })
  );
  console.log(tracePies);
  const traces = [...tracePies];
  //const traces = [pie0, pie1];
  console.log(traces);
  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: 500,
    autosize: true,
    grid: { rows: 2, columns: 5 },
  };
  const config = {
    //responsive: true,
  };
  return { traces: traces, layout: layout, config };
}
