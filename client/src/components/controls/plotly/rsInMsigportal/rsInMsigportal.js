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
  const grouped = groupBy(humanData, (item) => `${item.profile}${item.matrix}`);
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
      texttemplate: '%{value}',
      direction: 'clockwise',
      name: key,
      domain: {
        // row: index <= 4 ? 0 : 1,
        // column: index <= 4 ? index : index - 5,
        x: [
          index < 4
            ? Math.round(index * (1 / 5) * 10) / 10
            : Math.round((index - 4) * (1 / 5) * 10) / 10,
          index < 4
            ? Math.round((index * (1 / 5) + 0.2) * 10) / 10
            : Math.round(((index - 4) * (1 / 5) + 0.2) * 10) / 10,
        ],
        y: [index < 4 ? 0 : 0.5, index < 4 ? 0.4 : 0.9],
      },
    })
  );
  const pieTitles = Object.entries(grouped).map(
    ([key, element], index, array) => ({
      xref: 'paper',
      yref: 'paper',
      xanchor: 'bottom',
      yanchor: 'bottom',
      showarrow: false,
      text: key,
      //   x: [
      //     index <= 4 ? index * (1 / 5) : (index - 5) * (1 / 5),
      //     index <= 4 ? index * (1 / 5) + 0.2 : (index - 5) * (1 / 5) + 0.2,
      //   ],
      //   y: [index <= 4 ? 0 : 0.6, index <= 4 ? 0.4 : 1],
      x: index < 4 ? index * (1 / 5) + 0.1 : (index - 4) * (1 / 5) + 0.1,

      y: index < 4 ? 0.4 : 0.9,
    })
  );

  console.log(pieTitles);
  console.log(tracePies);
  const traces = [...tracePies];
  //const traces = [pie0, pie1];
  console.log(traces);
  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: 500,
    autosize: true,
    grid: { rows: 2, columns: 5 },
    annotations: pieTitles,
  };
  const config = {
    //responsive: true,
  };
  return { traces: traces, layout: layout, config };
}
