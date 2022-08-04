export default function mutationalPatternBar(data) {
  console.log(data);
  var trace1 = {
    x: ['giraffes', 'orangutans', 'monkeys'],
    y: [20, 14, 23],
    name: 'SF Zoo',
    type: 'bar',
  };

  var trace2 = {
    x: ['giraffes', 'orangutans', 'monkeys'],
    y: [12, 18, 29],
    name: 'LA Zoo',
    type: 'bar',
  };

  var traces = [trace1, trace2];
  var layout = { barmode: 'group' };
  var config = {
    responsive: true,
  };

  return { traces: traces, layout: layout, config };
}
