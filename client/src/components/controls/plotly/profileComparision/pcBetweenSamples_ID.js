import { groupBy } from 'lodash';

export default function pcBetweenSamples(rawData, args) {
  console.log(rawData);
  console.log(args);
  const samples = args.sample.split(',');

  const groupBySample = groupBy(rawData, 'sample');
  const sample1 = groupBySample[samples[0]].flat();
  const sample2 = groupBySample[samples[1]].flat();
  console.log(sample1);
  console.log(sample2);

  const colors = {
    '1:Del:C': { shape: '#FBBD6F', text: 'black' },
    '1:Del:T': { shape: '#FE8002', text: 'white' },
    '1:Ins:C': { shape: '#AEDD8A', text: 'black' },
    '1:Ins:T': { shape: '#35A12E', text: 'white' },
    '2:Del:R': { shape: '#FCC9B4', text: 'black' },
    '3:Del:R': { shape: '#FB8969', text: 'black' },
    '4:Del:R': { shape: '#F04432', text: 'black' },
    '5:Del:R': { shape: '#BB1A1A', text: 'white' },
    '2:Ins:R': { shape: '#CFDFF0', text: 'black' },
    '3:Ins:R': { shape: '#93C3DE', text: 'black' },
    '4:Ins:R': { shape: '#4B97C7', text: 'black' },
    '5:Ins:R': { shape: '#1863AA', text: 'white' },
    '2:Del:M': { shape: '#E1E1EE', text: 'blacl' },
    '3:Del:M': { shape: '#B5B5D6', text: 'black' },
    '4:Del:M': { shape: '#8482BC', text: 'black' },
    '5:Del:M': { shape: '#62409A', text: 'white' },
  };

  const groupByMutation1 = sample1.reduce((acc, e, i) => {
    const mutationRegex = /^(.{7})/;
    const mutation = e.mutationType.match(mutationRegex)[1];
    acc[mutation] = acc[mutation] ? [...acc[mutation], e] : [e];
    return acc;
  }, {});

  console.log(groupByMutation1);
  const sample1data = Object.entries(groupByMutation1).map(
    ([mutation, data]) => ({
      mutation,
      data,
    })
  );
  console.log(sample1data);

  const groupByMutation2 = sample2.reduce((acc, e, i) => {
    const mutationRegex = /^(.{7})/;
    const mutation = e.mutationType.match(mutationRegex)[1];
    acc[mutation] = acc[mutation] ? [...acc[mutation], e] : [e];
    return acc;
  }, {});
  console.log(groupByMutation2);
  const sample2data = Object.entries(groupByMutation2).map(
    ([mutation, data]) => ({
      mutation,
      data,
    })
  );
  console.log(sample2data);
  const indelNames = sample1data
    .map((indel) =>
      indel.data.map((e) => ({
        indel: indel.indel,
        index:
          indel.indel.substring(2, 5) == 'Del'
            ? +e.mutationType.slice(-1) + 1
            : e.mutationType.slice(-1),
        //index: e.mutationType.slice(-1),
      }))
    )
    .flat();
  console.log(indelNames);
}
