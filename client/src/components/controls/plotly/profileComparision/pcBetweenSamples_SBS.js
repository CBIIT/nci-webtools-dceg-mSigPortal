import { groupBy } from 'lodash';

export default function pcBetweenSamples(rawData, args) {
  console.log(rawData);
  console.log(args);
  const samples = args.sample.split(',');
  console.log(samples);
  const groupBySample = groupBy(rawData, 'sample');
  console.log(groupBySample);
  // const maxMutation = Math.max(
  //   rawData.map((mutation) => mutation.data.map((e) => e.mutations)).flat()
  // );
  // console.log(maxMutation);
  const sample1 = groupBySample[samples[0]];
  const sample2 = groupBySample[samples[1]];
  console.log(sample1);
  console.log(sample2);
  const colors = {
    'C>A': '#03BCEE',
    'C>G': 'black',
    'C>T': '#E32926',
    'T>A': '#CAC9C9',
    'T>C': '#A1CE63',
    'T>G': '#EBC6C4',
  };
}
