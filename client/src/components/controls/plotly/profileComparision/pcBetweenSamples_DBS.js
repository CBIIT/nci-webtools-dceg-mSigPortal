import { groupBy } from 'lodash';

export default function pcBetweenSamples(rawData, args) {
  console.log(rawData);
  console.log(args);
  const groupBySample = groupBy(rawData, 'sample');
  console.log(groupBySample);
  // const maxMutation = Math.max(
  //   ...rawData.map((mutation) => mutation.data.map((e) => e.mutations)).flat()
  // );
  // console.log(maxMutation);
  const colors = {
    AC: '#09BCED',
    AT: '#0266CA',
    CC: '#9FCE62',
    CG: '#006501',
    CT: '#FF9898',
    GC: '#E22925',
    TA: '#FEB065',
    TC: '#FD8000',
    TG: '#CB98FD',
    TT: '#4C0299',
  };
}
