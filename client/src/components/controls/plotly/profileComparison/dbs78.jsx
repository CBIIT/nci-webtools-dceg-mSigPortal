import { compareProfiles } from './profileComparison';
import { MsIndividualComparison } from '../msIndividual/msIndividual';

export default function dbs78(data1, data2, tab) {
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

  const mutationRegex = /^(.{2})/;
  const mutationLabels = (e) => `<b>${e.mutation}>NN</b>`;
  const formatTickLabels = (mutationGroups) =>
    mutationGroups
      .map(({ data }) => data.map((e) => e.mutationType.slice(-2)))
      .flat();

  if (tab === 'pc') {
    return compareProfiles(
      data1,
      data2,
      colors,
      mutationRegex,
      mutationLabels,
      formatTickLabels
    );
  } else {
    return MsIndividualComparison(
      data1,
      data2,
      colors,
      mutationRegex,
      mutationLabels,
      formatTickLabels
    );
  }
}
