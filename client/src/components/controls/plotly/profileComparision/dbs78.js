import { compareProfiles } from './profileComparison';

export default function dbs78(samples, apiData) {
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

  return compareProfiles(
    samples,
    apiData,
    colors,
    mutationRegex,
    mutationLabels,
    formatTickLabels
  );
}
