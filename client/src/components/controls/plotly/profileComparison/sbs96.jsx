import { compareProfiles } from './profileComparison';
import { MsIndividualComparison } from '../msIndividual/msIndividual';

export default function sbs96(data1, data2, tab) {
  const colors = {
    'C>A': '#03BCEE',
    'C>G': 'black',
    'C>T': '#E32926',
    'T>A': '#CAC9C9',
    'T>C': '#A1CE63',
    'T>G': '#EBC6C4',
  };

  const mutationRegex = /\[(.*)\]/;
  const mutationLabels = (e) => `<b>${e.mutation}</b>`;
  const formatTickLabels = (mutationGroups) =>
    mutationGroups
      .map(({ mutation, data }) =>
        data.map((e) => {
          const color = colors[mutation];
          const regex = /^(.)\[(.).{2}\](.)$/;
          const match = e.mutationType.match(regex);

          return `${match[1]}<span style="color:${color}"><b>${match[2]}</b></span>${match[3]}`;
        })
      )
      .flat();
  if (tab === 'msIndividual') {
    return MsIndividualComparison(
      data1,
      data2,
      colors,
      mutationRegex,
      mutationLabels,
      formatTickLabels
    );
  } else {
    return compareProfiles(
      data1,
      data2,
      colors,
      mutationRegex,
      mutationLabels,
      formatTickLabels
    );
  }
}
