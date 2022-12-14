import { compareProfiles } from './profileComparison';

export default function sbs96(samples, apiData) {
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

  return compareProfiles(
    samples,
    apiData,
    colors,
    mutationRegex,
    mutationLabels,
    formatTickLabels
  );
}
