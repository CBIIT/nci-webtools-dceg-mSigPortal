import { compareProfiles } from './profileComparison';
import { MsIndividualComparison } from '../msIndividual/msIndividual';

export default function id83(data1, data2, tab) {
  const colors = {
    '1:Del:C': '#FBBD6F',
    '1:Del:T': '#FE8002',
    '1:Ins:C': '#AEDD8A',
    '1:Ins:T': '#35A12E',
    '2:Del:R': '#FCC9B4',
    '3:Del:R': '#FB8969',
    '4:Del:R': '#F04432',
    '5:Del:R': '#BB1A1A',
    '2:Ins:R': '#CFDFF0',
    '3:Ins:R': '#93C3DE',
    '4:Ins:R': '#4B97C7',
    '5:Ins:R': '#1863AA',
    '2:Del:M': '#E1E1EE',
    '3:Del:M': '#B5B5D6',
    '4:Del:M': '#8482BC',
    '5:Del:M': '#62409A',
  };

  const mutationRegex = /^(.{7})/;
  const mutationLabels = (e) =>
    e.data.length > 3 ? e.mutation : e.mutation[0];
  const formatTickLabels = (mutationGroups) =>
    mutationGroups
      .map(({ data }) =>
        data.map((_, index) => (index >= 5 ? index + 1 + '+' : index + 1))
      )
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
