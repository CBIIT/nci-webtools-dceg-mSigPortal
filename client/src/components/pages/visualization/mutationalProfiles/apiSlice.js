import { apiSlice } from '../../../../services/apiSlice';

export const mutationalProfilesApiSlice = apiSlice.injectEndpoints({
  endpoints: (builder) => ({
    mutationalProfiles: builder.query({
      query: (params) => ({ url: 'mutationalProfiles', params }),
      transformResponse: (data, meta, args) => {
        if (args.profile == 'SBS' || args.profile == 'DBS') {
          return baseSubstitution(data, args.profile, args.matrix);
        } else if (args.profile == 'ID') {
          return indel(data, args.matrix);
        }
      },
    }),
  }),
});

export const { useMutationalProfilesQuery } = mutationalProfilesApiSlice;

// organize seqmatrix data into a format suitable for graphing in plotly
const regexMap = {
  SBS96: /\[(.*)\]/,
  DBS78: /^(.{2})/,
  ID83: /^(.{7})/,
};

// SBS/DBS
export function baseSubstitution(data, profile, matrix) {
  const groupByMutation = data.reduce((acc, e, i) => {
    const mutationRegex = regexMap[profile + matrix];
    const mutation = e.mutationType.match(mutationRegex)[1];

    acc[mutation] = acc[mutation] ? [...acc[mutation], e] : [e];
    return acc;
  }, {});

  const transform = Object.entries(groupByMutation).map(([mutation, data]) => ({
    mutation,
    data,
  }));

  return transform;
}

// ID
export function indel(data, matrix) {
  const groupByIndel = data.reduce((acc, e, i) => {
    const indel = e.mutationType.match(regexMap['ID' + matrix])[1];

    acc[indel] = acc[indel] ? [...acc[indel], e] : [e];
    return acc;
  }, {});

  const transform = Object.entries(groupByIndel).map(([indel, data]) => ({
    indel,
    data,
  }));

  return transform;
}
