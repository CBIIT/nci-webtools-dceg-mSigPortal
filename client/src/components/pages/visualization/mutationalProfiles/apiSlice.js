import { apiSlice } from '../../../../services/apiSlice';
import SBS6 from '../../../controls/plotly/mutationalSignature/sbs6';
import SBS24 from '../../../controls/plotly/mutationalSignature/sbs24';
import SBS96 from '../../../controls/plotly/mutationalSignature/sbs96';
import SBS192 from '../../../controls/plotly/mutationalSignature/sbs192';
import SBS288 from '../../../controls/plotly/mutationalSignature/sbs288';
import SBS384 from '../../../controls/plotly/mutationalSignature/sbs384';
import SBS1536 from '../../../controls/plotly/mutationalSignature/sbs1536';
import DBS78 from '../../../controls/plotly/mutationalSignature/dbs78';
import DBS186 from '../../../controls/plotly/mutationalSignature/dbs186';
import ID83 from '../../../controls/plotly/mutationalSignature/id83';
import ID28 from '../../../controls/plotly/mutationalSignature/id28';
import ID415 from '../../../controls/plotly/mutationalSignature/id415';

export const mutationalProfilesApiSlice = apiSlice.injectEndpoints({
  endpoints: (builder) => ({
    mutationalProfiles: builder.query({
      query: (params) => ({ url: 'mutationalProfiles', params }),
      transformResponse: (data, meta, args) => {
        const { profile, matrix, sample } = args;
        const profileMatrix = profile + matrix;

        if (profileMatrix == 'SBS6') {
          return SBS6(data, sample);
        } else if (profileMatrix == 'SBS24') {
          return SBS24(data, sample);
        } else if (profileMatrix == 'SBS96') {
          const transform = baseSubstitution(data, profileMatrix);
          return SBS96(transform, sample);
        } else if (profileMatrix == 'SBS192') {
          return SBS192(data, sample);
        } else if (profileMatrix == 'SBS288') {
          return SBS288(data, sample);
        } else if (profileMatrix == 'SBS384') {
          return SBS384(data, sample);
        } else if (profileMatrix == 'SBS1536') {
          return SBS1536(data, sample);
        } else if (profileMatrix == 'DBS78') {
          const transform = baseSubstitution(data, profileMatrix);
          return DBS78(transform, sample);
        } else if (profileMatrix == 'DBS186') {
          return DBS186(data, sample);
        } else if (profileMatrix == 'ID83') {
          const transform = indel(data, profileMatrix);
          return ID83(transform, sample);
        } else if (profileMatrix == 'ID28') {
          return ID28(data, sample);
        } else if (profileMatrix == 'ID415') {
          return ID415(data, sample);
        } else {
          throw `Unsupported profile and matrix: ${profile} - ${matrix}`;
        }
      },
    }),
  }),
});

export const { useMutationalProfilesQuery } = mutationalProfilesApiSlice;

// organize seqmatrix data into a format suitable for graphing in plotly
const regexMap = {
  SBS24: /^.{2}(.*)$/,
  SBS96: /\[(.*)\]/,
  DBS78: /^(.{2})/,
  ID83: /^(.{7})/,
};

// SBS/DBS
export function baseSubstitution(data, profileMatrix) {
  const groupByMutation = data.reduce((acc, e, i) => {
    const mutationRegex = regexMap[profileMatrix];
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
export function indel(data, profileMatrix) {
  const groupByIndel = data.reduce((acc, e, i) => {
    const indel = e.mutationType.match(regexMap[profileMatrix])[1];

    acc[indel] = acc[indel] ? [...acc[indel], e] : [e];
    return acc;
  }, {});

  const transform = Object.entries(groupByIndel).map(([indel, data]) => ({
    indel,
    data,
  }));

  return transform;
}
