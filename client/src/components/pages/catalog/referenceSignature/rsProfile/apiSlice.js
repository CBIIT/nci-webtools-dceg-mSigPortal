import { catalogApiSlice } from '../../../../../services/store/rootApi';
import { groupBy } from 'lodash';

export const rsProfileApiSlice = catalogApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    rsProfileOptions: builder.query({
      query: (params) => ({
        url: 'mutational_signature_options',
        params,
      }),
      transformResponse: (data, meta, arg) => {
        console.log(data);
        return { data };
      },
    }),
    rsProfileData: builder.query({
      query: (params) => ({
        url: 'mutational_signature',
        params,
      }),
      transformResponse: (data, meta, args) => {
        const { profile, matrix } = args;
        const profileMatrix = profile + matrix;
        console.log(profileMatrix);
        console.log(data);
        return { data };
      },
    }),
  }),
});

export const { useRsProfileOptionsQuery, useRsProfileDataQuery } =
  rsProfileApiSlice;
