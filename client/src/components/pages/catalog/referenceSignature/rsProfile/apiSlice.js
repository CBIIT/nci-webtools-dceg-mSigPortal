import { catalogApiSlice } from '../../../../../services/store/rootApi';
import { groupBy } from 'lodash';

export const rsProfileApiSlice = catalogApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    rsProfile: builder.query({
      query: (params) => ({
        url: 'mutational_signature',
        params,
      }),
      transformResponse: (data, meta, arg) => {
        console.log(data);
        return { data };
      },
    }),
  }),
});

export const { useRsProfileQuery } = rsProfileApiSlice;
