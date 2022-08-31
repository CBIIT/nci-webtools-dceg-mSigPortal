import { catalogApiSlice } from '../../../../services/store/rootApi';

export const etiologyApiSlice = catalogApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    etiologyOptions: builder.query({
      query: (params) => ({
        url: 'etiologyOptions',
        params,
      }),
    }),
  }),
});

export const { useEtiologyOptionsQuery } = etiologyApiSlice;
