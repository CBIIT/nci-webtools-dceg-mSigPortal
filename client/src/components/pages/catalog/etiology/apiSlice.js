import { catalogApiSlice } from '../../../../services/store/rootApi';

export const etiologyApiSlice = catalogApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    etiology: builder.query({
      query: (params) => ({
        url: 'etiology',
        params,
      }),
    }),
  }),
});

export const { useEtiologyQuery } = etiologyApiSlice;
