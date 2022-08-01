import { visualizationApiSlice } from '../../../../services/store/rootApi';

export const mutationalPatternApiSlice = visualizationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    mutationalPattern: builder.query({
      query: (params) => ({
        url: 'visualizationWrapper',
        method: 'POST',
        body: params,
      }),
    }),
  }),
});

export const { useMutationalPatternQuery } = mutationalPatternApiSlice;
