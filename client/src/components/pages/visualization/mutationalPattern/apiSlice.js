import { apiSlice } from '../../../../services/apiSlice';

export const mutationalPatternApiSlice = apiSlice.injectEndpoints({
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
