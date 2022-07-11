import { apiSlice } from '../../../../services/apiSlice';

export const profilerSummaryApiSlice = apiSlice.injectEndpoints({
  endpoints: (builder) => ({
    profilerSummary: builder.query({
      query: (params) => ({
        url: 'visualizationWrapper',
        method: 'POST',
        body: params,
      }),
    }),
    profilerSummary2: builder.query({
      query: (params) => ({
        url: 'seqmatrix',
        params,
      }),
    }),
  }),
});

export const { useProfilerSummaryQuery, useProfilerSummary2Query } =
  profilerSummaryApiSlice;
