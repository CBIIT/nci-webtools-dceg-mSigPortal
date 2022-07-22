import { visualizationApiSlice } from '../../../../services/store/rootApi';

export const profilerSummaryApiSlice = visualizationApiSlice.injectEndpoints({
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
        url: 'profilerSummary',
        params,
      }),
    }),
  }),
});

export const { useProfilerSummaryQuery, useProfilerSummary2Query } =
  profilerSummaryApiSlice;
