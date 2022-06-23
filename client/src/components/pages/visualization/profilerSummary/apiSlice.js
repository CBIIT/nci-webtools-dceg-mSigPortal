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
  }),
});

export const { useProfilerSummaryQuery } = profilerSummaryApiSlice;
