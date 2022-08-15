import { visualizationApiSlice } from '../../../../services/store/rootApi';
import profilerSummary from '../../../controls/plotly/profilerSummary/profilerSummary';

export const profilerSummaryApiSlice = visualizationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    profilerSummary: builder.query({
      query: (params) => ({
        url: 'seqmatrixSummary',
        params,
      }),
      transformResponse: (data) => profilerSummary(data),
    }),
  }),
});

export const { useProfilerSummaryQuery } = profilerSummaryApiSlice;
