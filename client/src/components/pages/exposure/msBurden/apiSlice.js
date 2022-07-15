import { apiSlice } from '../../../../services/apiSlice';
import MsBurden from '../../../controls/plotly//msBurden/msBurden';

export const msBurdenApiSlice = apiSlice.injectEndpoints({
  endpoints: (builder) => ({
    tmbSignaturesPlot: builder.query({
      query: (params) => ({
        url: 'tmb',
        params,
      }),
      transformResponse: (data, meta, arg) => {
        return MsBurden(data);
      },
    }),
  }),
});

export const { useTmbPlotQuery } = msBurdenApiSlice;
