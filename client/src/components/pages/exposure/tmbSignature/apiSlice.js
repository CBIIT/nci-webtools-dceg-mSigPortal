import { apiSlice } from '../../../../services/apiSlice';
import TMBSignature from '../../../controls/plotly/tmb/tmbSignature';

export const tmbApiSlice = apiSlice.injectEndpoints({
  endpoints: (builder) => ({
    tmbSignaturesPlot: builder.query({
      query: (params) => ({
        url: 'tmb',
        params,
      }),
      transformResponse: (data, meta, arg) => {
        return TMBSignature(data, arg);
      },
    }),
  }),
});

export const { useTmbPlotQuery } = tmbApiSlice;
