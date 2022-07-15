import { apiSlice } from '../../../../services/apiSlice';
import TMBSignature from '../../../controls/plotly/tmbsignature/tmbSignature';

export const tmbApiSlice = apiSlice.injectEndpoints({
  endpoints: (builder) => ({
    tmbSignaturesPlot: builder.query({
      query: (cancer) => ({
        url: 'tmb',
        params: cancer,
      }),
      transformResponse: (data, meta, arg) => {
        return TMBSignature(data);
      },
    }),
  }),
});

export const { useTmbPlotQuery } = tmbApiSlice;
