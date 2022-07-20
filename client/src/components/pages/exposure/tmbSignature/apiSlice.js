import { apiSlice } from '../../../../services/apiSlice';
import TMBSignature from '../../../controls/plotly/tmbsignature/tmbSignature';

export const tmbSignatureApiSlice = apiSlice.injectEndpoints({
  endpoints: (builder) => ({
    tmbSignaturesPlot: builder.query({
      query: (params) => ({
        url: 'tmbSignature',
        params,
      }),
      transformResponse: (data, meta, arg) => {
        return TMBSignature(data);
      },
    }),
  }),
});

export const { useTmbSignaturesPlotQuery } = tmbSignatureApiSlice;
