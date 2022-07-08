import { apiSlice } from '../../../../services/apiSlice';
import TMB from '../../../controls/plotly/tmb/tmb';

export const tmbApiSlice = apiSlice.injectEndpoints({
  endpoints: (builder) => ({
    tmbPlot: builder.query({
      query: (params) => ({
        url: 'explorationTmbData',
        params,
      }),
      transformResponse: (data) => {
        return TMB(data);
      },
    }),
  }),
});

export const { useTmbPlotQuery } = tmbApiSlice;
