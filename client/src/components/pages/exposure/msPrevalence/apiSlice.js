import { explorationApiSlice } from '../../../../services/store/rootApi';
import MsPrevalence from '../../../controls/plotly/msPrevalence/msPrevalence';

export const msPrevalenceApiSlice = explorationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    msPrevelencePlot: builder.query({
      // query: ({ mutation, ...params }) => ({
      query: ({ minimum, ...params }) => ({
        url: 'mutational_activity',
        params,
      }),
      transformResponse: (data, meta, arg) => {
        const { minimum } = arg;
        return MsPrevalence(data, minimum);
      },
    }),
  }),
});

export const { useMsPrevelencePlotQuery } = msPrevalenceApiSlice;
