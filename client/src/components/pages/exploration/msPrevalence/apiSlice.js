import { explorationApiSlice } from '../../../../services/store/rootApi';
import MsPrevalence from '../../../controls/plotly/msPrevalence/msPrevalence';

export const msPrevalenceApiSlice = explorationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    msPrevelencePlot: builder.query({
      query: ({ minimum, ...params }) => ({
        url: 'signature_activity',
        params: { ...params, limit: 1000000 },
      }),
      transformResponse: (data, meta, arg) => {
        const { minimum } = arg;
        return MsPrevalence(data, minimum);
      },
    }),
  }),
});

export const { useMsPrevelencePlotQuery } = msPrevalenceApiSlice;
