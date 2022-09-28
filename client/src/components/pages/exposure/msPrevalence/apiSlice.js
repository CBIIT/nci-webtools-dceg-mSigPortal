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
        console.log('prevalence api');
        console.log(data);
        console.log(arg);
        const { minimum } = arg;
        console.log(minimum);

        return MsPrevalence(data, minimum);
        //return MsPrevalence(transform, groupBySample);
      },
    }),
  }),
});

export const { useMsPrevelencePlotQuery } = msPrevalenceApiSlice;
