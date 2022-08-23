import { explorationApiSlice } from '../../../../services/store/rootApi';
import MsPrevalence from '../../../controls/plotly/msPrevalence/msPrevalence';
import { groupBy } from 'lodash';

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
        const { minimum } = arg;
        console.log(minimum);

        // calculate median burden across cancer types
        const groupBySignature = groupBy(data, 'signatureName');
        const groupBySample = groupBy(data, 'sample');

        const transform = Object.entries(groupBySignature)
          .map(([signatureName, data]) => {
            const samples = data
              .filter((e) => e.exposure)
              .sort((a, b) => a.burden - b.burden);

            const burdens = samples
              .filter((e) => e.burden)
              .map((e) => e.burden);

            const medianBurden =
              burdens.length % 2 == 0
                ? (burdens[burdens.length / 2] +
                    burdens[burdens.length / 2 - 1]) /
                  2
                : burdens[Math.floor(burdens.length / 2)];

            return {
              signatureName,
              samples,
              medianBurden,
              totalSamples: data.length,
            };
          })
          .filter((e) => e.medianBurden)
          .sort((a, b) => a.medianBurden - b.medianBurden);

        return MsPrevalence(transform, minimum);
        //return MsPrevalence(transform, groupBySample);
      },
    }),
  }),
});

export const { useMsPrevelencePlotQuery } = msPrevalenceApiSlice;
