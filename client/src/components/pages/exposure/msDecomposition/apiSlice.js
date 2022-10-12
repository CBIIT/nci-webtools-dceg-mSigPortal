import { explorationApiSlice } from '../../../../services/store/rootApi';
import MsDecomposition from '../../../controls/plotly/msDecomposition/msDecomposition';
import { groupBy } from 'lodash';

export const msDecompositionApiSlice = explorationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    msDecomposition: builder.query({
      query: (params) => ({
        url: 'mutational_activity',
        params,
      }),
      transformResponse: (data, meta, arg) => {
        // calculate median burden across cancer types
        const groupByCancer = groupBy(data, 'cancer');
        const transform = Object.entries(groupByCancer)
          .map(([cancer, data]) => {
            const samples = Object.values(groupBy(data, 'sample'))
              .map((e) => ({
                sample: e[0].sample,
                burden: e[0].cancerBurden,
              }))
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

            console.log(samples);

            return {
              cancer,
              samples,
              medianBurden,
              totalSamples: samples.length,
            };
          })
          .filter((e) => e.medianBurden)
          .sort((a, b) => a.medianBurden - b.medianBurden);

        return MsDecomposition(transform);
      },
    }),
  }),
});

export const { useMsDecompositionQuery } = msDecompositionApiSlice;
