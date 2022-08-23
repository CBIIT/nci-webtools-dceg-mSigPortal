import { explorationApiSlice } from '../../../../services/store/rootApi';
import MsBurden from '../../../controls/plotly/tmb/tmb';
import { groupBy } from 'lodash';

export const msBurdenApiSlice = explorationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    msBurden: builder.query({
      query: (params) => ({
        url: 'mutational_activity',
        params,
      }),
      transformResponse: (data, meta, arg) => {
        const { signatureName } = arg;
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

            return {
              cancer,
              samples,
              medianBurden,
              totalSamples: samples.length,
            };
          })
          .filter((e) => e.medianBurden)
          .sort((a, b) => a.medianBurden - b.medianBurden);

        return MsBurden(transform, 'MSBurden', signatureName);
      },
    }),
  }),
});

export const { useMsBurdenQuery } = msBurdenApiSlice;
