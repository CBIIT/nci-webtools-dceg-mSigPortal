import { explorationApiSlice } from '../../../../services/store/rootApi';
import TMBSignature from '../../../controls/plotly/tmb/tmb';
import { groupBy } from 'lodash';

export const tmbSignatureApiSlice = explorationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    tmbSignaturesPlot: builder.query({
      query: (params) => ({
        url: 'signature_activity',
        params,
      }),
      transformResponse: (data, meta, arg) => {
        // calculate median burden across samples
        const groupBySignature = groupBy(data, 'signatureName');
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
              sampleSize: data.length,
              signatureSize: burdens.length,
            };
          })
          .filter((e) => e.medianBurden)
          .sort((a, b) => a.medianBurden - b.medianBurden);
        return TMBSignature(transform, 'TMBSignature');
      },
    }),
  }),
});

export const { useTmbSignaturesPlotQuery } = tmbSignatureApiSlice;
