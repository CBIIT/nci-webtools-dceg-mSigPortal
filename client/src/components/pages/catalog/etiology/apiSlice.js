import { catalogApiSlice } from '../../../../services/store/rootApi';
import { groupBy } from 'lodash';
import TMB from '../../../controls/plotly/tmb/tmb';
import SBS96 from '../../../controls/plotly/mutationalProfiles/sbs96';
import DBS78 from '../../../controls/plotly/mutationalProfiles/dbs78';
import ID83 from '../../../controls/plotly/mutationalProfiles/id83';

export const etiologyApiSlice = catalogApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    etiologyOptions: builder.query({
      query: (params) => ({
        url: 'signatureEtiologyOptions',
        params,
      }),
    }),
    etiologyData: builder.query({
      query: (params) => ({
        url: 'signatureEtiology',
        params,
      }),
      transformResponse: (data, _, arg) => {
        const groupByCancer = groupBy(data, 'cancer');
        const transform = Object.entries(groupByCancer)
          .map(([cancer, data]) => {
            const samples = Object.values(groupBy(data, 'sample'))
              .map((e) => ({
                sample: e[0].sample,
                burden: e[0].burden,
                mutations: e[0].mutations,
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

        return {
          distributionPlot: TMB(transform, '', arg.signatureName),
          data: transform,
        };
      },
    }),
    etiologyProfile: builder.query({
      query: (params) => ({
        url: 'signatureEtiology',
        params,
      }),
      transformResponse: (data, _, arg) => {
        if (arg.signatureName.includes(/sbs/gi)) return SBS96(data);
        else if (arg.signatureName.includes(/dbs/gi)) return DBS78(data);
        else if (arg.signatureName.includes(/id/gi)) return ID83(data);
        else throw 'No supported profiles for ' + arg.signatureName;
      },
    }),
    thumbnails: builder.query({
      query: (body) => ({
        url: 'getImageS3Batch',
        method: 'POST',
        body,
      }),
    }),
  }),
});

export const {
  useEtiologyOptionsQuery,
  useEtiologyDataQuery,
  useEtiologyProfileQuery,
  useThumbnailsQuery,
} = etiologyApiSlice;
