import { explorationApiSlice } from '../../../../services/store/rootApi';
import MsLandscape from '../../../controls/plotly/msLandscape/msLandscape';
import { groupBy } from 'lodash';

export const msLandscapeApiSlice = explorationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    // msLandscapePlot: builder.query({
    //   query: (params) => ({
    //     url: 'mutational_activity',
    //     params,
    //   }),
    //   transformResponse: (data, meta, arg) => {
    //     console.log(data);
    //     // calculate median burden across cancer types
    //     const groupByCancer = groupBy(data, 'cancer');
    //     const transform = Object.entries(groupByCancer)
    //       .map(([cancer, data]) => {
    //         const samples = Object.values(groupBy(data, 'sample'))
    //           .map((e) => ({
    //             sample: e[0].sample,
    //             burden: e[0].cancerBurden,
    //           }))
    //           .sort((a, b) => a.burden - b.burden);

    //         const burdens = samples
    //           .filter((e) => e.burden)
    //           .map((e) => e.burden);

    //         const medianBurden =
    //           burdens.length % 2 == 0
    //             ? (burdens[burdens.length / 2] +
    //                 burdens[burdens.length / 2 - 1]) /
    //               2
    //             : burdens[Math.floor(burdens.length / 2)];

    //         return {
    //           cancer,
    //           samples,
    //           medianBurden,
    //           totalSamples: samples.length,
    //         };
    //       })
    //       .filter((e) => e.medianBurden)
    //       .sort((a, b) => a.medianBurden - b.medianBurden);

    //     return MsLandscape(transform, 'MsLandscape');
    //   },
    // }),
    msLandscapePlot: builder.query({
      async queryFn(_arg, _queryApi, _extraOptions, fetchWithBQ) {
        try {
          const res = await Promise.all([
            fetchWithBQ(
              '/mutational_activity?' +
                new URLSearchParams(_arg.params_activity)
            ), //exposure
            fetchWithBQ(
              '/mutational_signature?' +
                new URLSearchParams(_arg.params_signature)
            ), //signature
            fetchWithBQ(
              '/mutational_spectrum?' +
                new URLSearchParams(_arg.params_spectrum)
            ), //seqmatrix
          ]);

          console.log(res);
          console.log(_arg);
          //return MsLandscape(res, _arg);
          return { data: MsLandscape(res, _arg) };
        } catch (error) {
          return { error };
        }
      },
    }),
  }),
});

export const { useMsLandscapePlotQuery } = msLandscapeApiSlice;
