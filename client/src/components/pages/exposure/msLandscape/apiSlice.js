import { explorationApiSlice } from '../../../../services/store/rootApi';
import MsLandscape from '../../../controls/plotly/msLandscape/msLandscape';

export const msLandscapeApiSlice = explorationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    msLandscapePlot: builder.query({
      query: (params) => ({
        url: 'signature_landscape',
        params,
      }),
      transformResponse: (data) => {
        const { cosineData, exposureData } = data.output;
        console.log(cosineData);
        return MsLandscape(cosineData, exposureData);
      },
    }),

    // msLandscapePlot: builder.query({
    //   async queryFn(_arg, _queryApi, _extraOptions, fetchWithBQ) {
    //     try {
    //       const res = await Promise.all([
    //         fetchWithBQ(
    //           '/signature_activity?' +
    //             new URLSearchParams(_arg.params_activity)
    //         ), //exposure
    //         fetchWithBQ(
    //           '/mutational_signature?' +
    //             new URLSearchParams(_arg.params_signature)
    //         ), //signature
    //         fetchWithBQ(
    //           '/mutational_spectrum?' +
    //             new URLSearchParams(_arg.params_spectrum)
    //         ), //seqmatrix
    //       ]);

    //       console.log(res);
    //       console.log(_arg);
    //       //return MsLandscape(res, _arg);
    //       return { data: MsLandscape(res, _arg) };
    //     } catch (error) {
    //       return { error };
    //     }
    //   },
    // }),
  }),
});

export const { useMsLandscapePlotQuery, useMsLandscapePlot2Query } =
  msLandscapeApiSlice;
